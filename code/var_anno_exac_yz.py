# Python 3 
# coding: utf-8


import requests, json
import pandas as pd
import subprocess

import sys

args = sys.argv
if(len(args)<3):
    print("Missing input...")
    print("Format: python var_anno_exac_yz.py <vcf_file> <out_file>")
    print("Example: python var_anno_exac_yz.py ../data/Challenge_data_clean.vcf ../data/Challenge_data_clean.vcf.annotated.txt")
    sys.exit()

[vcf_file, out_file] = sys.argv[1:1+2]

#vcf_file = '../data/Challenge_data_clean.vcf'

tmp_DIR = './TMP/'

gencode_file = '../data/gencode.v19.annotation.gtf'

tmp_prefix = "out.tmp"

print("Input vcf: ")
print(vcf_file)

print("TMP DIR:")
print(tmp_DIR)

print("Output table:")
print(out_file)

## Intermediate files: 

vcf_anno_tmp_file = tmp_DIR + tmp_prefix +'.INFO' ## vcf INFO file 
exac_tmp_file = tmp_DIR + tmp_prefix +'.json' ## exac query result 
vcf_tmp_bed_file = tmp_DIR + tmp_prefix+'.input.bed' ## vcf bed file 


def exac_api_query(query_str):
    """
    Bulk query of ExAC API for variants
    """
    var = query_str
    URL="http://exac.hms.harvard.edu/rest/bulk/variant/variant"
    response = requests.post(URL, data = var)
    if(response.status_code!=200):
        print('ERROR from the request')
        raise KeyboardInterrupt
    response_json_dict = response.json()
    with open( exac_tmp_file , "w") as wf:
        json.dump(response_json_dict, wf, indent=4, separators=(',', ': '))

    exacpd = pd.DataFrame.from_dict(response_json_dict, orient = 'index')
    return(exacpd)

def define_del_order():
    """
    define an ordered list of deleterious snp types
    """
    print("Using this file for ranking variation deleterious effect:")
    print('../data/VEP_variant_function_scores_Koscielny17NAR.cleaned.txt')
    print("")
    vep_var_scores_pd = pd.read_csv('../data/VEP_variant_function_scores_Koscielny17NAR.cleaned.txt', sep='\t')
    my_order = list(vep_var_scores_pd['label'].values)
    my_order_map = {key: i for i, key in enumerate(my_order)}
    return(my_order_map)

def ordered_del_type(snp_type, my_order_map):
    """
    return a ordered list of deleterious snp type
    """
    a = [t for t in snp_type if t not in my_order_map.keys()]
    if(len(a)>0):
        print("Warning: var type not found: ", a)
    return(sorted(snp_type, key=lambda d: my_order_map[d]))


def most_del_type(snpline, my_order_map ):
    """
    pick the most deleterious snp annotation for a exac line."""

    snp_type = {}
    # collect all annotations and return the lowest rank one
    for anno in snpline['vep_annotations']:
        snp_type = set(anno['Consequence'].split("&")).union(set(snp_type))
    snp_type = list(set(snp_type))
    ordered_snp = ordered_del_type(snp_type, my_order_map = my_order_map)
    if(len(ordered_snp)==0):
        print("No ordered variant type obtained ",snp_type)
        return("")
    else:
        return(ordered_snp[0])
    
def in_exac_gene(snpline):
    """
    record the gene annotation from exac. 
    """
    gene = {}
    for anno in snpline['vep_annotations']:
        gene = set([anno['SYMBOL']]).union(set(gene))
    gene = ";".join(list(gene))
    return(gene)


# ## Get INFO from vcf



print("Getting INFO from vcf : DP, SRF, SRR, SAF, SAR")
print("")
GET_INFO_OPTION =" --get-INFO DP --get-INFO SRF --get-INFO SRR --get-INFO SAF --get-INFO SAR --max-alleles 2"
#subprocess.getoutput(['vcftools --vcf '+vcf_file + GET_INFO_OPTION + ' --chr 19 --from-bp 1000000 --to-bp 19000000'])
subprocess.getoutput(['vcftools --vcf '+vcf_file + GET_INFO_OPTION + ' --out '+ tmp_DIR + tmp_prefix ]) 




vcf = pd.read_csv( vcf_anno_tmp_file , sep='\t')
print("Number of variants: ", len(vcf))
print("")


# ## Querying ExAC API



vcf['Add_snpname'] = vcf['CHROM'].astype(str) + '-' + vcf['POS'].astype(str) + '-' + vcf['REF'] + '-' + vcf['ALT']
query_str = "["+ "\"" + "\",\"".join(list(vcf['Add_snpname'].values)) + "\"" + "]"

print("Querying ExAC API...")

print("")


exacpd = exac_api_query(query_str)




print("Parsing ExAC results...")
print("")
mask_has_vep = ( exacpd['vep_annotations'].str.len()>0 ) 
my_order_map = define_del_order()
exacpd.loc[ mask_has_vep , 'var_type'] = exacpd[ mask_has_vep ][['vep_annotations']].apply(most_del_type, axis=1 , my_order_map = my_order_map)
exacpd.loc[ mask_has_vep , 'gene'] = exacpd[ mask_has_vep ][['vep_annotations']].apply(in_exac_gene, axis=1)


# ## Gene annotation



## computing ratio
vcf['ratio_A_R'] = (vcf['SRF']+vcf['SRR'])/(vcf['SAF']+vcf['SAR'])
vcf['ratio_A_R'] = vcf['ratio_A_R'].map("{:.2f}".format)

# annotate intergenic variants
print("Annotating variants using the gene annotation file...")
print("")
vcf['POS0'] = vcf['POS'] - vcf['REF'].str.len()
vcf['chr'] = "chr" + vcf['CHROM']
vcf[['chr','POS0','POS','Add_snpname']].to_csv( vcf_tmp_bed_file ,sep='\t',header=False, index=False)
intergenic = subprocess.getoutput(['bedtools intersect -v -b ' + gencode_file + ' -a '+ vcf_tmp_bed_file])
intergenic_var_list = [t.split('\t')[3] for t in intergenic.split('\n')]
exacpd.loc[intergenic_var_list]['var_type'] = 'intergenic'


# ## Merge vcf and ExAC info



## Merge exAc and vcf information

exacpd_keep_col_list = ['Add_snpname','var_type','gene']
vcf_annotated = pd.merge(vcf, exacpd.reset_index().rename(columns={"index":"Add_snpname"})[exacpd_keep_col_list], on='Add_snpname',how='left')

annotate_col_list = ['CHROM','POS','REF','ALT','var_type','DP','n_reads_var','ratio_A_R','gene']
vcf_annotated['n_reads_var'] = vcf_annotated[['SRF','SRR','SAF','SAR']].sum(axis=1)
vcf_annotated[annotate_col_list].to_csv(out_file, sep='\t',index=False)

print("DONE! Have a nice day! ")
print("Check out: "+ out_file)
