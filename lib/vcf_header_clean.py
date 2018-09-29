
import re

def vcfheader_clean(IN_vcf_file, OUT_vcf_header_file):
    """
    Check whether a vcf file contains description with comma in them, which may raise warnings for some vcftools version such as vcftools-0.1.16
    """
    #IN_vcf_file = vcf_file
    #OUT_vcf_header_file = './header_commaclear.tmp.vcfheader'
    
    FIX_COMMA = 0
    new_header_list = []
    with open(IN_vcf_file, 'r') as vcff:
        for line in vcff:
            if(line[0]!="#"):
                break
            if(line[0:2] == '##' and line.split("=")[0] in ["##INFO", "##FORMAT"]):
                ## check comma in the description part of INFO.
                description_after = re.search(r'Description\=.+\>', line)
                if(description_after):
                    if("," in description_after.group()):
                        FIX_COMMA = 1
                    description_after_commaclear = description_after.group().replace(",", " ").strip()
                    description_before = re.search(r'.+Description', line).group().replace("Description","")
                    line_comma_clear = description_before + description_after_commaclear
                    line = line_comma_clear
            new_header_list.append(line.strip())
    if(FIX_COMMA):
        print("Fixed vcf file header is in : \n"+ OUT_vcf_header_file)
        print("It can be used to reheader the vcf file. Example: ")
        print("bcftools reheader -h "+ OUT_vcf_header_file + " --output ./newheader.tmp.txt "+ IN_vcf_file + "\n")
        with open(OUT_vcf_header_file,'w') as OUTf:
            OUTf.write("\n".join(new_header_list))
    else:
        print("The vcf file is comma-free in INFO tags.")


if __name__ == "__main__":
    """
        Generate new comma-free headers for INFO and FORMAT tags in vcffile. 
        Example: python --vcf 
    """
    import argparse
    parser = argparse.ArgumentParser(description = "clean vcf header")
    parser.add_argument('--vcf', required=True, help="vcffile")
    parser.add_argument('--newheader', required = False, default = "./newheader.tmp.txt" , help="where to put the new header")
    args = parser.parse_args()
    vcfheader_clean( args.vcf, args.newheader )
    
