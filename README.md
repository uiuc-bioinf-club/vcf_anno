# VCF annotation with ExAC API. 

A VCF annotation tool (prototype) using [vcftools](https://vcftools.github.io/index.html), [ExAC REST API](http://exac.hms.harvard.edu/), and [gencode](https://www.gencodegenes.org/releases/19.html). 
The output file is a tab-separated text file with the information following: 
1. Type of variation 
2. Depth of sequence coverage 
3. Number of reads supporting the variant 
4. Ratio of reads supporting alternative vs reference allele 
5. gene inforamtion from ExAC.

## Quick start
*Please first download gencode annotation from [here](https://www.gencodegenes.org/releases/19.html) and uncompress it as: data/gencode.v19.annotation.gtf
*The following code runs this VCF annotation tool on a dataset in the repo.

    cd vcf_anno/code/    
	python var_anno_exac_yz.py ../data/Challenge_data_clean.vcf ../data/Challenge_data_clean.vcf.annotated.txt

## How to use

### Dependencies

1. Python3 and modules: pandas, requests, json.
2. vcftools, bedtools.

### The main variant annotation tool:
1. The main script: code/var_anno_exac_yz.py
2. Supporting data: 
	* A file ranking the deleterious effect of variants: 
		data/VEP_variant_function_scores_Koscielny17NAR.cleaned.txt
	* Gencode gene annotation: 
        Gencode v19 should be downloaded from [here](https://www.gencodegenes.org/releases/19.html) and uncompressed to:
		data/gencode.v19.annotation.gtf
3. Input data: 
	* A sample VCF file: 
		data/Challenge_data_clean.vcf
4. Output result:
	* The annotated variants after running the main script on the Input data.
		data/Challenge_data_clean.vcf.annotated.txt

### A small tool/resource for ranking deleterious effect of variants
	* A table for ranking deleterious effect of variants was found in [Koscielny et al., 2017, NAR, Supplementary Table 2.](https://academic.oup.com/nar/article/45/D1/D985/2605745#51199338)
	* I have copied the table here: data/VEP_variant_function_scores_Koscielny17NAR.cleaned.txt 
	* It is then cleaned it up using this script: lib/vep_variant_order_clean.py.
	* The cleaned table generated is here, good to use: data/VEP_variant_function_scores_Koscielny17NAR.cleaned.txt 

### A small tool to cleanup vcf header 
	* Some versions of vcftools will give Warnings for commas in header lines, while some versions accept them. A small tool here is to generate a new VCF header file to remove those comma.
	* Running lib/vcf_header_clean.py on the original data will generate the new VCF header file lib/newheader.tmp.txt, and I have used bcftools to reheader the original file ../data/Challenge_data\ \(1\).vcf to make the Warning-free input: data/Challenge_data_clean.vcf.

This is an edit 
