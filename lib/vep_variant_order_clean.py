import sys
print("### This is for cleaning up the variant function ranking file. Already ran.")
print("Check here: ../data/VEP_variant_function_scores_Koscielny17NAR.cleaned.txt")
sys.exit()
import pandas as pd
### This is for cleaning up the variant function ranking file. Already ran.
vep_var_scores_pd = pd.read_csv('../data/VEP_variant_function_scores_Koscielny17NAR.txt', sep='\t', header=0)
vep_var_scores_pd['label'] = vep_var_scores_pd['label'].str.replace(' ','_')
vep_var_scores_pd = vep_var_scores_pd.sort_values('functional_score', ascending=False).drop_duplicates(subset = 'label', keep = 'first')
#vep_var_scores_pd
vep_var_scores_pd.to_csv('../data/VEP_variant_function_scores_Koscielny17NAR.cleaned.txt', sep='\t',index=False)
