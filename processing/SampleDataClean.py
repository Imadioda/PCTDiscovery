import pandas as pd

df_sampl= pd.read_csv("TCGA_phenotype_denseDataOnlyDownload.tsv", sep="\t")
#print(df_sampl.columns)
#print(df_sampl["sample_type"].value_counts())
#print(df_sampl["_primary_disease"].value_counts())

#print(df_sampl.isna().sum())
#print("======================")
#print(df_sampl["_primary_disease"].isna().sum())
#print(df_sampl["sample_type"].isna().sum())
#print("Duplicated samples:", df_sampl["sample"].duplicated().sum())

dtf_ssna=df_sampl.dropna()
#print("==========verification============")
#print(dtf_ssna.isna().sum())

#print("======================")
dtf_clean=dtf_ssna.drop(columns=["sample_type_id"])
#print(dtf_clean.shape)
