import pandas as pd
from SampleDataClean import dtf_clean

dtf_exp=pd.read_csv("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",sep="\t",compression="gzip")

# Gènes en colonnes, échantillons en lignes
dtf_exp = dtf_exp.set_index("sample").T
dtf_exp.index.name = "sample"

# Conserver uniquement les échantillons communs
dtf_final = dtf_clean.merge(dtf_exp,
    left_on="sample",
    right_index=True,
    how="inner",
    validate="one_to_one")

dtf_numeric = dtf_final[dtf_exp.columns]

#print("Dataset final :", dtf_final.shape)
#print("Dataset numérique :", dtf_numeric.shape)

#print("Types de données :")
#print(dtf_numeric.dtypes.value_counts())

#print("Valeurs non numériques :",dtf_numeric.select_dtypes(exclude="number").shape[1])

#print("Nombre total de NaN :",dtf_numeric.isna().sum().sum())

#print("Échantillons communs :", len(dtf_final))

# Vérifier les gènes
#print("Gènes en double :", dtf_numeric.columns.duplicated().sum())

# Vérifier les types d'échantillons
#print(dtf_final["sample_type"].value_counts(dropna=False))

#doublons = dtf_numeric.columns[dtf_numeric.columns.duplicated(keep=False)]
#print(doublons.tolist())

# Repérer les 4 colonnes
positions = [i for i, col in enumerate(dtf_final.columns) if col == "SLC35E2"]

# Supprimer les colonnes 3 et 4 (copies)
dtf_final = dtf_final.iloc[:, [
    i for i in range(dtf_final.shape[1])
    if i not in positions[2:]
]]

# Renommer les deux colonnes conservées
colonnes = list(dtf_final.columns)
colonnes[positions[0]] = "SLC35E2_1"
colonnes[positions[1]] = "SLC35E2_2"

dtf_final.columns = colonnes

# Actualiser la matrice numérique
dtf_numeric = dtf_final[dtf_exp.columns.drop("SLC35E2") .tolist() + ["SLC35E2_1", "SLC35E2_2"]]

#print("Final :", dtf_final.shape)
#print("Numérique :", dtf_numeric.shape)
#print("Doublons :", dtf_numeric.columns.duplicated().sum())
