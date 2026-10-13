import pandas as pd
import numpy as np
from DataMerge import dtf_final, dtf_numeric

# Dimensions initiales
#print("Avant :", dtf_numeric.shape)
# Variance de chaque gène (colonnes)
variances = dtf_numeric.var(axis=0)

# Identifier les gènes à variance nulle
genes_constants = variances[variances == 0].index
#print("Gènes à variance nulle :", len(genes_constants))

# Supprimer ces gènes
dtf_numeric = dtf_numeric.drop(columns=genes_constants)

# Vérification
#print("Après :", dtf_numeric.shape)
#print("NaN restants :", dtf_numeric.isna().sum().sum())

# Nombre de NaN pour chaque échantillon
na_samples = dtf_numeric.isna().sum(axis=1)
#print(na_samples.describe(percentiles=[0.50, 0.75, 0.85, 0.90, 0.95, 0.99]))

# Les 20 échantillons avec le plus de NaN
#print(na_samples.nlargest(20))

# dtf_numeric : échantillons en lignes, gènes en colonnes
# dtf_final   : mêmes échantillons + phénotypes

# 1. Matrice des valeurs manquantes
masque_na = dtf_numeric.isna().to_numpy()

# Nombre de NaN par échantillon
nb_na = masque_na.sum(axis=1)

# 2. Compresser les profils pour économiser la RAM
profils_compactes = np.packbits(masque_na, axis=1)
del masque_na

# 3. Identifier les profils strictement identiques
_, profil_id, effectifs = np.unique(
    profils_compactes,
    axis=0,
    return_inverse=True,
    return_counts=True
)

del profils_compactes

# 4. Construire le tableau de diagnostic
diagnostic = pd.DataFrame({
    "sample": dtf_final.loc[dtf_numeric.index, "sample"].to_numpy(),
    "profil": profil_id,
    "nb_na": nb_na,
    "sample_type": dtf_final.loc[dtf_numeric.index, "sample_type"].to_numpy(),
    "disease": dtf_final.loc[dtf_numeric.index, "_primary_disease"].to_numpy()
})

# 5. Résumer chaque groupe de profils
#resume = diagnostic.groupby("profil").agg(nb_echantillons=("sample", "size"),na_par_echantillon=("nb_na", "first"))

#resume["na_total"] = (resume["nb_echantillons"] * resume["na_par_echantillon"])

#resume["part_na_pct"] = (resume["na_total"] / diagnostic["nb_na"].sum() * 100)

#resume["part_echantillons_pct"] = (resume["nb_echantillons"] / len(diagnostic) * 100)

# Trier selon le nombre total de NaN concentrés
#resume = resume.sort_values("na_total", ascending=False)

#print("Nombre de profils distincts :", len(resume))
#print("Nombre total de NaN :", diagnostic["nb_na"].sum())

#print("\nTop 15 des profils concentrant le plus de NaN :")
#print(resume.head(15).round(2).to_string())

# 6. Étudier les phénotypes des 5 profils les plus problématiques
#for profil in resume.head(5).index:

    #groupe = diagnostic[diagnostic["profil"] == profil]

    #print(f"\n===== PROFIL {profil} =====")
    #print("Échantillons :", len(groupe))
    #print("NaN par échantillon :", groupe["nb_na"].iloc[0])

    #print("\nTypes d'échantillons :")
    #print(groupe["sample_type"].value_counts(dropna=False).to_string())

    #print("\nTypes de cancers :")
    #print(groupe["disease"].value_counts(dropna=False).head(15).to_string())


total = diagnostic["disease"].value_counts()

sans_na = diagnostic.loc[
    diagnostic["nb_na"] == 0, "disease"
].value_counts()

bilan = pd.DataFrame({
    "total": total,
    "conserves": sans_na
}).fillna(0).astype(int)

bilan["perdus"] = bilan["total"] - bilan["conserves"]
bilan["perte_pct"] = 100 * bilan["perdus"] / bilan["total"]

#print(bilan.sort_values("perte_pct", ascending=False).round(2).to_string())

genes_manquants = {}

for profil in sorted(diagnostic["profil"].unique()):
    # Position du premier échantillon de ce profil
    position = np.flatnonzero(
        diagnostic["profil"].to_numpy() == profil
    )[0]

    echantillon = dtf_numeric.iloc[position]

    genes_manquants[profil] = set(
        dtf_numeric.columns[echantillon.isna()]
    )

# Nombre de gènes manquants communs à chaque paire de profils
profils = [p for p in genes_manquants if genes_manquants[p]]

chevauchements = pd.DataFrame(0, index=profils, columns=profils)

for p1 in profils:
    for p2 in profils:
        chevauchements.loc[p1, p2] = len(
            genes_manquants[p1] & genes_manquants[p2]
        )

#print("Gènes manquants communs entre profils :")
#print(chevauchements)

# Ensemble de tous les gènes concernés, sans compter les doublons
ensemble_manquants = set().union(*genes_manquants.values())

#print("\nGènes avec au moins un NaN :", len(ensemble_manquants))
#print("Gènes sans aucun NaN :", dtf_numeric.shape[1] - len(ensemble_manquants))
#print("Pourcentage de gènes concernés :",round(100 * len(ensemble_manquants) / dtf_numeric.shape[1], 2))

# Gènes manquants exclusivement dans un profil
#print("\nGènes manquants uniquement dans chaque profil :")

#for profil in profils:
 #   autres = set().union(*[genes_manquants[p] for p in profils if p != profil])

    #exclusifs = genes_manquants[profil] - autres
    #print(f"Profil {profil} : {len(exclusifs)}")



# Conserver uniquement les gènes sans NaN
dtf_complet = dtf_numeric.dropna(axis=1).copy()

# Réassocier les identifiants et les phénotypes
phenotypes = dtf_final.loc[
    dtf_complet.index,
    ["sample", "sample_type", "_primary_disease"]
]

dtf_modele = pd.concat([phenotypes, dtf_complet], axis=1)

#print("Expression :", dtf_complet.shape)
#print("Avec phénotypes :", dtf_modele.shape)
#print("NaN dans l'expression :", dtf_complet.isna().sum().sum())