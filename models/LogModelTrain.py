import joblib
import numpy as np
import pandas as pd
from ExpressionDataClean import dtf_modele

from pathlib import Path
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    balanced_accuracy_score,
    roc_auc_score,
)


def entrainer_binaire(dtf_modele, dossier="resultats_binaire"):
    dossier = Path(dossier)
    dossier.mkdir(parents=True, exist_ok=True)

    # 1. Sélection des prélèvements et définition de la cible
    df = dtf_modele.loc[
        dtf_modele["sample_type"].isin(
            ["Solid Tissue Normal", "Primary Tumor"]
        )
    ].copy()

    y = df["sample_type"].map({
        "Solid Tissue Normal": 0,
        "Primary Tumor": 1,
    })

    X = df.drop(columns=["sample", "sample_type", "_primary_disease"])

    if not X.columns.is_unique:
        raise ValueError("Les noms des gènes ne sont pas uniques.")

    if not np.isfinite(X.to_numpy()).all():
        raise ValueError("L'expression contient des NaN ou des infinis.")

    # Identifiant patient TCGA : les trois premiers blocs du barcode.
    samples = df["sample"].astype("string")
    format_valide = samples.str.match(
        r"^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-"
    )

    if not format_valide.fillna(False).all():
        raise ValueError("Certains identifiants ne sont pas au format TCGA attendu.")

    patients = samples.str.split("-").str[:3].str.join("-")

    # 2. Réserver environ 20 % pour le test final
    separation_test = StratifiedGroupKFold(
        n_splits=5, shuffle=True, random_state=42
    )

    idx_dev, idx_test = next(
        separation_test.split(X, y, groups=patients)
    )

    # Parmi les 80 % restants : environ 60 % train et 20 % validation
    separation_val = StratifiedGroupKFold(
        n_splits=4, shuffle=True, random_state=43
    )

    train_local, val_local = next(
        separation_val.split(
            X.iloc[idx_dev],
            y.iloc[idx_dev],
            groups=patients.iloc[idx_dev],
        )
    )

    idx_train = idx_dev[train_local]
    idx_val = idx_dev[val_local]

    groupes_train = set(patients.iloc[idx_train])
    groupes_val = set(patients.iloc[idx_val])
    groupes_test = set(patients.iloc[idx_test])

    assert groupes_train.isdisjoint(groupes_val)
    assert groupes_train.isdisjoint(groupes_test)
    assert groupes_val.isdisjoint(groupes_test)

    partitions = {
        "train": idx_train,
        "validation": idx_val,
        "test": idx_test,
    }

    for nom, indices in partitions.items():
        if y.iloc[indices].nunique() != 2:
            raise ValueError(f"Une classe manque dans la partition {nom}.")

        print(
            f"{nom} : {len(indices)} échantillons | "
            f"normaux={(y.iloc[indices] == 0).sum()} | "
            f"tumoraux={(y.iloc[indices] == 1).sum()}"
        )

    # Enregistrer les partitions pour les réutiliser avec le modèle 2
    repartition = df[
        ["sample", "sample_type", "_primary_disease"]
    ].copy()

    repartition["patient"] = patients
    repartition["partition"] = ""

    for nom, indices in partitions.items():
        repartition.iloc[
            indices, repartition.columns.get_loc("partition")
        ] = nom

    repartition.to_csv(dossier / "partitions.csv", index=False)

    X_train, y_train = X.iloc[idx_train], y.iloc[idx_train]
    X_val, y_val = X.iloc[idx_val], y.iloc[idx_val]

    # 3. Comparer trois forces de régularisation
    meilleur_score = -np.inf
    meilleur_modele = None
    meilleur_C = None

    for C in [0.01, 0.1, 1.0]:
        modele = Pipeline([
            ("variance", VarianceThreshold(threshold=0)),
            ("standardisation", StandardScaler()),
            ("classification", LogisticRegression(
                C=C,
                solver="lbfgs",
                class_weight="balanced",
                max_iter=3000,
            )),
        ])

        modele.fit(X_train, y_train)

        predictions = modele.predict(X_val)
        score = balanced_accuracy_score(y_val, predictions)

        print(f"C={C} | balanced accuracy validation={score:.4f}")

        if score > meilleur_score:
            meilleur_score = score
            meilleur_modele = modele
            meilleur_C = C

    # 4. Évaluation sur la validation ; test final encore intact
    predictions = meilleur_modele.predict(X_val)
    probabilites = meilleur_modele.predict_proba(X_val)[:, 1]

    print("\nMeilleur C :", meilleur_C)
    print(classification_report(
        y_val,
        predictions,
        labels=[0, 1],
        target_names=["Normal", "Tumoral"],
        digits=3,
        zero_division=0,
    ))

    print("Matrice de confusion : lignes réelles, colonnes prédites")
    print(confusion_matrix(y_val, predictions, labels=[0, 1]))
    print("ROC-AUC validation :", round(roc_auc_score(y_val, probabilites), 4))

    # 5. Sauvegarde du pipeline et des informations nécessaires
    joblib.dump({
        "pipeline": meilleur_modele,
        "genes": X.columns.tolist(),
        "classes": {0: "Normal", 1: "Tumoral"},
        "seuil": 0.5,
        "C": meilleur_C,
    }, dossier / "modele_binaire.joblib")

    print("\nModèle sauvegardé. Le test final n'a pas été évalué.")

    return meilleur_modele, repartition

dossier = Path("resultats_binaire")

sauvegarde = joblib.load(dossier / "modele_binaire.joblib")
modele = sauvegarde["pipeline"]
genes = sauvegarde["genes"]

partitions = pd.read_csv(
    dossier / "partitions.csv",
    dtype={"sample": str, "patient": str},
)

# 2. Retrouver les échantillons de validation uniquement
samples_validation = partitions.loc[
    partitions["partition"] == "validation", "sample"
]

if dtf_modele["sample"].duplicated().any():
    raise ValueError("Des identifiants sample sont dupliqués.")

donnees = dtf_modele.set_index("sample")

samples_absents = set(samples_validation) - set(donnees.index)
if samples_absents:
    raise ValueError(
        f"{len(samples_absents)} échantillons de validation sont absents."
    )

validation = donnees.loc[samples_validation].copy()

# Les gènes sont sélectionnés dans l'ordre attendu par le modèle
X_val = validation.loc[:, genes]

y_val = validation["sample_type"].map({
    "Solid Tissue Normal": 0,
    "Primary Tumor": 1,
})

if y_val.isna().any():
    raise ValueError("Type de prélèvement inattendu dans la validation.")

# 3. Prédire sans réentraîner
predictions = modele.predict(X_val)

position_classe_tumorale = list(modele.classes_).index(1)
scores = modele.predict_proba(X_val)[:, position_classe_tumorale]

# 4. Vérifier que l'on retrouve les résultats précédents
print("\nRÉSULTATS DE VALIDATION")

print(classification_report(
    y_val,
    predictions,
    labels=[0, 1],
    target_names=["Normal", "Tumoral"],
    digits=3,
    zero_division=0,
))

print("Matrice : lignes réelles, colonnes prédites")
print(confusion_matrix(y_val, predictions, labels=[0, 1]))

print(
    "Balanced accuracy :",
    round(balanced_accuracy_score(y_val, predictions), 4),
)
print("ROC-AUC :", round(roc_auc_score(y_val, scores), 4))

# 5. Construire le détail des prédictions
details = validation[
    ["_primary_disease", "sample_type"]
].copy()

details["classe_reelle"] = y_val.to_numpy()
details["classe_predite"] = predictions
details["score_tumoral"] = scores
details["erreur"] = details["classe_reelle"] != details["classe_predite"]

details["type_erreur"] = "Correct"

details.loc[
    (details["classe_reelle"] == 1)
    & (details["classe_predite"] == 0),
    "type_erreur"
] = "Tumeur manquée"

details.loc[
    (details["classe_reelle"] == 0)
    & (details["classe_predite"] == 1),
    "type_erreur"
] = "Fausse alerte"

# 6. Compter les erreurs par maladie ET type de prélèvement
bilan = details.groupby(
    ["_primary_disease", "sample_type"],
    dropna=False,
).agg(
    effectif=("erreur", "size"),
    erreurs=("erreur", "sum"),
)

bilan["erreur_pct"] = 100 * bilan["erreurs"] / bilan["effectif"]

bilan = bilan.sort_values(
    ["erreur_pct", "erreurs"],
    ascending=False,
)

print("\nBILAN PAR MALADIE ET TYPE DE PRÉLÈVEMENT")
print(bilan.round(2).to_string())

# 7. Afficher les échantillons mal classés
erreurs = details.loc[details["erreur"]].copy()

print(f"\nÉCHANTILLONS MAL CLASSÉS : {len(erreurs)}")
print(erreurs.round(4).to_string())

# 8. Sauvegarder les résultats
details.to_csv(dossier / "predictions_validation.csv")
bilan.to_csv(dossier / "bilan_erreurs_validation.csv")
erreurs.to_csv(dossier / "erreurs_validation.csv")

print("\nRésultats sauvegardés dans :", dossier.resolve())