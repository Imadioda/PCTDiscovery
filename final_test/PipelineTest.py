import os
import subprocess
from pathlib import Path
from datetime import datetime
from time import perf_counter

import numpy as np
import pandas as pd

from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    balanced_accuracy_score,
    roc_auc_score,
)


if __name__ == "__main__":
    debut = perf_counter()
    succes = False

    # Tous les chemins partent du dossier contenant ce script.
    dossier = Path(__file__).resolve().parent
    os.chdir(dossier)

    try:
        print("Début :", datetime.now().strftime("%H:%M:%S"), flush=True)
        print("Dossier :", dossier, flush=True)

        # Vérifier les fichiers avant de charger les grosses données.
        fichier_partitions = dossier / "resultats_binaire/partitions.csv"
        fichier_binaire = dossier / "resultats_binaire/modele_binaire.joblib"
        fichier_cancer = dossier / "resultats_cancer/modele_cancer.pth"

        for fichier in [
            fichier_partitions,
            fichier_binaire,
            fichier_cancer,
        ]:
            if not fichier.is_file():
                raise FileNotFoundError(f"Fichier introuvable : {fichier}")

        # Imports après le positionnement dans le bon dossier.
        from PipelineModel import ChaineCancer
        from ExpressionDataClean import dtf_modele

        # 1. Récupérer uniquement le test réservé
        partitions = pd.read_csv(
            fichier_partitions,
            dtype={"sample": str, "patient": str},
        )

        if partitions[["sample", "patient", "partition"]].isna().any().any():
            raise ValueError("Informations manquantes dans les partitions.")

        if partitions["sample"].duplicated().any():
            raise ValueError("Échantillons dupliqués dans les partitions.")

        if partitions.groupby("patient")["partition"].nunique().gt(1).any():
            raise ValueError("Un patient appartient à plusieurs partitions.")

        if dtf_modele["sample"].duplicated().any():
            raise ValueError("Échantillons dupliqués dans les données.")

        ids_test = partitions.loc[
            partitions["partition"] == "test", "sample"
        ]

        if ids_test.empty:
            raise ValueError("La partition test est vide.")

        test = (
            dtf_modele.set_index("sample")
            .loc[ids_test]
            .reset_index()
        )

        y = test["sample_type"].map({
            "Solid Tissue Normal": 0,
            "Primary Tumor": 1,
        })

        if y.isna().any() or y.nunique() != 2:
            raise ValueError("Le test doit contenir des normaux et des tumeurs.")

        y = y.to_numpy(dtype=int)
        tumeurs = y == 1

        if test.loc[tumeurs, "_primary_disease"].isna().any():
            raise ValueError("Diagnostic manquant pour une tumeur.")

        print(
            f"\nTest : {len(test)} échantillons | "
            f"{(~tumeurs).sum()} normaux | {tumeurs.sum()} tumeurs",
            flush=True,
        )

        # 2. Charger et appliquer la chaîne, sans entraînement
        chaine = ChaineCancer(
            chemin_binaire=fichier_binaire,
            chemin_cancer=fichier_cancer,
        )

        inconnues = (
            set(test.loc[tumeurs, "_primary_disease"])
            - set(chaine.classes_cancer)
        )
        if inconnues:
            raise ValueError(f"Cancers non appris présents dans le test : {inconnues}")

        print("Prédictions en cours...", flush=True)
        resultats = chaine.predire(test)

        pred_binaire = (
            resultats["prediction_tissu"] == "Tumoral"
        ).to_numpy(dtype=int)

        # 3. Résultats binaires
        rapport_binaire = classification_report(
            y,
            pred_binaire,
            labels=[0, 1],
            target_names=["Normal", "Tumoral"],
            digits=4,
            zero_division=0,
        )

        matrice_binaire = confusion_matrix(y, pred_binaire, labels=[0, 1])
        balanced_acc = balanced_accuracy_score(y, pred_binaire)
        auc = roc_auc_score(y, resultats["score_tumoral"])

        # 4. Résultats finaux : normal OU type de cancer
        verite = np.full(len(test), "Normal", dtype=object)
        verite[tumeurs] = test.loc[tumeurs, "_primary_disease"].to_numpy()

        prediction = resultats["cancer_1"].fillna("Normal").to_numpy()
        labels = ["Normal"] + chaine.classes_cancer.tolist()

        rapport_final = classification_report(
            verite,
            prediction,
            labels=labels,
            digits=4,
            zero_division=0,
        )

        # Une tumeur manquée au premier étage compte comme une erreur.
        correct_tumeurs = prediction[tumeurs] == verite[tumeurs]

        top_k = min(3, len(chaine.classes_cancer))
        choix = resultats.loc[
            tumeurs,
            [f"cancer_{rang}" for rang in range(1, top_k + 1)],
        ].fillna("").to_numpy()

        correct_top = (choix == verite[tumeurs, None]).any(axis=1)

        rapport = (
            f"Échantillons : {len(test)}\n"
            f"Normaux : {(~tumeurs).sum()} | Tumeurs : {tumeurs.sum()}\n\n"
            f"=== NORMAL / TUMORAL ===\n{rapport_binaire}\n"
            f"Matrice (lignes réelles, colonnes prédites) :\n"
            f"{matrice_binaire}\n"
            f"Balanced accuracy : {balanced_acc:.4f}\n"
            f"ROC-AUC : {auc:.4f}\n\n"
            f"=== CHAÎNE COMPLÈTE ===\n"
            f"Tumeurs manquées au premier étage : "
            f"{((y == 1) & (pred_binaire == 0)).sum()}\n"
            f"Fausses alertes : "
            f"{((y == 0) & (pred_binaire == 1)).sum()}\n"
            f"Tumeurs détectées et correctement identifiées : "
            f"{correct_tumeurs.sum()}/{tumeurs.sum()} "
            f"({100 * correct_tumeurs.mean():.2f} %)\n"
            f"Tumeurs détectées avec vrai cancer dans le top-{top_k} : "
            f"{correct_top.sum()}/{tumeurs.sum()} "
            f"({100 * correct_top.mean():.2f} %)\n\n"
            f"{rapport_final}"
        )

        print("\n" + rapport)

        # 5. Sauvegarder et afficher quelques exemples
        sortie = dossier / "resultats_test_final"
        sortie.mkdir(exist_ok=True)

        resultats["classe_reelle"] = verite
        resultats["prediction_finale"] = prediction
        resultats["correct"] = verite == prediction

        resultats.to_csv(sortie / "predictions_test.csv", index=False)

        (sortie / "rapport_test_final.txt").write_text(
            rapport, encoding="utf-8"
        )

        pd.DataFrame(
            confusion_matrix(verite, prediction, labels=labels),
            index=labels,
            columns=labels,
        ).to_csv(sortie / "confusion_chaine_complete.csv")

        print("\nExemples aléatoires :")
        print(
            resultats.sample(
                n=min(10, len(resultats)), random_state=42
            ).round(4).to_string(index=False)
        )

        print("\nRésultats enregistrés dans :", sortie)
        succes = True

    finally:
        print("\nFin :", datetime.now().strftime("%H:%M:%S"))
        print(f"Durée : {(perf_counter() - debut) / 60:.2f} minutes")

        message = (
            "Le test est terminé."
            if succes
            else "Erreur pendant le test. Consulte la console."
        )

        try:
            subprocess.run(
                [
                    "/usr/bin/osascript",
                    "-e",
                    f'display notification "{message}" with title "Projet Cancer"',
                ],
                check=False,
                timeout=10,
            )
        except (OSError, subprocess.TimeoutExpired):
            print("Notification non disponible.")