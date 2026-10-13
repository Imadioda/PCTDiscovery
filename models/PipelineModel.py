from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch

from NeuralNtwr import ReseauCancer


class ChaineCancer:
    def __init__(
        self,
        chemin_binaire=None,
        chemin_cancer=None,
    ):
        # Chemins relatifs au dossier de ce script
        dossier = Path(__file__).resolve().parent

        if chemin_binaire is None:
            chemin_binaire = (
                dossier / "resultats_binaire" / "modele_binaire.joblib"
            )

        if chemin_cancer is None:
            chemin_cancer = (
                dossier / "resultats_cancer" / "modele_cancer.pth"
            )

        # 1. Charger le pipeline binaire
        sauvegarde_binaire = joblib.load(chemin_binaire)

        self.modele_binaire = sauvegarde_binaire["pipeline"]
        self.genes_binaires = sauvegarde_binaire["genes"]
        self.seuil = float(sauvegarde_binaire["seuil"])

        self.colonne_tumorale = list(
            self.modele_binaire.classes_
        ).index(1)

        # 2. Charger le réseau multiclasse
        sauvegarde_cancer = torch.load(
            chemin_cancer,
            map_location="cpu",
            weights_only=True,
        )

        self.genes_cancer = sauvegarde_cancer["genes"]
        self.classes_cancer = np.array(sauvegarde_cancer["classes"])

        self.modele_cancer = ReseauCancer(
            nombre_genes=sauvegarde_cancer["nombre_genes"],
            nombre_classes=sauvegarde_cancer["nombre_classes"],
        )

        self.modele_cancer.load_state_dict(
            sauvegarde_cancer["model_state_dict"]
        )
        self.modele_cancer.eval()

        # Standardisation apprise pendant l'entraînement du réseau
        self.moyennes = (
            sauvegarde_cancer["scaler_mean"].cpu().numpy()
        )
        self.ecarts_types = (
            sauvegarde_cancer["scaler_scale"].cpu().numpy()
        )

    def predire(self, donnees, batch_size=256):
        """
        donnees :
            DataFrame avec une ligne par échantillon
            et une colonne par gène.

            Les colonnes sample, sample_type et _primary_disease
            peuvent être présentes : elles ne servent pas à prédire.

        Retour :
            DataFrame avec prédiction binaire et top-3 des cancers
            uniquement pour les échantillons prédits tumoraux.
        """
        if not isinstance(donnees, pd.DataFrame):
            raise TypeError("donnees doit être un DataFrame pandas.")

        if donnees.empty:
            raise ValueError("Le tableau fourni est vide.")

        if not donnees.columns.is_unique:
            raise ValueError("Les noms des colonnes doivent être uniques.")

        if batch_size < 1:
            raise ValueError("batch_size doit être supérieur ou égal à 1.")

        # 3. Vérifier la présence des gènes attendus
        genes_requis = set(self.genes_binaires) | set(self.genes_cancer)
        genes_absents = genes_requis - set(donnees.columns)

        if genes_absents:
            exemples = sorted(genes_absents)[:10]
            raise ValueError(
                f"{len(genes_absents)} gènes requis sont absents. "
                f"Exemples : {exemples}"
            )

        # Conserver l'identifiant et l'ordre des échantillons
        if "sample" in donnees.columns:
            identifiants = donnees["sample"].to_numpy()
        else:
            identifiants = donnees.index.to_numpy()

        resultats = pd.DataFrame({"sample": identifiants})

        scores_binaires = []

        # 4. Prédiction binaire par lots
        for debut in range(0, len(donnees), batch_size):
            lot = donnees.iloc[debut:debut + batch_size]

            # Sélection ET remise dans l'ordre attendu
            X_binaire = lot.loc[:, self.genes_binaires]

            if not np.isfinite(
                X_binaire.to_numpy(dtype=np.float64)
            ).all():
                raise ValueError(
                    "Les expressions contiennent des NaN ou des infinis."
                )

            # Le pipeline contient déjà son StandardScaler
            scores = self.modele_binaire.predict_proba(X_binaire)[
                :, self.colonne_tumorale
            ]

            scores_binaires.append(scores)

        scores_binaires = np.concatenate(scores_binaires)
        masque_tumoral = scores_binaires >= self.seuil

        resultats["score_tumoral"] = scores_binaires
        resultats["prediction_tissu"] = np.where(
            masque_tumoral, "Tumoral", "Normal"
        )

        # 5. Préparer les colonnes du top-3
        top_k = min(3, len(self.classes_cancer))

        for rang in range(1, top_k + 1):
            resultats[f"cancer_{rang}"] = pd.Series(
                pd.NA,
                index=resultats.index,
                dtype="string",
            )
            resultats[f"score_cancer_{rang}"] = np.nan

        # Positions des seuls échantillons prédits tumoraux
        positions_tumorales = np.flatnonzero(masque_tumoral)

        # 6. Classification du cancer pour ces échantillons
        with torch.no_grad():
            for debut in range(0, len(positions_tumorales), batch_size):
                positions = positions_tumorales[
                    debut:debut + batch_size
                ]

                X_cancer = donnees.iloc[positions].loc[
                    :, self.genes_cancer
                ].to_numpy(dtype=np.float32, copy=True)

                if not np.isfinite(X_cancer).all():
                    raise ValueError(
                        "Expression non finie pour le modèle multiclasse."
                    )

                # Même standardisation que lors de l'entraînement :
                # opérations en place sur une copie float32.
                X_cancer -= self.moyennes
                X_cancer /= self.ecarts_types

                tenseur = torch.from_numpy(X_cancer)

                logits = self.modele_cancer(tenseur)
                probabilites = torch.softmax(logits, dim=1)

                scores_top, indices_top = torch.topk(
                    probabilites,
                    k=top_k,
                    dim=1,
                )

                scores_top = scores_top.numpy()
                indices_top = indices_top.numpy()

                for rang in range(top_k):
                    resultats.loc[
                        positions, f"cancer_{rang + 1}"
                    ] = self.classes_cancer[indices_top[:, rang]]

                    resultats.loc[
                        positions, f"score_cancer_{rang + 1}"
                    ] = scores_top[:, rang]

        return resultats