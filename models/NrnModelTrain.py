from pathlib import Path
import random

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import TensorDataset, DataLoader

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
)

from NeuralNtwr import ReseauCancer
from ExpressionDataClean import dtf_modele


def entrainer_cancer(
    dtf_modele,
    fichier_partitions="resultats_binaire/partitions.csv",
    dossier="resultats_cancer",
    epochs=100,
    batch_size=64,
    patience=10,
):
    # ---------------------------------------------------------
    # 1. Paramètres et reproductibilité
    # ---------------------------------------------------------
    graine = 42

    random.seed(graine)
    np.random.seed(graine)
    torch.manual_seed(graine)

    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(graine)

    # CPU pour cette première version, notamment sur ton Mac.
    device = torch.device("cpu")

    dossier = Path(dossier)
    dossier.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------
    # 2. Récupérer les partitions du modèle binaire
    # ---------------------------------------------------------
    partitions = pd.read_csv(
        fichier_partitions,
        dtype={"sample": str, "patient": str},
    )

    if partitions["sample"].duplicated().any():
        raise ValueError("Identifiants dupliqués dans partitions.csv.")

    if partitions.groupby("patient")["partition"].nunique().gt(1).any():
        raise ValueError("Un patient appartient à plusieurs partitions.")

    if dtf_modele["sample"].duplicated().any():
        raise ValueError("Identifiants sample dupliqués dans dtf_modele.")

    # On entraîne sur les vraies tumeurs primaires,
    # indépendamment des prédictions du modèle binaire.
    df = dtf_modele.loc[
        dtf_modele["sample_type"] == "Primary Tumor"
    ].copy()

    genes = [
        colonne for colonne in dtf_modele.columns
        if colonne not in ["sample", "sample_type", "_primary_disease"]
    ]

    if len(genes) != len(set(genes)):
        raise ValueError("Des noms de gènes sont dupliqués.")

    df = df.merge(
        partitions[["sample", "patient", "partition"]],
        on="sample",
        how="left",
        validate="one_to_one",
    )

    if df["partition"].isna().any():
        raise ValueError("Certaines tumeurs n'ont pas de partition.")

    if not df["partition"].isin(["train", "validation", "test"]).all():
        raise ValueError("Nom de partition inattendu.")

    train = df.loc[df["partition"] == "train"].copy()
    val = df.loc[df["partition"] == "validation"].copy()

    if train.empty or val.empty:
        raise ValueError("Partition train ou validation vide.")

    if train["_primary_disease"].isna().any():
        raise ValueError("Diagnostic manquant dans le train.")

    if val["_primary_disease"].isna().any():
        raise ValueError("Diagnostic manquant dans la validation.")

    # ---------------------------------------------------------
    # 3. Encoder les maladies : nom -> entier
    # ---------------------------------------------------------
    encodeur = LabelEncoder()

    y_train = encodeur.fit_transform(train["_primary_disease"])

    classes_inconnues = (
        set(val["_primary_disease"]) - set(encodeur.classes_)
    )

    if classes_inconnues:
        raise ValueError(
            f"Maladies présentes en validation mais absentes du train : "
            f"{classes_inconnues}"
        )

    y_val = encodeur.transform(val["_primary_disease"])

    classes = encodeur.classes_.tolist()
    nombre_classes = len(classes)

    effectifs = pd.DataFrame({
        "train": train["_primary_disease"].value_counts(),
        "validation": val["_primary_disease"].value_counts(),
    }).fillna(0).astype(int).sort_index()

    print("\nEffectifs par maladie :")
    print(effectifs.to_string())

    if (effectifs["validation"] == 0).any():
        print(
            "\nAttention : certaines maladies du train "
            "ne sont pas évaluables sur cette validation."
        )

    effectifs.to_csv(dossier / "effectifs_classes.csv")

    # ---------------------------------------------------------
    # 4. Préparer et standardiser les expressions
    # ---------------------------------------------------------
    X_train = train[genes].to_numpy(dtype=np.float32)
    X_val = val[genes].to_numpy(dtype=np.float32)

    if not np.isfinite(X_train).all():
        raise ValueError("Valeur non finie dans le train.")

    if not np.isfinite(X_val).all():
        raise ValueError("Valeur non finie dans la validation.")

    scaler = StandardScaler()

    # Apprentissage de la standardisation sur le train uniquement.
    X_train = scaler.fit_transform(X_train).astype(np.float32)
    X_val = scaler.transform(X_val).astype(np.float32)

    # ---------------------------------------------------------
    # 5. Conversion en tenseurs et mini-batchs
    # ---------------------------------------------------------
    X_train_tensor = torch.from_numpy(X_train)
    X_val_tensor = torch.from_numpy(X_val)

    y_train_tensor = torch.tensor(y_train, dtype=torch.long)
    y_val_tensor = torch.tensor(y_val, dtype=torch.long)

    train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
    val_dataset = TensorDataset(X_val_tensor, y_val_tensor)

    generateur = torch.Generator().manual_seed(graine)

    # BatchNorm ne peut pas apprendre sur un batch d'un seul élément.
    # On retire le dernier batch uniquement s'il contient 1 échantillon.
    retirer_dernier = len(train_dataset) % batch_size == 1

    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=retirer_dernier,
        num_workers=0,
        generator=generateur,
    )

    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0,
    )

    # ---------------------------------------------------------
    # 6. Modèle, loss et optimiseur
    # ---------------------------------------------------------
    modele = ReseauCancer(
        nombre_genes=len(genes),
        nombre_classes=nombre_classes,
    ).to(device)

    # Première baseline sans pondération des classes.
    # Le macro-F1 et les rappels par maladie permettront
    # d'identifier les éventuelles difficultés sur les classes rares.
    critere = nn.CrossEntropyLoss()

    optimiseur = torch.optim.AdamW(
        modele.parameters(),
        lr=0.001,
        weight_decay=0.0001,
    )

    print(f"\nTrain : {len(train)} tumeurs")
    print(f"Validation : {len(val)} tumeurs")
    print(f"Gènes : {len(genes)}")
    print(f"Classes : {nombre_classes}")
    print(f"Appareil : {device}")

    # ---------------------------------------------------------
    # 7. Entraînement et arrêt anticipé
    # ---------------------------------------------------------
    meilleure_loss = float("inf")
    sans_amelioration = 0
    historique = []
    chemin_modele = dossier / "modele_cancer.pth"

    for epoch in range(1, epochs + 1):
        modele.train()

        somme_loss_train = 0.0
        nombre_train = 0

        for X_batch, y_batch in train_loader:
            X_batch = X_batch.to(device)
            y_batch = y_batch.to(device)

            optimiseur.zero_grad()

            logits = modele(X_batch)
            loss = critere(logits, y_batch)

            loss.backward()
            optimiseur.step()

            somme_loss_train += loss.item() * len(y_batch)
            nombre_train += len(y_batch)

        loss_train = somme_loss_train / nombre_train

        # Validation : aucun apprentissage des poids.
        modele.eval()

        somme_loss_val = 0.0
        predictions_val = []

        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch = X_batch.to(device)
                y_batch = y_batch.to(device)

                logits = modele(X_batch)
                loss = critere(logits, y_batch)

                somme_loss_val += loss.item() * len(y_batch)

                predictions_val.extend(
                    logits.argmax(dim=1).cpu().tolist()
                )

        loss_val = somme_loss_val / len(val_dataset)

        macro_f1 = f1_score(
            y_val,
            predictions_val,
            average="macro",
            labels=np.arange(nombre_classes),
            zero_division=0,
        )

        historique.append({
            "epoch": epoch,
            "loss_train": loss_train,
            "loss_validation": loss_val,
            "macro_f1_validation": macro_f1,
        })

        print(
            f"Epoch {epoch:03d} | "
            f"train loss={loss_train:.4f} | "
            f"val loss={loss_val:.4f} | "
            f"macro-F1={macro_f1:.4f}"
        )

        if loss_val < meilleure_loss:
            meilleure_loss = loss_val
            sans_amelioration = 0

            # Le state_dict contient les poids et les statistiques BatchNorm.
            # Les autres champs permettent de reconstruire le pipeline.
            torch.save({
                "model_state_dict": {
                    nom: valeur.detach().cpu().clone()
                    for nom, valeur in modele.state_dict().items()
                },
                "nombre_genes": len(genes),
                "nombre_classes": nombre_classes,
                "hidden_sizes": [128, 64],
                "dropout": [0.3, 0.2],
                "genes": genes,
                "classes": classes,
                "scaler_mean": torch.tensor(
                    scaler.mean_, dtype=torch.float64
                ),
                "scaler_scale": torch.tensor(
                    scaler.scale_, dtype=torch.float64
                ),
                "epoch": epoch,
                "validation_loss": loss_val,
                "graine": graine,
            }, chemin_modele)

        else:
            sans_amelioration += 1

            if sans_amelioration >= patience:
                print(f"\nArrêt anticipé à l'epoch {epoch}.")
                break

    pd.DataFrame(historique).to_csv(
        dossier / "historique_entrainement.csv",
        index=False,
    )

    # ---------------------------------------------------------
    # 8. Recharger et évaluer le meilleur état sauvegardé
    # ---------------------------------------------------------
    sauvegarde = torch.load(
        chemin_modele,
        map_location=device,
        weights_only=True,
    )

    modele.load_state_dict(sauvegarde["model_state_dict"])
    modele.eval()

    probabilites = []

    with torch.no_grad():
        for X_batch, _ in val_loader:
            logits = modele(X_batch.to(device))

            # Softmax uniquement pour produire les probabilités.
            proba_batch = torch.softmax(logits, dim=1)
            probabilites.append(proba_batch.cpu().numpy())

    probabilites = np.concatenate(probabilites, axis=0)
    predictions = probabilites.argmax(axis=1)

    top_k = min(3, nombre_classes)
    top_indices = np.argsort(-probabilites, axis=1)[:, :top_k]

    score_top_k = (
        top_indices == y_val[:, None]
    ).any(axis=1).mean()

    print("\nMEILLEUR MODÈLE — VALIDATION")
    print("Epoch conservée :", sauvegarde["epoch"])
    print("Accuracy :", round(accuracy_score(y_val, predictions), 4))
    print(
        "Balanced accuracy :",
        round(balanced_accuracy_score(y_val, predictions), 4),
    )
    print(
        "Macro-F1 :",
        round(f1_score(
            y_val,
            predictions,
            average="macro",
            labels=np.arange(nombre_classes),
            zero_division=0,
        ), 4),
    )
    print(f"Top-{top_k} accuracy :", round(float(score_top_k), 4))

    rapport = classification_report(
        y_val,
        predictions,
        labels=np.arange(nombre_classes),
        target_names=classes,
        digits=3,
        zero_division=0,
    )

    print(rapport)

    (dossier / "rapport_validation.txt").write_text(
        rapport,
        encoding="utf-8",
    )

    matrice = confusion_matrix(
        y_val,
        predictions,
        labels=np.arange(nombre_classes),
    )

    pd.DataFrame(
        matrice,
        index=classes,
        columns=classes,
    ).to_csv(dossier / "confusion_validation.csv")

    # Détail des prédictions et des trois premiers choix
    details = val[["sample", "_primary_disease"]].reset_index(drop=True)
    noms_classes = np.array(classes)

    details["prediction"] = noms_classes[predictions]
    details["erreur"] = predictions != y_val

    for rang in range(top_k):
        indices = top_indices[:, rang]

        details[f"choix_{rang + 1}"] = noms_classes[indices]
        details[f"proba_{rang + 1}"] = probabilites[
            np.arange(len(val)), indices
        ]

    details.to_csv(
        dossier / "predictions_validation.csv",
        index=False,
    )

    print("\nModèle sauvegardé :", chemin_modele.resolve())
    print("Le test final n'a pas été utilisé.")

    return modele, sauvegarde


if __name__ == "__main__":
    from ExpressionDataClean import dtf_modele

    modele_cancer, sauvegarde_cancer = entrainer_cancer(dtf_modele)