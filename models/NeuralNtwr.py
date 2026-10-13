import torch.nn as nn

class ReseauCancer(nn.Module):
    def __init__(self, nombre_genes, nombre_classes):
        super().__init__()

        self.reseau = nn.Sequential(
            nn.Linear(nombre_genes, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(0.3),

            nn.Linear(128, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Dropout(0.2),

            nn.Linear(64, nombre_classes),
        )

    def forward(self, x):
        return self.reseau(x)

