"""Independent residual methylation MLP; no external model code or weights."""

from torch import nn


class MethylationMLP(nn.Module):
    def __init__(
        self,
        features: int,
        classes: int,
        hidden: int = 256,
        bottleneck: int = 128,
        dropout: float = 0.5,
    ):
        super().__init__()
        self.input = nn.Sequential(
            nn.Linear(features, hidden), nn.LayerNorm(hidden), nn.SiLU(), nn.Dropout(dropout)
        )
        self.residual = nn.Sequential(
            nn.LayerNorm(hidden), nn.Linear(hidden, hidden), nn.SiLU(), nn.Dropout(dropout)
        )
        self.output = nn.Sequential(
            nn.Linear(hidden, bottleneck),
            nn.LayerNorm(bottleneck),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(bottleneck, classes),
        )

    def forward(self, x):
        x = self.input(x)
        return self.output(x + self.residual(x))
