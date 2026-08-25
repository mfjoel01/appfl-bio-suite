"""Model, loss, and metric for Fed-Heart-Disease. SHIPPED TO WORKERS.

Same import rule as ``dataset.py``: this file's source is read on the driver and executed
on a partner's cluster, so it may import only what a partner has installed -- here, torch
and numpy. No ``appfl_bio_suite`` imports. Enforced by tests/test_shipped_modules.py.

The model is deliberately the FLamby paper's own baseline -- one linear layer -- and not
something better. The scientific claim of this experiment is about the infrastructure:
that a genuine multi-round federated job runs across independently administered clusters
on different continents and reproduces a published number. A 13-parameter logistic
regression on 40 KB of tabular data means every failure is unambiguously an
infrastructure failure and never a compute failure, which is exactly what makes it a
useful vehicle for that claim.
"""

import numpy as np
import torch
import torch.nn as nn
from torch.nn.modules.loss import _Loss


class Baseline(nn.Module):
    """FLamby's Fed-Heart-Disease baseline: logistic regression on 13 features."""

    def __init__(self, input_dim: int = 13, output_dim: int = 1):
        super().__init__()
        self.linear = nn.Linear(input_dim, output_dim)

    def forward(self, x):
        return torch.sigmoid(self.linear(x))


class BaselineLoss(_Loss):
    """Binary cross-entropy.

    Takes probabilities rather than logits, because the model applies the sigmoid
    itself -- so this is BCELoss, not BCEWithLogitsLoss. Swapping in the latter without
    also removing the sigmoid would train against a doubly-squashed output and quietly
    produce a worse model rather than an error.
    """

    def __init__(self):
        super().__init__()
        self.bce = nn.BCELoss()

    def forward(self, input: torch.Tensor, target: torch.Tensor):
        return self.bce(input, target)


def metric(y_true, y_pred):
    """Binary accuracy.

    Accuracy rather than AUC deliberately. The centers are small and highly unbalanced --
    one has a single class present in most batches -- and ``roc_auc_score`` raises on a
    batch with one class, which turns a metric into a crash partway through validation.
    """
    y_true = np.asarray(y_true).astype("uint8")
    try:
        return ((np.asarray(y_pred) > 0.5) == y_true).mean()
    except ValueError:
        return np.nan
