"""Covariate adjustment of federated moments (Frisch-Waugh-Lovell). Coordinator side.

A pooled eQTL analysis adjusts for covariates ``Z`` (sex, genotype PCs, PEER factors) by
regressing them out of **both** the genotypes and the phenotype,

    X_perp = X - Z1 B,   B = (Z1'Z1)^-1 Z1'X
    y_perp = y - Z1 g,   g = (Z1'Z1)^-1 Z1'y

with ``Z1 = [1, Z]``, and fine-mapping ``(X_perp, y_perp)``. Regressing ``Z`` out of ``y``
alone is not the same model: wherever a variant's dosage tracks a covariate (any variant
whose frequency differs across ancestries tracks the genotype PCs), the part of its effect
lying in span(Z) is removed from ``y`` and left in ``X`` as noise. On MAGE that weakened the
published lead variant at ERN1 from p = 9e-10 to p = 1e-3.

Every quantity the fit needs is a sum over samples, so the adjustment federates exactly.
A site that opts in ships, per (locus, ancestry) block, ``Z'X`` (k x M), ``Z'Z`` (k x k) and
``1'Z`` (k), and per instance ``Z'y`` (k), alongside its usual ``G, u, c, q, w, n``. The
coordinator sums the blocks to get the global ``B`` and ``g``, then rewrites each block's
moments as those of its rows of ``(X_perp, y_perp)``:

    G' = G - B'(Z1'X) - (Z1'X)'B + B'(Z1'Z1)B        u' = u - (1'Z1) B
    c' = c - (Z1'X)'g - B'(Z1'y) + B'(Z1'Z1) g        q' = q - (1'Z1) g
    w' = w - 2 g'(Z1'y) + g'(Z1'Z1) g

Because ``B`` and ``g`` are fit on every block the coordinator receives, the result is the
centralized adjustment of exactly that cohort: the full federation reproduces the pooled
FWL fit, and a one-site run reproduces that site adjusting on its own. Downstream code is
unchanged -- it centres and standardizes these moments as before -- except that allele
frequencies (the MAF filter and SuSiEx's ``.frq``) must come from the unadjusted dosages,
so each block keeps its original ``u`` as ``u_raw``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np


@dataclass
class CovariateGeno:
    """Covariate cross-moments of one (site, locus, ancestry) block."""

    zx: np.ndarray  # (k, M)  Z'X
    zz: np.ndarray  # (k, k)  Z'Z
    zu: np.ndarray  # (k,)    1'Z


def _augmented(n: int, block_u: np.ndarray, cov: CovariateGeno):
    """``Z1'Z1`` and ``Z1'X`` for one block, with ``Z1 = [1, Z]``."""
    k = cov.zz.shape[0]
    z1z1 = np.empty((k + 1, k + 1), dtype=np.float64)
    z1z1[0, 0] = float(n)
    z1z1[0, 1:] = cov.zu
    z1z1[1:, 0] = cov.zu
    z1z1[1:, 1:] = cov.zz
    z1x = np.vstack([block_u[None, :], cov.zx])
    return z1z1, z1x


def fwl_adjust_locus(
    geno_blocks: Sequence,
    pheno_by_instance: Mapping[str, Sequence],
    cov_geno: Mapping[tuple[str, str], CovariateGeno],
    cov_pheno: Mapping[tuple[str, str, str], np.ndarray],
) -> None:
    """Rewrite one locus's blocks in place as moments of ``(X_perp, y_perp)``.

    ``geno_blocks`` are the locus's ``GenoAggregate`` blocks, ``pheno_by_instance`` maps each
    of its instances to that instance's ``PhenoAggregate`` blocks. ``cov_geno`` is keyed by
    ``(site, pop)`` and ``cov_pheno`` by ``(site, pop, instance)``.
    """
    if not geno_blocks:
        return
    ids = geno_blocks[0].snp_ids
    for b in geno_blocks:
        if not np.array_equal(b.snp_ids, ids):
            # B is fit per variant over every row; a variant some block lacks has no
            # global coefficient, and dropping it silently would bias the others' rows.
            raise ValueError(
                f"covariate adjustment needs one variant list across blocks; locus "
                f"{b.locus_id} differs at site {b.site}/{b.pop}"
            )
        if (b.site, b.pop) not in cov_geno:
            raise ValueError(f"site {b.site}/{b.pop} sent no covariate moments")

    augmented = {}
    A = None
    XZ = None
    for b in geno_blocks:
        z1z1, z1x = _augmented(b.n, b.u, cov_geno[(b.site, b.pop)])
        augmented[(b.site, b.pop)] = (z1z1, z1x)
        A = z1z1.copy() if A is None else A + z1z1
        XZ = z1x.copy() if XZ is None else XZ + z1x
    # lstsq rather than solve: a site holding one sex, say, makes a covariate column
    # constant there, and the minimum-norm solution is then still the projection.
    B = np.linalg.lstsq(A, XZ, rcond=None)[0]  # (k+1, M)

    gammas = {}
    for inst, blocks in pheno_by_instance.items():
        zy_total = None
        for p in blocks:
            zy1 = np.concatenate([[p.q], cov_pheno[(p.site, p.pop, inst)]])
            zy_total = zy1 if zy_total is None else zy_total + zy1
        gammas[inst] = np.linalg.lstsq(A, zy_total, rcond=None)[0]

    # Phenotype blocks first: their rewrite needs each block's unadjusted Z1'X.
    for inst, blocks in pheno_by_instance.items():
        g = gammas[inst]
        for p in blocks:
            z1z1, z1x = augmented[(p.site, p.pop)]
            if not np.array_equal(p.snp_ids, ids):
                raise ValueError(f"phenotype block {p.site}/{p.pop}/{inst} variant order differs")
            zy1 = np.concatenate([[p.q], cov_pheno[(p.site, p.pop, inst)]])
            z1z1_g = z1z1 @ g
            p.c = p.c - z1x.T @ g - B.T @ zy1 + B.T @ z1z1_g
            p.q = float(p.q - z1z1_g[0])
            p.w = float(p.w - 2.0 * g @ zy1 + g @ z1z1_g)

    for b in geno_blocks:
        z1z1, z1x = augmented[(b.site, b.pop)]
        cross = B.T @ z1x  # (M, M)
        G = b.G - cross - cross.T + B.T @ (z1z1 @ B)
        b.G = 0.5 * (G + G.T)  # exact in arithmetic; symmetrize the rounding
        b.u_raw = b.u
        b.u = b.u - z1z1[0] @ B
