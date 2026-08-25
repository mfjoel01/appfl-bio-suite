# The SuSiEx algorithm and its federated formulation

**Companion to** [`design.md`](design.md) (simulator).

**Scope.** Part I reconstructs the original SuSiEx model and its summary-statistic
inference exactly as published (Yuan et al., *Nat. Genet.* 2024) and implemented in the
reference code (`getian107/SuSiEx`, C++ core). Part II derives the federated formulation
we use in FedFM, with the correctness proofs that make federation *exact* rather than an
approximation. Formulas in Part I are cross-checked against the reference source and cite
it by file/line (`model.cpp`, `data.cpp`); nothing here is invented.

> **Status.** Part I describes published, shipping code. Part II is the federation design;
> [`src/fed_fine_mapping.py`](../../../../src/appfl_bio_suite/experiments/fine_mapping/fedfm/fed_fine_mapping.py) implements it against this
> simulator, and [`src/fine_mapping.py`](../../../../src/appfl_bio_suite/experiments/fine_mapping/fedfm/fine_mapping.py) is the centralized
> comparator of §8.1 that Corollary 1 says it must reproduce. The theorems are statements
> about the estimator; they hold for any correct implementation of the SuSiEx recursion,
> including the reference binary run on federated sufficient statistics. §11.1 reports what
> the two paths actually agree to when both are run.

---

## Notation

| symbol | meaning |
|---|---|
| $s = 1,\dots,S$ | populations / ancestries (EUR, AFR, …). In the reference paper one population = one cohort. |
| $k = 1,\dots,K$ | federated **sites** (ANL, Covenant, MBZUAI). Distinct from populations — a site is ancestry-mixed. |
| $M$ | number of variants in the fine-mapping locus |
| $n_s$ | number of individuals in population $s$ |
| $X_s \in \mathbb{R}^{n_s \times M}$ | **standardized** genotype matrix (each column mean-0, variance-1) |
| $y_s \in \mathbb{R}^{n_s}$ | standardized phenotype (mean-0, variance-1) |
| $L$ | number of single effects (upper bound on causal variants), typically 5 or 10 |
| $R_s = \tfrac1{n_s} X_s^\top X_s$ | in-sample LD (correlation) matrix, $\operatorname{diag} = 1$ |
| $\hat\beta_s = \tfrac1{n_s} X_s^\top y_s$ | standardized marginal effect estimates (one GWAS coefficient per SNP) |
| $\tau_s^2$ | prior variance of a single causal effect in population $s$ |
| $\sigma_s^2$ | residual variance in population $s$ |

The reference code stores $\hat\beta_{sj}$ as `beta = (effect/se)/√n = z/√n`
(`data.cpp:311`), $R_s$ as `dat.ld[s]`, $n_s$ as `ngwas[s]`, and
$\tau_s^2 = \max_j \hat\beta_{sj}^2$ as `tau_sq[s]` (`data.cpp:482`).

---

# Part I — The original SuSiEx

## 1. Single-population foundation: the Sum of Single Effects (SuSiE)

SuSiEx is a cross-population generalization of SuSiE (Wang et al., *JRSS-B* 2020). SuSiE
writes the regression coefficient vector of a locus as a **sum of $L$ single-effect
vectors**, each nonzero in exactly one coordinate:

$$
y = X\beta + \epsilon, \qquad
\beta = \sum_{l=1}^{L} b_l, \qquad
b_l = \gamma_l\, b_l, \quad
\gamma_l \sim \mathrm{Mult}(1, \pi), \quad
b_l \sim \mathcal N(0, \tau^2).
$$

Here $\gamma_l \in \{0,1\}^M$ is a one-hot indicator ("which SNP does effect $l$ act
through") and $\pi = (\pi_1,\dots,\pi_M)$ is the prior probability that each SNP is the
causal one for that effect (uniform $\pi_j = 1/M$ by default). Because each $b_l$ has a
*single* nonzero entry, the posterior on $\gamma_l$ is a categorical distribution
$\alpha_l = (\alpha_{l1},\dots,\alpha_{lM})$ — the per-effect **posterior inclusion
probabilities (PIPs)** — and a level-$\rho$ **credible set** is the smallest set of SNPs
whose $\alpha_l$ mass reaches $\rho$.

The model is fit by **Iterative Bayesian Stepwise Selection (IBSS)**, a coordinate-ascent
variational inference that cycles over the $L$ effects. For effect $l$ it (i) forms the
residual by subtracting the fitted contribution of all *other* effects, (ii) solves a
**Single Effect Regression (SER)** — a Bayesian one-variable regression over all $M$ SNPs
— and (iii) adds the new effect back. IBSS maximizes the Evidence Lower BOund (ELBO) and
converges in a handful of passes.

## 2. The cross-population extension (SuSiEx model)

SuSiEx couples $S$ population-specific SuSiE models by **sharing the causal indicators
$\gamma_l$ across populations while letting the effect magnitudes differ** (paper, Methods,
"The Cross-population Sum of Single Effects model"):

$$
y_s = X_s\beta_s + \epsilon_s, \qquad \epsilon_s \sim \mathcal N(0,\sigma_s^2 I),
\quad s = 1,\dots,S,
$$

$$
\beta_s = \sum_{l=1}^{L} b_{sl}, \qquad
b_{sl} = \gamma_l\, b_{sl}, \qquad
\gamma_l \sim \mathrm{Mult}(1,\pi), \qquad
b_{sl} \sim \mathcal N(0, \tau_{sl}^2).
$$

The single, load-bearing modelling assumption:

> $\gamma_l$ **does not depend on $s$** (the causal SNP is *shared* across populations),
> but $b_{sl}$ **does depend on $s$** (the effect size, including a null effect, may vary
> by ancestry).

This is why SuSiEx reports **one PIP per variant** across all populations (a single
credible set, not population-specific credible sets), yet still estimates a
**population-specific causal probability** per credible set. It also permits a variant to
be *absent* in a population (e.g. monomorphic / below MAF); that population simply does not
contribute to the variant's evidence.

## 3. Summary-statistic inference (the actual recursion)

The decisive practical property — the one that makes federation almost free — is stated in
the paper:

> *"Both the IBSS algorithm and the ELBO can be computed using GWAS summary statistics, and
> thus the SuSiEx model can be fitted without access to any individual-level data."*

Concretely, the fit consumes only $\{\hat\beta_s, R_s, n_s\}_{s=1}^S$ (plus derived
$\tau_s^2$). Below is the exact recursion as implemented in `model.cpp`.

### 3.1 Residualization (IBSS outer loop)

Maintain the fitted marginal contribution of all effects,
$\widehat{m}_{sj} = \big(R_s \textstyle\sum_{l} b_{sl}\big)_j$ (code: `bhat`). For the
current effect $l$, remove its own contribution to form the **residualized marginal
association**

$$
r_{sj}^{(l)} \;=\; \hat\beta_{sj} \;-\; \Big(R_s \sum_{l'\neq l} b_{sl'}\Big)_j
\qquad\text{(code: \texttt{beta\_ll = beta - bhat}, model.cpp:152–153).}
$$

This is the summary-statistic analogue of "subtract the other effects' fitted values" — the
LD matrix $R_s$ plays the role of $\tfrac1{n_s}X_s^\top X_s$.

### 3.2 Single Effect Regression (per population, per SNP)

For each population $s$ and SNP $j$, with sampling variance $v_{sj}^2 = \sigma_s^2/n_s$
(set to $10^6$, i.e. "no information", when the SNP is absent in $s$ — code `ind[i]`):

$$
v_{sj}^2 = \frac{\sigma_s^2}{n_s}, \qquad
z_{sj} = \frac{r_{sj}^{(l)}}{v_{sj}},
$$

$$
\boxed{\;\lambda_{sj} \;=\; \underbrace{\tfrac12 \log\!\frac{v_{sj}^2}{v_{sj}^2 + \tau_s^2}
\;+\; \tfrac12\, z_{sj}^2\,\frac{\tau_s^2}{v_{sj}^2 + \tau_s^2}}_{\text{Wakefield log approximate Bayes factor}}\;}
\qquad\text{(model.cpp:264)}
$$

The single-effect **posterior** for SNP $j$ in population $s$ (Gaussian conjugate update):

$$
\phi_{sj}^2 = \Big(\tfrac1{v_{sj}^2} + \tfrac1{\tau_s^2}\Big)^{-1}, \qquad
\mu_{sj} = \frac{\phi_{sj}^2}{v_{sj}^2}\, r_{sj}^{(l)}, \qquad
\mathbb E[b_{sj}^2] = \phi_{sj}^2 + \mu_{sj}^2
\qquad\text{(model.cpp:267–269).}
$$

### 3.3 Cross-population coupling — the SuSiEx step

This is where the $S$ populations are *linked*. Because the populations are conditionally
independent given the shared causal SNP, the **combined evidence** that SNP $j$ carries
effect $l$ is the **product of per-population Bayes factors**, i.e. the **sum of log-ABFs**:

$$
\boxed{\;\Lambda_j \;=\; \sum_{s=1}^{S} \lambda_{sj}\;}
\qquad\text{(code: \texttt{w[j] += curLogBF}, model.cpp:272).}
$$

The posterior inclusion probability for effect $l$ is then a prior-weighted softmax over the
combined evidence:

$$
\alpha_{lj} \;=\; \frac{\pi_j\, e^{\Lambda_j}}{\sum_{k=1}^{M} \pi_k\, e^{\Lambda_k}}
\qquad\text{(model.cpp:288–297, computed in a max-stabilized form).}
$$

The single-effect posterior mean and second moment that feed the next IBSS pass are the
PIP-weighted per-population quantities:

$$
b_{slj} = \alpha_{lj}\,\mu_{sj}, \qquad
\mathbb E[b_{slj}^2] = \alpha_{lj}\,\mathbb E[b_{sj}^2]
\qquad\text{(model.cpp:305–307).}
$$

Note the split of roles: **$\alpha_{lj}$ is shared across populations** (one PIP — the
combined $\Lambda_j$), while **$\mu_{sj}$ is population-specific**. This is the algorithmic
realization of "shared $\gamma_l$, population-specific $b_{sl}$".

### 3.4 Variance updates and the ELBO

The residual variances are re-estimated each pass from the expected residual sum of squares
(ERSS), expressible entirely in summary statistics:

$$
\mathrm{ERSS}_s = n_s\Big(1 - 2\, {\textstyle\sum_l} b_{sl}^\top \hat\beta_s
+ {\textstyle\sum_{l,l'}} b_{sl}^\top R_s b_{sl'} - {\textstyle\sum_l} b_{sl}^\top R_s b_{sl}
+ {\textstyle\sum_l}\mathbb E[b_{sl}^2]\Big), \qquad
\sigma_s^2 = \frac{\mathrm{ERSS}_s}{n_s}
$$

(code `erss`, model.cpp:205,212). Convergence is monitored on the ELBO,

$$
\mathrm{ELBO} = \sum_{l=1}^{L}\big(\mathcal L_l - \tilde{\mathcal L}_l\big)
+ \sum_{s=1}^{S}\Big(-\tfrac12 n_s\log(2\pi\sigma_s^2) - \tfrac{\mathrm{ERSS}_s}{2\sigma_s^2}\Big),
$$

with the per-effect log-likelihood term $\mathcal L_l = \big(\max_s\lambda\big) + \log\sum_j w_j$
and its posterior counterpart $\tilde{\mathcal L}_l$ (model.cpp:214–239, 312–320); IBSS stops
when $|\Delta\mathrm{ELBO}| < $ `tol`.

### 3.5 Prior variance $\tau_s^2$

Rather than an EM update, the reference sets the per-population prior variance to the largest
squared standardized marginal effect in the locus, $\tau_s^2 = \max_j \hat\beta_{sj}^2$
(`data.cpp:482`) — a data-driven, conservative "the strongest signal bounds the effect
scale" choice. (This detail matters for federation: see §9.)

## 4. Post-processing: credible sets, purity, causal probability

Given a converged fit (`model.cpp: cal_pip`):

- **Credible sets.** For each effect $l$, sort SNPs by $\alpha_{lj}$ descending and admit
  until cumulative mass exceeds the level $\rho$ (`par.level`, e.g. 0.95 or 0.99).
- **Purity filter.** A credible set's **purity** is the *minimum absolute pairwise LD* among
  its SNPs, taken as the **worst (max) over populations**:
  $\text{purity} = \max_s \min_{j,j' \in \mathrm{CS}} |R_{s,jj'}|$. Sets below `min_purity`
  (default 0.5) — diffuse, inferentially useless effects that appear when $L$ exceeds the
  true number of signals — are discarded, as are sets with no genome-wide-significant SNP
  (`minP > pth`).
- **Overall PIP** of a SNP aggregates across the retained effects,
  $\mathrm{PIP}_j = 1 - \prod_l (1 - \alpha_{lj})$ (model.cpp:25).
- **Population-specific causal probability.** For each retained credible set, SuSiEx
  estimates, per population, the probability that the set has a nonzero effect in that
  population (thresholded at 0.8 for a binary "causal in population $s$" call). A small value
  is a power statement, not evidence of no effect.
- **Choosing $L$.** A heuristic multi-step schedule: start $L=5$; if it fails to converge,
  decrement until it does; if it converges *and* returns 5 sets (suggesting more exist), jump
  to $L=10$ and repeat.

---

# Part II — The federated formulation (FedFM)

## 5. Setting and threat model

In FedFM the individuals live at $K=3$ institutions (ANL, Covenant, MBZUAI; see
[`design.md` §3](design.md)). Data-residency rules forbid pooling individual-level genotypes
or phenotypes. Each site may compute statistics on its own individuals and send **only**
those statistics to a coordinator, which runs the SuSiEx recursion and returns credible sets.
No raw genotype row ever leaves a site.

The pivotal observation is that **SuSiEx is already a summary-statistic method** (§3): its
entire recursion touches the data only through $\{\hat\beta_s, R_s, n_s, \tau_s^2\}$. So
federation is not a new inference algorithm — it is a statement about *which quantities are
computed where*, plus a proof that the coordinator's fit on aggregated statistics equals the
fit on pooled data. We make that precise.

## 6. Lemma (sufficiency): the fit factors through per-population summary statistics

**Lemma 1.** Every quantity produced by a SuSiEx iteration —
$r_{sj}^{(l)}, v_{sj}^2, z_{sj}, \lambda_{sj}, \Lambda_j, \phi_{sj}^2, \mu_{sj}, \alpha_{lj},
b_{slj}, \sigma_s^2$, and the ELBO — is a deterministic function of

$$
T_s \;=\; \big(R_s,\ \hat\beta_s,\ n_s,\ \tau_s^2\big),\quad s=1,\dots,S,
$$

together with the current parameter state $\{b_{sl},\sigma_s^2\}$. The individual-level
$(X_s, y_s)$ never appear otherwise.

*Proof.* Structural induction over the recursion of §3. The only data arrays read in
`iter()`/`ser()` are `dat.ld` $= R_s$, `beta` $= \hat\beta_s$, `ngwas` $= n_s$, and `tau_sq`
$= \tau_s^2$; the ERSS (§3.4) expands to $n_s$, $\hat\beta_s$, and $R_s$ contracted with the
current $b_{sl}$. No step forms $X_s$ or $y_s$. $\square$

Lemma 1 is the whole game: to federate, it suffices to produce each $T_s$ from
site-resident data and ship $T_s$ to the coordinator.

## 7. Theorem (exact horizontal federation): summary statistics are additive over sites

Populations must be assembled from individuals scattered across sites. Write the individuals
of population $s$ as a disjoint union over sites, $\mathcal D_s = \bigsqcup_{k} \mathcal
D_{s,k}$, with $X_{s,k}, y_{s,k}$ the (raw-dosage) genotype/phenotype blocks held by site
$k$. Define the **raw second-moment aggregates**

$$
G_{s,k} = X_{s,k}^\top X_{s,k}\ (M\times M), \quad
c_{s,k} = X_{s,k}^\top y_{s,k}\ (M), \quad
u_{s,k} = \mathbf 1^\top X_{s,k}\ (M),
$$
$$
n_{s,k} = |\mathcal D_{s,k}|, \quad
q_{s,k} = \mathbf 1^\top y_{s,k}, \quad
w_{s,k} = y_{s,k}^\top y_{s,k}.
$$

**Theorem 1 (exactness).** The pooled sufficient statistics are the **sums** of the
per-site aggregates,

$$
X_s^\top X_s = \sum_k G_{s,k},\quad
X_s^\top y_s = \sum_k c_{s,k},\quad
\mathbf 1^\top X_s = \sum_k u_{s,k},\quad
n_s = \sum_k n_{s,k},\ \dots
$$

and $T_s = (R_s, \hat\beta_s, n_s, \tau_s^2)$ is an exact algebraic function of these sums
(after column standardization; §9). Consequently, **under the harmonization conditions
below**, the coordinator's SuSiEx fit on $\{T_s\}$ is **identical** (up to floating-point
summation order) to the fit obtained if all individual-level data had been pooled centrally.

**Conditions.** The additivity of the sums is unconditional linear algebra; equality
*against a centralized comparator* additionally requires that the pooled matrix it is
compared to be well-defined identically at every site:
(i) a common ordered, allele-harmonized variant list — a flipped reference allele negates a
standardized column and corrupts $X_s^\top y_s$ and the off-diagonals of $X_s^\top X_s$;
(ii) a common complete-case cohort per ancestry, *or* identical missing-data handling —
sporadic per-SNP missingness makes a single $n_p$ insufficient unless columns are imputed to
a consistent (global) value, as the FedFM simulator does by mean imputation
([`design.md` §4.2](design.md));
(iii) identical phenotype transform, dosage coding, QC filtering, and SuSiEx initialization to
the comparator;
(iv) if covariates are adjusted, a globally-equivalent adjustment — site-specific
residualization against differing covariate designs does not reproduce pooled regression;
(v) nondegeneracy — every retained variant has positive pooled variance, with monomorphic
variants excluded or carried under the reference "absent / no-information" convention (§3.2).
These are data-harmonization requirements, not extra assumptions on the estimator.

**Conditions in practice.** Two of them are *not* free on the FedFM package, and both were
found by running the federated path against the comparator rather than by reading the data:

* **(i) is violated at rest.** The per-site filesets were cut with `plink --keep --make-bed`
  and *without* `--keep-allele-order` ([`design.md` §3](design.md)), so PLINK 1.9 assigned A1
  by each site's *own* minor allele. Against the HAPNEST coding, 9,398 / 31,031 / 9,777 of
  533,532 chr1 variants are reversed at ANL / Covenant / MBZUAI — about 6% of variants differ
  in orientation between ANL and Covenant. The centralized path never sees this, because it
  cuts its window straight from the HAPNEST fileset; the federated path reads the site
  filesets, where the flipped site contributes $2-x$ where the others contribute $x$. This is
  trap 2 of §A.10, and it sums in silently. `fed_fine_mapping.reference_variants` publishes a
  canonical variant coding, each site recodes its own dosages to it before forming any
  moment, and `pool_geno` re-checks rather than trusting the sites. The upstream fix is to
  re-cut the filesets with `--keep-allele-order`.
* **(ii) has a residual cost.** The comparator's missing-data handling is *per-variant* for
  the marginals (`plink2 --glm` regresses on that variant's complete cases) and *pairwise*
  for LD (`plink --r` correlates on pairwise-complete observations). Neither is expressible
  in a single $(G_{p,k}, u_{p,k}, n_{p,k})$ triple — reproducing them exactly would take four
  $M\times M$ matrices per block instead of one. HAPNEST's missing-call rate is ~2×10⁻⁷ and
  concentrated in a handful of variants, so the implementation instead drops any window
  variant a block cannot observe completely, keeping every retained variant's moments exact
  over the block's *full* cohort, and reports the count. At L0000 that dropped 6 / 1 / 3 / 1
  window variants from the EUR / AMR / CSA / MID columns, of which 2 / 1 / 1 / 1 would
  otherwise have survived the MAF filter. Where the count is nonzero the two paths fine-map
  slightly different variant sets and Corollary 1's equality is no longer exact — the
  alternative, dropping the ~0.15% of *individuals* who carry a missing call anywhere in the
  window, perturbs *every* variant's $\hat\beta_p$ instead and is strictly worse.

*Proof.* A matrix product over a set of rows decomposes over any partition of those rows:
$X_s^\top X_s = \sum_i x_i x_i^\top = \sum_k \sum_{i\in\mathcal D_{s,k}} x_i x_i^\top = \sum_k
G_{s,k}$, and identically for $X_s^\top y_s$, $\mathbf 1^\top X_s$, $\mathbf 1^\top y_s$,
$y_s^\top y_s$, $n_s$. The standardized $R_s$ and $\hat\beta_s$ are the usual functions of
these central moments (§9, eqns 9.1–9.2). Compose with Lemma 1: the fit is a function of
$T_s$, and each $T_s$ is a function of $\sum_k(\cdot)$. Hence the federated fit equals the
pooled fit. $\square$

**Remark (one round, then silence).** The site $\to$ coordinator transfer happens **once**
(the aggregates are static). The IBSS iterations of §3 then run entirely on the coordinator
with **zero** further site communication. This is unlike federated gradient/EM schemes, which
require a synchronization round per iteration. SuSiEx federates in a *single* round.

## 8. Populations = ancestries, pooled across sites

Theorem 1 lets us assemble a population from many sites. But *what is a population* here?
FedFM's scientific invariant is that **effect sizes are indexed by superpopulation, not by
site** ([`design.md` §5.2](design.md)): an AFR individual at ANL and an AFR individual at
Covenant share the same causal effect $\beta$ by construction. So a population is an
**ancestry**, and its cohort is that ancestry's individuals *wherever they live* — the AFR
population is all 60,000 AFR individuals across the three sites, not three AFR sub-cohorts.

Take $S$ = number of superpopulations. Each site $k$ partitions *its own* individuals by
superpopulation and computes, **locally**, the per-(site, ancestry) aggregates
$G_{p,k}, c_{p,k}, n_{p,k}, \dots$ of §7 for each ancestry $p$ it holds. The coordinator sums
over sites,

$$
G_p = \sum_{k} G_{p,k}, \qquad c_p = \sum_k c_{p,k}, \qquad n_p = \sum_k n_{p,k},
$$

forms $(R_p, \hat\beta_p, n_p)$, and runs SuSiEx with **populations = ancestries**.

**Corollary 1.** The federated fit above exactly reconstructs the SuSiEx fit that would be
obtained by pooling *all* individuals of ancestry $p$ across every site into a single
per-ancestry GWAS — precisely the estimand matched to the simulator's per-superpopulation
effect model.

*Proof.* Apply Theorem 1 with the row partition of ancestry-$p$ individuals given by their
site of residence, $\mathcal D_p = \bigsqcup_k \mathcal D_{p,k}$. $\square$

This is exactly federatable, privacy-compatible (only per-ancestry aggregates leave a site),
and estimand-correct. It is also the formulation against which the ground-truth causal sets
(`causal_manifest.tsv`) are the right yardstick, since those are defined per superpopulation.

The alternative of taking **populations = sites** ($S = K = 3$, one GWAS per ancestry-mixed
site cohort) is trivially federated — it *is* the reference API with cohort = population — but
each $\hat\beta_k$ is then a *mixture* marginal and each $R_k$ *mixed-ancestry* LD, so the
single-effect model of §2 is misspecified against the per-ancestry ground truth and $\tau_k^2$
conflates ancestries. It is valid but discards the structure the simulator encodes, so FedFM
does not use it.

### 8.1 The centralized baseline (what is built today)

Corollary 1 has a useful consequence: because the federated fit *equals* the pooled
per-ancestry fit, the pooled fit can be computed directly — with no site boundaries at all —
and serves as the **centralized reference point** the federated implementation must reproduce.

That is what `src/fine_mapping.py` currently does. Rather than have sites emit aggregates and
sum them (§7, §9), it pools each ancestry's individuals across sites outright and runs one
GWAS + one in-sample LD panel per ancestry over the full pooled cohort:

| column | pooled $n_p$ | contributing sites |
|---|---:|---|
| EUR | 36,500 | ANL 30,000 + Covenant 1,500 + MBZUAI 5,000 |
| AFR | 60,000 | ANL 7,500 + Covenant 47,500 + MBZUAI 5,000 |
| AMR | 7,500 | ANL |
| EAS | 2,500 | ANL |
| CSA | 18,500 | ANL 2,500 + Covenant 1,000 + MBZUAI 15,000 |
| MID | 25,000 | MBZUAI |

$S = 6$ columns covering all 150,000 enrolled individuals exactly once. `min_gwas_n` drops any
ancestry whose *pooled* cohort is too small to inform a GWAS (at 1,000 all six are kept).

Because it is computed centrally, this run is the **upper bound** on what the federated
implementation can achieve: every ancestry gets its maximum attainable sample size, and each
LD matrix is the true in-sample LD of the full pooled cohort. By Corollary 1 an exact
federated implementation should *match* it, not merely approach it; any shortfall measures
implementation loss (or violated harmonization conditions from Theorem 1), which is precisely
what makes it the right baseline to report against.

Two honest caveats on "best possible". First, pooling only helps ancestries that actually span
sites — AMR (7,500) and EAS (2,500) live solely at ANL, so they gain nothing and remain the
weak columns. Second, the ceiling is a ceiling *for this estimand*: it assumes the shared
causal indicator $\gamma_l$ of §2, which the `ancestry_divergent_causal` instances
([`design.md` §5.4](design.md)) deliberately violate, so those instances are not guaranteed to
be best-case under any column construction.

## 9. Standardization reconciliation (the one real subtlety)

SuSiEx assumes **column-standardized** $X_s$ (variance 1) and standardized $y_s$. If each
site standardizes locally — using its *own* per-column mean $\mu_{jk}$ and SD $\sigma_{jk}$ —
then the locally standardized blocks are **not** on a common scale, and their Gram matrices
do **not** sum to the globally standardized Gram. Naively summing locally standardized
$R_{p,k}$ is therefore *biased*.

The fix is to ship **raw (unstandardized) second-moment aggregates** and standardize once, at
the coordinator, against the *pooled* moments. With pooled per-column count $n_p$, sum
$U_{pj} = \sum_k (u_{p,k})_j$, sum-of-squares $Q_{pj} = \sum_k (G_{p,k})_{jj}$, and pooled
cross-moments, the exact global quantities are

$$
\bar x_{pj} = \frac{U_{pj}}{n_p}, \qquad
\widehat{\operatorname{sd}}(x_{pj}) = \sqrt{\tfrac{Q_{pj}}{n_p} - \bar x_{pj}^2},
\tag{9.1}
$$

$$
R_{p,jj'} = \frac{(G_p)_{jj'} - n_p\,\bar x_{pj}\,\bar x_{pj'}}
{n_p\,\widehat{\operatorname{sd}}(x_{pj})\,\widehat{\operatorname{sd}}(x_{pj'})}, \qquad
\hat\beta_{pj} = \frac{(c_p)_j - n_p\,\bar x_{pj}\,\bar y_p}
{n_p\,\widehat{\operatorname{sd}}(x_{pj})\,\widehat{\operatorname{sd}}(y_p)}.
\tag{9.2}
$$

Every term on the right is a **sum of per-site aggregates** (Theorem 1), so (9.2) is exact
and requires no raw rows. This is the recommended primitive: sites emit
$(G_{p,k}, c_{p,k}, u_{p,k}, n_{p,k}, q_{p,k}, w_{p,k})$; the coordinator applies (9.1)–(9.2)
then §3.

Two practical notes:
1. **Prior variance.** $\tau_p^2 = \max_j \hat\beta_{pj}^2$ (§3.5) is computed at the
   coordinator *after* pooling, so it too is exact — sites need not agree on it in advance.
2. **Shared-reference shortcut (with a caveat).** If instead each site standardizes against a
   *common* external reference (same $\mu_j, \sigma_j$ everywhere), the standardized **Gram**
   matrices are additive, $Z^\top Z = \sum_k Z_k^\top Z_k$ — but a **correlation / LD** matrix
   $R_k = Z_k^\top Z_k / n_k$ is **not**: the coordinator must recombine it by the
   $n_k$-weighted average
   $R = \tfrac1n \sum_k n_k R_k$, $n = \sum_k n_k$, *not* a plain sum. So a site may ship
   either $Z_k^\top Z_k$ (a Gram — sum, then divide by $n$) or $R_k$ (an LD matrix — average
   with $n_k$ weights), but "ship an LD matrix and sum" is wrong. Even done correctly, the
   pooled diagonal is only approximately unit (the external $\sigma_j \neq$ pooled SD), so this
   route only *approximately* matches in-sample standardization; the FedFM simulator computes
   per-site in-sample LD with the population-variance convention
   ([`design.md` §4.2](design.md)), so the exact route (9.2) is preferred.

### 9.1 Feeding (9.1)–(9.2) through the reference binary

Eqs (9.1)–(9.2) are the quantities the *model* consumes. The reference implementation does
not accept them directly: it reads a GWAS summary-statistics file and reconstructs
$\hat\beta_{pj} = (\text{BETA}/\text{SE})/\sqrt{n_p}$ (`data.cpp:311`), and it reads an LD
panel as three files on disk rather than as a matrix. A coordinator holding only moments has
to produce both, and neither is a re-derivation of the estimator — only of its file
interface.

**Summary statistics.** The simple regression of $y_p$ on one dosage column with an
intercept has a closed form in exactly the moments the sites shipped. With
$S_{xx} = (G_p)_{jj} - u_{pj}^2/n_p$, $S_{xy} = (c_p)_j - u_{pj} q_p / n_p$ and
$S_{yy} = w_p - q_p^2/n_p$,

$$
\text{BETA}_j = \frac{S_{xy}}{S_{xx}},\qquad
\text{SE}_j = \sqrt{\frac{S_{yy} - \text{BETA}_j S_{xy}}{(n_p-2)\,S_{xx}}},\qquad
t_j = \frac{\text{BETA}_j}{\text{SE}_j}\ \text{on } n_p-2 \text{ df},
\tag{9.3}
$$

which is what `plink2 --glm` computes, to the last digit it prints. Note that (9.2) and
(9.3) are *not* the same number: $\hat\beta_{pj}$ of (9.2) is the Pearson correlation
$S_{xy}/\sqrt{S_{xx}S_{yy}}$, while the binary forms $t_j/\sqrt{n_p}$; they differ by the
usual $t$-vs-correlation factor of order $1/n_p$. That gap is the reference implementation's
own approximation and is present identically on both paths, so it does not affect the
comparison — but a coordinator that wrote (9.2) into the BETA column instead of (9.3) would
*not* be running the same estimator as the comparator.

**LD panel.** SuSiEx probes for `{prefix}.ld.bin`, `{prefix}_ref.bim` and
`{prefix}_frq.frq`; finding all three it takes its precomputed path and never touches a
genotype panel (`main.cpp:204–222`). The binary is $M\times M$ float32, row-major, and its
size is checked against the `_ref.bim` line count (`data.cpp:522`); the `.frq` must repeat
each bim line's alleles in order and carries the A1 frequency, which is $u_{pj}/2n_p$
(`data.cpp:163–183`). So the coordinator writes $R_p$ from (9.2) and the frequencies from
(9.1) directly, and the fit it obtains is the fit the centralized run obtains from
`plink --r square bin4` on the pooled genotypes.

Both are implemented in `fed_fine_mapping.pooled_sumstats` and
`fed_fine_mapping.write_ld_panel`. §11.1 reports how closely they reproduce the comparator's
own files.

## 10. Communication and privacy

Per fit, per ancestry block a site holds, the site sends an $M\times M$ Gram (or LD) matrix
plus $O(M)$ vectors:

$$
\text{uplink per (site, ancestry)} \;=\; O(M^2)\ \text{floats}, \qquad
\text{rounds} = 1 .
$$

The $O(M^2)$ LD matrix is the dominant term — but it is exactly what *any* summary-statistic
fine-mapper already transmits; federation adds no new class of disclosure. What crosses the
boundary is aggregate LD and marginal associations, i.e. standard GWAS-sharing artifacts.
These support the usual (bounded) inference risks — membership / LD-reconstruction attacks —
but **not** recovery of individual genotype rows, which are never formed off-site. The IBSS
iterations disclose nothing further, since they are coordinator-local (§7 remark).

If even aggregate LD sharing is disallowed, one can additionally mask $G_{p,k}$ with
calibrated noise (a differential-privacy Gaussian mechanism on the Gram) at the cost of
converting Theorem 1's *exact* equality into a controlled approximation — a natural extension,
deliberately out of scope here.

## 11. Mapping to the pipeline

The federation primitives line up with existing FedFM stages:

| SuSiEx / federation need | FedFM artifact |
|---|---|
| per-site, per-ancestry genotype blocks $X_{p,k}$ | `data/processed/{site}/{site}_chr1.{bed,bim,fam}` + per-individual superpop in `{site}_manifest.tsv` ([`design.md` §3](design.md)) |
| pooled per-ancestry LD $R_p$ | built at fine-mapping time — the locus window restricted to the ancestry's pooled individuals by `plink --keep`, then `plink --r square` (divisor-$n$ convention, [`design.md` §4.2](design.md)) |
| locus variant lists $M$ | in-window variants taken by bp range from `selected_loci.tsv` (`data/loci/tag_snps/` is unused) |
| per-(site, locus, arch, rep) phenotypes $y$ | `data/ground_truth/phenotypes/{site}/*.pheno`, concatenated across sites per instance |
| per-individual ancestry (to form the pooled columns) | `data/processed/{site}/{site}_manifest.tsv` |
| ground-truth causal sets to score PIPs against | `data/ground_truth/causal_manifest.tsv` |
| canonical variant coding (Thm 1 cond. (i)) | `data/raw/hapnest/chr1.bim`, read once by `fed_fine_mapping.reference_variants` |
| pooled per-ancestry LD $R_p$ **without genotypes** | written directly by `fed_fine_mapping.write_ld_panel` (§9.1) |

Two modules, one estimator.

`src/fine_mapping.py` implements the centralized baseline of §8.1 for each (locus,
architecture, replicate): it pools each ancestry across sites, runs one GWAS + one in-sample
LD panel per ancestry, and runs SuSiEx across the six ancestry columns. It computes the pooled
fit *directly* rather than by the aggregate-and-sum route of §7/§9 — legitimate because
Corollary 1 says the two coincide. The returned credible sets are scored against the
ground-truth causal SNPs to produce calibration curves and credible-set-size distributions as
a function of $(h^2, r_g,\ \text{ancestry composition},\ \text{LD-divergence stratum})$.

`src/fed_fine_mapping.py` implements §7 and §9. **Site stage:** each site reads only its own
`{site}_chr1.{bed,bim,fam}`, manifest and `.pheno` files, recodes to the canonical allele
order, and emits the six raw aggregates of §7 per (site, ancestry) — twelve blocks in all
(ANL holds five ancestries, Covenant three, MBZUAI four). Because only $c_{p,k}, q_{p,k},
w_{p,k}$ depend on the phenotype, the $O(M^2)$ half goes once per locus and the $O(M)$ half
once per instance; that is a transmission economy, not a change to what §7 says crosses the
boundary. **Coordinator stage:** sum over sites (§7), then standardize *once* against the
pooled moments (§9) — including the MAF filter, applied to pooled $u_{pj}/2n_p$ — then write
the summary statistics and the LD panel per §9.1 and invoke the same SuSiEx binary with the
same arguments. Result rows carry the schema `fm_results.tsv` carries, so the two tables are
directly comparable. Entry point `scripts/run_fed_fine_mapping.py`; outputs under
`reports/fed_fine_mapping/`.

### 11.1 Measured agreement

`tests/test_fed_fine_mapping.py` runs both paths over one miniature three-site package (three
ancestries, 60 variants, a causal variant flanked by two near-perfect LD proxies so a credible
set has to spread its mass) and asserts the credible sets and PIPs agree. The residual
decomposes, and each half is pinned separately:

* **The summary-statistic half is exact.** Feeding the coordinator's closed-form (9.3) to
  SuSiEx alongside the comparator's *own* `plink --r` panels returns PIPs bit-identical to
  the centralized fit — asserted with no tolerance at all.
* **The LD half agrees to $\le 2$ units in the last place.** Both paths store $R_p$ as
  float32, so single precision is the floor on agreement by construction.

What a couple of float32 ULP in $R_p$ can move is the last of the six significant digits
SuSiEx prints (observed worst case 0.0729835 against 0.0729836), which is the tolerance the
PIP assertion uses.

On the production package, one instance (L0000, `ncsl1_h2-0.0005_rg1`, rep 0) was run through
both paths: the summary statistics are **byte-identical** to `plink2 --glm` for all 9,848
variants across the six ancestry columns, the credible set is the same single variant
`chr1:150522977:A:T`, and $\max_j|\Delta\text{PIP}_j| = 6.9\times10^{-16}$ over 2,033
variants — despite four of the six panels losing one or two variants each to the
condition-(ii) drop rule of §7. That is one instance at one locus, not the full grid.

---

## Summary of the federated claim

1. **Sufficiency (Lemma 1).** SuSiEx reads data only through $T_s = (R_s,\hat\beta_s,n_s,\tau_s^2)$.
2. **Exactness (Theorem 1).** Those statistics are additive over the site partition, so —
   given a harmonized variant set and a common complete-case / identically-imputed cohort
   (Theorem 1 conditions) — the coordinator's fit on aggregated statistics equals the
   centralized fit: federation is *lossless*, not approximate.
3. **Estimand-correctness (Corollary 1).** Stratifying by ancestry and pooling per-ancestry
   across sites reconstructs exactly the per-superpopulation GWAS that FedFM's ground truth is
   defined against — and, computed centrally, is the baseline the federated implementation is
   measured against (§8.1).
4. **One-round, genotype-private.** A single uplink of $O(M^2)$ aggregates per (site,
   ancestry); IBSS then runs coordinator-side with no further communication and no
   individual-level data ever leaving a site.
5. **Demonstrated, not just derived.** Both paths are implemented and run against each other
   (§11.1): identical credible sets, PIPs equal to the reference binary's print resolution,
   and a summary-statistic half that is bit-exact. The two places the claim is *not* exact on
   this package are harmonization failures of the kind Theorem 1's conditions name, not
   properties of the estimator — see "Conditions in practice" in §7.

The reason all of this is clean rather than hard-won: SuSiEx was designed from the start to
run on GWAS summary statistics, and summary statistics are exactly the additive sufficient
statistics of the Gaussian linear model. Federation is the corollary.

---

### References

- Yuan, K. *et al.* Fine-mapping across diverse ancestries drives the discovery of putative
  causal variants underlying human complex traits and diseases. *Nat. Genet.* **56**,
  1841–1850 (2024). [Model & Methods; PDF: `paper/nihms-2047735.pdf`.]
- Wang, G., Sarkar, A., Carbonetto, P. & Stephens, M. A simple new approach to variable
  selection in regression, with application to genetic fine mapping. *JRSS-B* **82**,
  1273–1300 (2020). [SuSiE / IBSS / SER.]
- Reference implementation: `getian107/SuSiEx` (v1.1.1). C++ core — `src/model.cpp` (IBSS,
  SER, ELBO, credible sets), `src/data.cpp` (summary-statistic ingestion, $\tau^2$,
  standardization).
- Wakefield, J. Bayes factors for genome-wide association studies (approximate Bayes factor
  underlying $\lambda_{sj}$ in §3.2).

---

# Appendix A — The math from scratch (beginner's walkthrough)

**Who this is for.** Someone who knows what a matrix is and roughly what a normal
distribution is, and wants to understand what SuSiEx computes and what changes when we
federate it. Nothing here is new — it is Parts I and II retold with the smallest possible
worked examples. Section numbers in parentheses point back to the rigorous statement.

---

## A.1 The raw materials: three objects

Everything starts with three things measured on a group of people.

**1. The genotype matrix $X$.** One row per person, one column per genetic variant (SNP).
Each entry is a **dosage**: how many copies of the alternate allele that person carries at
that SNP — 0, 1, or 2. With $n$ people and $M$ SNPs, $X$ is $n \times M$:

$$
X = \begin{pmatrix}
2 & 0 & 1 & 1\\
1 & 0 & 1 & 2\\
0 & 1 & 0 & 0\\
\vdots & & & \vdots
\end{pmatrix}
\quad
\begin{array}{l}
\leftarrow \text{person 1}\\
\leftarrow \text{person 2}\\
\leftarrow \text{person 3}
\end{array}
$$

**2. The phenotype vector $y$.** One number per person — the trait being studied (LDL
cholesterol, height, a disease liability score). Length $n$.

**3. The count $n$.** How many people. That is it. Sample size will turn out to be a
first-class input, not bookkeeping.

In FedFM (see [`design.md`](design.md)): 150,000 simulated people, split across 3 sites, and
a "locus" is a 1 Mb window of chromosome 1. $M$ is every variant in that window — thousands
of columns, not the $\approx$500 tag SNPs that [`design.md` §4.2](design.md) uses for the
separate job of measuring LD divergence between sites.

---

## A.2 The question, and why it is hard

The scientific question: **which SNP in this window actually causes the trait?**

The naive approach — **GWAS** — tests each SNP one at a time: regress $y$ on column $j$
alone, get a coefficient and a $p$-value, repeat $M$ times. This is fast and it works, in the
sense that the causal SNP will light up. The problem is that so will its neighbours.

The reason is **linkage disequilibrium (LD)**: nearby SNPs are inherited together, so their
columns in $X$ are *correlated*. If SNP 1 is causal and SNP 2 is 98% correlated with SNP 1,
then SNP 2 has essentially the same relationship with $y$ that SNP 1 does. A one-at-a-time
test cannot tell them apart. A typical GWAS hit is not one SNP; it is a smear of 50 SNPs with
near-identical $p$-values.

**Fine-mapping** is the job of taking that smear and asking which SNP is the real one — or,
honestly, which *small set* of SNPs the data cannot distinguish. The whole game is to use the
correlation structure explicitly instead of ignoring it.

That correlation structure is one matrix:

$$
R = \frac{1}{n} X^\top X \quad (M \times M), \qquad R_{jk} = \operatorname{corr}(\text{SNP } j,\ \text{SNP } k)
$$

when the columns of $X$ have been standardized (mean 0, variance 1). $R$ is the **LD matrix**.
Its diagonal is all 1s (every SNP correlates perfectly with itself) and its off-diagonals say
how confusable each pair of SNPs is.

---

## A.3 Compressing the data: from $(X, y)$ to $(R, \hat\beta, n)$

Here is the fact the entire method — and all of the federation — rests on.

**Standardize first.** Replace each genotype column by $\tilde x_j = (x_j - \bar x_j)/\mathrm{sd}(x_j)$,
and likewise $\tilde y$. Now every column has mean 0 and variance 1, and dosages are on a
common scale regardless of how common the allele is.

**Then two summaries are all you need:**

$$
\hat\beta = \frac{1}{n}\tilde X^\top \tilde y \quad (\text{length } M), \qquad
R = \frac{1}{n}\tilde X^\top \tilde X \quad (M \times M).
$$

- $\hat\beta_j$ is exactly the **correlation between SNP $j$ and the trait** — the
  standardized GWAS effect estimate for that SNP. One number per SNP.
- $R$ is the LD matrix above. One number per *pair* of SNPs.

These are related to the $z$-scores you see in GWAS summary files by

$$
z_j = \hat\beta_j \sqrt{n} \quad\Longleftrightarrow\quad \hat\beta_j = \frac{z_j}{\sqrt n},
$$

which is literally the line the reference implementation runs on ingestion
(`beta = (effect/se)/√n`, `data.cpp:311`).

**Why this matters.** $X$ has $n \times M$ = 36,500 × 2,000 ≈ 73 million entries and is
personally identifying. $(\hat\beta, R, n)$ has $M + M^2$ ≈ 4 million entries and is
per-population aggregate. SuSiEx never touches $X$ or $y$ again — it reads only
$(\hat\beta, R, n)$. That is Lemma 1 (§6), and it is the reason federation is going to be
easy rather than hard.

---

## A.4 One causal SNP: the single-effect regression

Start with a deliberately simple model: **exactly one SNP in this window is causal, and we
don't know which.** Let $\gamma$ be a one-hot vector of length $M$ marking the causal one,
and let its effect size be drawn from a normal prior, $b \sim \mathcal N(0, \tau^2)$.

Bayes' rule then gives a posterior over *which* SNP it is. For each candidate $j$ we need one
number: how much better the data look if $j$ is the causal SNP than if nothing is. That number
is the **log Bayes factor**, and for this model it has a closed form (the Wakefield ABF, §3.2):

$$
\lambda_j \;=\; \underbrace{\tfrac12 \log\frac{v^2}{v^2 + \tau^2}}_{\text{(a) complexity penalty}}
\;+\; \underbrace{\tfrac12\, z_j^2\, \frac{\tau^2}{v^2 + \tau^2}}_{\text{(b) signal strength}},
\qquad v^2 = \frac{\sigma^2}{n},\quad z_j = \frac{\hat\beta_j}{v}.
$$

Read it as: **(b)** rewards a big $z$-score — strong marginal association; **(a)** is a fixed
cost paid by every SNP for the privilege of having a free parameter, so a SNP has to earn its
keep. $v^2 = \sigma^2/n$ is the sampling variance of $\hat\beta_j$ — this is where $n$ enters,
and why more people means sharper answers.

**Worked example.** Take EUR with $n = 36{,}500$, residual variance $\sigma^2 = 1$, so
$v^2 = 1/36{,}500 = 2.74\times10^{-5}$ and $v = 0.00523$. Suppose the strongest SNP in the
window has $\hat\beta = 0.03$, which sets the prior variance $\tau^2 = 0.03^2 = 9\times10^{-4}$
(§3.5). Then the shrinkage factor is $\tau^2/(v^2+\tau^2) = 0.970$ and the penalty term (a) is
$\tfrac12\log(0.0295) = -1.76$. So:

| SNP | $\hat\beta$ | $z = \hat\beta/v$ | (a) penalty | (b) signal | $\lambda$ |
|---|---:|---:|---:|---:|---:|
| a strong hit | 0.030 | 5.73 | $-1.76$ | $+15.94$ | $\mathbf{+14.18}$ |
| a null SNP | 0.003 | 0.57 | $-1.76$ | $+0.16$ | $\mathbf{-1.60}$ |

The strong SNP has a log Bayes factor of $+14$ (the data are $e^{14} \approx 1.2$ million times
more likely under "this SNP is causal"); the null SNP scores *negative* — the penalty exceeds
the evidence. Good.

**Turning scores into probabilities.** Normalize with a softmax, weighted by the prior $\pi_j$
that SNP $j$ is causal (uniform $1/M$ by default):

$$
\alpha_j = \frac{\pi_j\, e^{\lambda_j}}{\sum_k \pi_k\, e^{\lambda_k}}.
$$

$\alpha_j$ is the **posterior inclusion probability (PIP)** — "probability SNP $j$ is the
causal one". The $\alpha$ vector sums to 1 by construction.

**Worked example, four SNPs in tight LD.** Say the log-BFs come out

$$
\lambda = (14.2,\; 13.9,\; 10.2,\; 9.4).
$$

Only *differences* matter (the denominator cancels the rest), so subtract the max:
$(0, -0.3, -4.0, -4.8)$, exponentiate to $(1,\ 0.741,\ 0.0183,\ 0.00823)$, divide by the sum
$1.767$:

$$
\alpha = (0.566,\; 0.419,\; 0.010,\; 0.005).
$$

SNP 1 is the best guess at 57% — but SNP 2 is right behind it at 42%, because it is in near-perfect
LD with SNP 1 and the data genuinely cannot separate them. Hold onto this example; §A.7 fixes it.

---

## A.5 Several causal SNPs: the IBSS loop

Real loci often have more than one causal variant. SuSiE's trick is not to fit a complicated
multi-variant model but to **stack $L$ copies of the simple one** (typically $L = 5$ or 10):

$$
\beta = b_1 + b_2 + \dots + b_L, \qquad \text{each } b_l \text{ nonzero at exactly one SNP}.
$$

Each $b_l$ is called a **single effect**. Fitting is done by **IBSS** (Iterative Bayesian
Stepwise Selection) — a loop that is close to backfitting:

```
initialise all L effects to zero
repeat until the ELBO stops changing (usually <10 passes):
    for l = 1 .. L:
        # 1. RESIDUALIZE: pretend the other L-1 effects are known and subtract them
        r = β̂ − R · (sum of all effects except l)         # §3.1

        # 2. SOLVE: run the one-causal-SNP analysis of A.4 on the residual
        λ, α, μ = single_effect_regression(r, τ², σ², n)   # §3.2

        # 3. PUT BACK: effect l is now the PIP-weighted estimate
        b_l = α ⊙ μ                                       # §3.3
    re-estimate the residual variance σ²                   # §3.4
```

Step 1 is the only line that needs comment. In individual-level regression, "subtract the
other effects" means $y - X\sum_{l'\neq l} b_{l'}$. In summary-statistic space the same
operation is $\hat\beta - R \sum_{l'\neq l} b_{l'}$ — the LD matrix $R$ plays the role of
$\tfrac1n X^\top X$. This is why the whole loop runs without $X$: **$R$ is the only thing the
algorithm ever needed $X$ for.**

If the locus really has only two signals and you set $L = 5$, the extra three effects just
spread their $\alpha$ mass thinly over everything — they are detected and dropped in
post-processing (§A.6).

---

## A.6 Reading the output

IBSS converges and hands back $L$ vectors $\alpha_1, \dots, \alpha_L$, each a probability
distribution over the $M$ SNPs. Three things are derived from them (§4):

**Credible sets.** For effect $l$, sort SNPs by $\alpha_{lj}$ descending and take them until
the cumulative probability crosses 95%. That set is the **95% credible set**: "we are 95% sure
the causal SNP for this signal is one of these." A credible set of size 1 is a resolved
variant; size 40 means the locus is unresolvable with this data.

**Overall PIP.** A SNP's total probability of being causal for *any* of the signals:
$\mathrm{PIP}_j = 1 - \prod_l (1 - \alpha_{lj})$ — i.e. 1 minus the probability that every
effect passed it over.

**Purity filter.** A junk credible set (from an unneeded $l$) contains SNPs scattered across
the window with no LD between them. So compute the set's **purity** = the smallest absolute
LD among its members, and throw the set away if purity < 0.5. A genuine credible set is a
tight LD clump; a diffuse one is an artefact of $L$ being too large.

---

## A.7 Many ancestries: why SuSiE**x** exists

Here is the payoff, and the reason the method has an "x" (cross-population) on the end.

LD is **ancestry-specific**. Two SNPs 98% correlated in Europeans may be only 40% correlated
in Africans, because African populations are older and their LD blocks are shorter. The pair of
SNPs that is *indistinguishable* in one ancestry can be *easily separated* in another.

SuSiEx exploits this with one modelling assumption (§2):

> **The causal SNP is the same in every ancestry; its effect size may differ.**

Formally, the one-hot indicator $\gamma_l$ is shared across populations, while the magnitude
$b_{sl}$ is population-specific. So the algorithm runs the single-effect regression of §A.4
*separately per ancestry* — each with its own $\hat\beta_s$, $R_s$, $n_s$, $\tau_s^2$ — and
then combines the evidence in the only way independent evidence can be combined: **multiply the
Bayes factors, i.e. add the logs.**

$$
\Lambda_j = \sum_{s=1}^{S} \lambda_{sj}, \qquad
\alpha_j = \frac{\pi_j e^{\Lambda_j}}{\sum_k \pi_k e^{\Lambda_k}}.
$$

**Worked example, continued.** From §A.4, EUR alone gave

$$
\lambda^{\mathrm{EUR}} = (14.2,\; 13.9,\; 10.2,\; 9.4)
\;\Longrightarrow\;
\alpha = (0.566,\; 0.419,\; 0.010,\; 0.005),
$$

a 95% credible set of **{SNP 1, SNP 2}** — size 2, unresolved. Now add AFR, where SNP 1 and
SNP 2 sit in different LD blocks so SNP 2's marginal association is much weaker:

$$
\lambda^{\mathrm{AFR}} = (13.0,\; 8.0,\; 9.5,\; 9.1).
$$

Add them SNP-by-SNP:

$$
\Lambda = (27.2,\; 21.9,\; 19.7,\; 18.5).
$$

Subtract the max: $(0, -5.3, -7.5, -8.7)$; exponentiate: $(1,\ 0.00499,\ 0.00055,\ 0.00017)$;
normalize by $1.0057$:

$$
\alpha = (\mathbf{0.994},\; 0.005,\; 0.0005,\; 0.0002).
$$

The 95% credible set is now **{SNP 1}** — size 1. *This is the entire value proposition of
multi-ancestry fine-mapping*, and it is why the FedFM simulator stratifies loci by inter-site
LD divergence ([`design.md` §4.3](design.md)): the benefit is largest exactly where the LD
matrices disagree most.

Note the division of labour, which will matter in a moment: $\alpha_j$ is **one shared vector**
across all ancestries (there is one causal SNP), while the effect sizes $\mu_{sj}$ stay
**per-ancestry** (the effect may be bigger in one group than another, or zero).

---

## A.8 The federation problem

Now the constraint that makes this a federated-learning problem. In FedFM the 150,000 people
do not live in one place. They live at $K = 3$ institutions (ANL, Covenant, MBZUAI), and
data-residency rules mean **no genotype row may leave the institution that holds it**. Nobody
is allowed to build the pooled $X$.

Two distinct groupings are now in play, and keeping them apart is essential:

| | what it is | count |
|---|---|---|
| **site** $k$ | a physical institution holding a machine with data on it | $K = 3$ |
| **population** $s$ | an ancestry (EUR, AFR, AMR, EAS, CSA, MID) | $S = 6$ |

Every site holds a **mixture** of ancestries, and every ancestry is **spread across** sites.
SuSiEx's columns are ancestries, not sites — an AFR person at ANL and an AFR person at Covenant
share the same causal effect by construction ([`design.md` §5.2](design.md)), so they belong in
the same column. This gives a 3 × 6 grid of (site, ancestry) blocks, e.g. AFR = 7,500 at ANL +
47,500 at Covenant + 5,000 at MBZUAI = 60,000 people (§8.1). Our job is to fit SuSiEx over the
6 ancestry columns without ever assembling any of them in one place.

---

## A.9 Why it's easy: everything is a sum over people

Recall from §A.3 that SuSiEx only ever reads $(\hat\beta_s, R_s, n_s)$, and those come from
$X^\top X$ and $X^\top y$. Look at what $X^\top X$ actually is:

$$
X^\top X = \sum_{i=1}^{n} x_i x_i^\top
$$

— a **sum with one term per person**. And a sum does not care how you group its terms. Split
the people into site A and site B and you can add up each group separately, then add the two
answers:

$$
X^\top X = \underbrace{\sum_{i \in A} x_i x_i^\top}_{G_A \text{, computed at site A}} + \underbrace{\sum_{i \in B} x_i x_i^\top}_{G_B \text{, computed at site B}}.
$$

That is Theorem 1 (§7), and it is just associativity of addition. The same holds for
$X^\top y$, for the column sums $\mathbf 1^\top X$, and for $n$ itself.

**Worked example.** Six people, two SNPs $j$ and $k$, split 3/3 across two sites:

| | person | $x_j$ | $x_k$ |
|---|---|---:|---:|
| site A | 1 | 2 | 2 |
| | 2 | 2 | 1 |
| | 3 | 0 | 0 |
| site B | 4 | 0 | 0 |
| | 5 | 0 | 1 |
| | 6 | 2 | 2 |

Site A computes, from its own rows only:

$$
G_A = X_A^\top X_A = \begin{pmatrix} 8 & 6\\ 6 & 5\end{pmatrix}, \qquad
u_A = \mathbf 1^\top X_A = (4,\ 3), \qquad n_A = 3
$$

(e.g. $8 = 2^2+2^2+0^2$, $6 = 2\!\cdot\!2 + 2\!\cdot\!1 + 0\!\cdot\!0$). Site B computes

$$
G_B = \begin{pmatrix} 4 & 4\\ 4 & 5\end{pmatrix}, \qquad u_B = (2,\ 3), \qquad n_B = 3.
$$

The coordinator adds them:

$$
G = G_A + G_B = \begin{pmatrix} 12 & 10\\ 10 & 10\end{pmatrix}, \qquad u = (6,\ 6), \qquad n = 6,
$$

which is **exactly** what you would get from the pooled 6-row matrix — check
$12 = 4+4+0+0+0+4$. No approximation, no averaging, no iteration. The coordinator now has
everything it needs and no genotype row ever moved.

**One round, then silence.** Because the aggregates are static — they do not depend on the
current parameter values — the sites transmit **once**, and the entire IBSS loop of §A.5 then
runs on the coordinator with **zero** further communication. This is what makes SuSiEx unusual
among federated methods: federated SGD or federated EM need a synchronization round per
iteration; this needs one round, ever.

---

## A.10 The one real trap: standardize *last*

There is exactly one place to get this wrong, and it is worth doing the arithmetic (§9).

§A.3 said to standardize the columns before computing $R$ and $\hat\beta$. The tempting
implementation is: each site standardizes its own columns, computes its own LD matrix $R_k$,
and ships that. **This is biased**, because each site standardizes against *its own* mean and
SD — the blocks end up on different scales and no longer describe a common population.

**Worked example.** Same six people as §A.9. Note that SNP $j$'s allele frequency differs
between the sites (site A mean $4/3$, site B mean $2/3$) while SNP $k$'s does not (both mean 1).

*The right answer* — pool the raw sums, then standardize once:

$$
\bar x_j = \tfrac{6}{6} = 1, \quad \operatorname{var}(x_j) = \tfrac{12}{6} - 1^2 = 1, \quad
\bar x_k = 1, \quad \operatorname{var}(x_k) = \tfrac{10}{6} - 1^2 = 0.667,
$$
$$
\operatorname{cov}(x_j, x_k) = \tfrac{10}{6} - 1\cdot 1 = 0.667
\quad\Longrightarrow\quad
R_{jk} = \frac{0.667}{\sqrt{1}\cdot\sqrt{0.667}} = \mathbf{0.816}.
$$

*The wrong answer* — standardize locally, then combine. Site A gets
$\operatorname{corr}_A(j,k) = 0.667/(0.943 \times 0.816) = 0.866$; site B, by symmetry, also
$0.866$. Any weighted average of the two is $\mathbf{0.866}$.

$0.866 \neq 0.816$, and the gap is not floating-point noise. The mechanism is visible in the
numbers: pooling adds the *between-site* variance in SNP $j$ (its frequency really does differ
between sites) to the denominator, but adds nothing to the numerator, since SNP $k$ has no
between-site difference to co-vary with. Local centering deletes exactly that between-site
information before the coordinator can see it.

**The fix.** Sites ship **raw, unstandardized** sums; the coordinator standardizes once
against the *pooled* moments. Per (site, ancestry) block, a site emits six things:

| symbol | what it is | size |
|---|---|---|
| $G_{p,k} = X^\top X$ | genotype cross-products | $M \times M$ |
| $c_{p,k} = X^\top y$ | genotype–phenotype cross-products | $M$ |
| $u_{p,k} = \mathbf 1^\top X$ | column sums (for the means) | $M$ |
| $n_{p,k}$ | number of people | 1 |
| $q_{p,k} = \mathbf 1^\top y$ | phenotype sum | 1 |
| $w_{p,k} = y^\top y$ | phenotype sum of squares | 1 |

and the coordinator adds them across sites and applies equations (9.1)–(9.2) to get
$R_p$ and $\hat\beta_p$ exactly. The prior variance $\tau_p^2 = \max_j \hat\beta_{pj}^2$ is
also computed after pooling, so sites never need to agree on it in advance.

**Two smaller traps in the same family:**

1. **Never sum LD matrices.** $R_A + R_B$ has 2s on the diagonal, which is not a correlation
   matrix. If you are given per-site $R_k$ (rather than Grams), the correct recombination is
   the $n$-weighted *average* $R = \frac{1}{n}\sum_k n_k R_k$ — and even that is only exact if
   every site standardized against a common external reference.
2. **Allele flips.** If site A codes a SNP relative to allele `A` and site B relative to `G`,
   site B's column is $2 - x$ rather than $x$: its sign flips after centering, corrupting every
   off-diagonal it touches. Harmonizing the variant list and reference alleles across sites is
   condition (i) of Theorem 1 and is a real engineering task, not a formality.

   **This one is not hypothetical here.** FedFM's own per-site filesets were cut without
   `plink --keep-allele-order`, so PLINK picked A1 by each site's own minor allele and ~6% of
   chr1 variants are coded oppositely at ANL and Covenant (§7, "Conditions in practice"). The
   centralized baseline never notices, because it never opens a site fileset. The federated
   path had to recode against a published canonical variant list before touching a single
   moment — which is exactly what "a real engineering task" means in practice.

---

## A.11 The whole thing in one page

```
AT EACH SITE k  (data never leaves)
  ├─ split my people by ancestry p using {site}_manifest.tsv
  └─ for each ancestry p I hold:
        read my genotype block X for the locus window (M variants)
        emit G[p,k]=XᵀX,  c[p,k]=Xᵀy,  u[p,k]=1ᵀX,  n[p,k],  q[p,k],  w[p,k]
                                    │
                                    │  ONE transmission, O(M²) floats
                                    ▼
AT THE COORDINATOR
  ├─ for each ancestry p:  G[p] = Σ_k G[p,k],  c[p] = Σ_k c[p,k],  n[p] = Σ_k n[p,k], …
  ├─ standardize ONCE against pooled moments  →  R[p] (M×M),  β̂[p] (M),  n[p],  τ²[p]
  └─ run IBSS over the S=6 ancestry columns:
        repeat until ELBO converges:
          for l = 1..L:
            residualize:  r[p] = β̂[p] − R[p]·(other effects)        per ancestry
            score:        λ[p,j] = Wakefield log-BF                  per ancestry, per SNP
            COUPLE:       Λ[j] = Σ_p λ[p,j]                          ← ancestries meet here
            PIP:          α[j] = softmax(Λ)[j]                       shared across ancestries
            update:       b[p,l] = α ⊙ μ[p]                          per ancestry
  └─ credible sets (95% mass) → purity filter → PIPs
```

Cost: one round of communication; $M^2$ float64s per (site, ancestry) block, which is 32 MB at
$M = 2{,}000$. What crosses the boundary is an LD matrix and a set of marginal
associations — the same artifacts any public GWAS summary-statistic release already contains
(§10).

And the result is not an approximation of the centralized answer — by Theorem 1 it *is* the
centralized answer, bit-for-bit up to summation order. That is what makes
[`src/fine_mapping.py`](../../../../src/appfl_bio_suite/experiments/fine_mapping/fedfm/fine_mapping.py)'s centralized baseline (§8.1) a legitimate
reference: it computes the same estimator by pooling directly, so any discrepancy the federated
implementation shows is an implementation bug or a violated harmonization condition — never an
inherent cost of federating.

The box above is [`src/fed_fine_mapping.py`](../../../../src/appfl_bio_suite/experiments/fine_mapping/fedfm/fed_fine_mapping.py), roughly function for
function, and the two paths have been run against each other: same credible sets, PIPs equal
to the last digit SuSiEx prints, and a bit-exact summary-statistic half (§11.1). Both
discrepancies that did turn up were harmonization conditions — the allele flips of §A.10 and
the missing-genotype handling of Theorem 1 condition (ii) — never the federation itself.

---

## A.12 Symbol glossary

| symbol | plain English | where it lives |
|---|---|---|
| $n$ | number of people | a scalar per (ancestry, site) |
| $M$ | number of SNPs in the locus window | ~thousands, 1 Mb of chr1 |
| $X$ | genotype matrix, dosages 0/1/2 | **never leaves a site** |
| $y$ | the trait value per person | **never leaves a site** |
| $\hat\beta_j$ | GWAS effect for SNP $j$ = corr(SNP $j$, trait) | coordinator, per ancestry |
| $z_j$ | $\hat\beta_j\sqrt n$, the familiar GWAS $z$-score | derived |
| $R$ | LD matrix: correlation between every pair of SNPs | coordinator, per ancestry |
| $\tau^2$ | prior guess at how big a causal effect is | $\max_j \hat\beta_j^2$ |
| $\sigma^2$ | residual (unexplained) variance | re-estimated each IBSS pass |
| $L$ | how many causal SNPs we allow | 5, or 10 if 5 isn't enough |
| $\lambda_{sj}$ | evidence (log Bayes factor) for SNP $j$ in ancestry $s$ | per ancestry |
| $\Lambda_j$ | total evidence for SNP $j$ = $\sum_s \lambda_{sj}$ | **the coupling step** |
| $\alpha_{lj}$ | PIP: P(SNP $j$ is causal for effect $l$) | shared across ancestries |
| $\gamma_l$ | which SNP effect $l$ acts through (one-hot) | shared across ancestries |
| $b_{sl}$ | how big that effect is in ancestry $s$ | per ancestry |
| $G, c, u, q, w$ | the raw sums a site ships | site → coordinator, once |
| $k$ vs $s$ | site (3 institutions) vs population (6 ancestries) | do not confuse these |
</content>
