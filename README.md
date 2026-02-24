# rust-fastqtl

A Rust reimplementation of the core FastQTL cis-QTL mapping logic.

**Original C++ source:** https://github.com/francois-a/fastqtl

**Reimplemented by:** GPT-5.3-codex and Claude-sonnet-4.6

> **Warning:** The numerical outputs of `rust-fastqtl` have not been formally evaluated against the original FastQTL. Results may differ due to implementation differences described in the [Known differences](#known-differences) section below. Use with caution and validate before drawing conclusions.

---

## Features

- Parse phenotype BED (`.bed` / `.bed.gz`) within a target `--region`
- Parse genotype VCF (`.vcf` / `.vcf.gz`) within the cis-window (`--window`)
- MAF and minor-allele-sample filtering (`--maf-threshold`, `--ma-sample-threshold`)
- Missing genotype/phenotype imputation by per-feature mean
- Optional phenotype rank-normal transformation (`--normal`)
- Optional covariate correction by OLS residualization (`--cov`)
- Mean-centering + L2 normalization before correlation computation
- **Nominal mode** (default): scan all cis variants per phenotype
- **Permutation mode** (`--permute`): empirical p-value with adaptive stopping, Beta distribution fitting (`beta_shape1`, `beta_shape2`), effective degrees-of-freedom learning, and beta-adjusted p-value — matching FastQTL's permutation pass output

---

## Installation

```bash
git clone https://github.com/huangnengCSU/rust-fastqtl.git
cd rust-fastqtl
cargo build --release
# binary at target/release/rust_fastqtl
```

Requires Rust 1.85+ (edition 2024).

---

## Usage

```
Usage: rust_fastqtl [OPTIONS] --vcf <VCF> --bed <BED> --out <OUT> --region <REGION>

Options:
  -v, --vcf <VCF>                            Input VCF/BCF file (may be gzip-compressed)
  -b, --bed <BED>                            Input BED phenotype file (may be gzip-compressed)
  -o, --out <OUT>                            Output file path
  -c, --cov <COV>                            Covariate file (optional)
  -r, --region <REGION>                      Genomic region to analyse (chr:start-end or chr)
  -w, --window <WINDOW>                      Cis-window size in bp [default: 1000000]
      --threshold <THRESHOLD>                P-value threshold for nominal-pass output [default: 1]
      --maf-threshold <MAF_THRESHOLD>        Minor allele frequency filter [default: 0]
      --ma-sample-threshold <N>              Minimum samples carrying the minor allele [default: 0]
  -p, --permute <PERMUTE>...                 Permutation counts (1–3 integers, see below)
      --seed <SEED>                          Random seed [default: 12345]
  -n, --normal                               Apply rank-normal transformation to phenotypes
  -h, --help                                 Print help
```

### `--permute` modes

| Invocation | Behaviour |
|---|---|
| `--permute N` | Fixed N permutations |
| `--permute H M` | Adaptive: stop when `n_better ≥ H` or `n_perms ≥ M` |
| `--permute N H M` | Adaptive with guard: stop when `n_perms ≥ N` and (`n_better ≥ H` or `n_perms ≥ M`) |

### Nominal mode

```bash
rust_fastqtl \
  -v example/genotypes.vcf.gz \
  -b example/phenotypes.bed.gz \
  -c example/covariates.txt.gz \
  -r 22:17000000-18000000 \
  -w 1000000 \
  --threshold 1.0 \
  --maf-threshold 0.01 \
  --ma-sample-threshold 1 \
  -o nominal.txt
```

Output columns:
```
phenotype_id  variant_id  distance  ma_samples  ma_count  maf  pval  beta  beta_se
```

### Permutation mode

```bash
rust_fastqtl \
  -v example/genotypes.vcf.gz \
  -b example/phenotypes.bed.gz \
  -c example/covariates.txt.gz \
  -r 22:17000000-18000000 \
  -p 100 1000 \
  --seed 12345 \
  -o permutation.txt
```

Output columns:
```
phenotype_id  n_variants  beta_shape1  beta_shape2  eff_df  p_nom_eff
variant_id  distance  ma_samples  ma_count  maf  p_nom  beta  beta_se
p_empirical  p_beta_adjusted
```

`p_beta_adjusted` is the key column for downstream FDR control (e.g. Storey's q-value or Benjamini-Hochberg).

---

## Known differences from original FastQTL

The following differences from https://github.com/francois-a/fastqtl are known and have **not** been benchmarked:

| Aspect | Original FastQTL (C++) | rust-fastqtl |
|---|---|---|
| Nominal p-value | Exact F(1, df) CDF via GSL `gsl_cdf_fdist_Q()` | Exact F(1, df) CDF via custom regularized incomplete beta (`pbeta`) — numerically equivalent |
| Internal numeric precision | `float` (32-bit) for genotype/phenotype storage; `double` for statistics | `f64` (64-bit) throughout — small numerical differences possible |
| Permutation shuffle | Shuffles raw (pre-residualization) phenotype, then re-residualizes each permutation | Shuffles already-residualized phenotype; no re-residualization per permutation |
| Covariate residualization | Eigen `ColPivHouseholderQR` (QR decomposition); detects and drops linearly dependent covariates | Normal equations (X'X)⁻¹X'y via Gauss-Jordan inversion — less numerically stable for ill-conditioned matrices; no rank-deficiency handling |
| Categorical covariates | Automatically converts factor covariates to binary dummy variables; drops constant covariates | Only numeric covariates supported |
| `learnDF` optimizer | 1-D Nelder-Mead via GSL (`nmsimplex2`), ≤20 iterations, tol=0.01 | Golden-section search over [1, max(3×df, 100)], 50 iterations, tol=0.01 |
| PRNG | `std::random_shuffle` backed by `rand()`; supports `--seed` (default: `time(NULL)`) | XorShift64; supports `--seed` (default: 12345) — algorithms and seed spaces are incompatible |
| `ref_factor` output | Reported in permutation output (REF/ALT orientation) | Not reported |
| Conditional mapping | Implemented (`--mapping`) | Not implemented |
| Interaction testing | Implemented (`--interaction`) | Not implemented |
| Group-aware permutation | Implemented (`--grp`) | Not implemented |
| Chunk/parallelisation | Via Python wrapper | Not implemented (run per-region) |

---

## Dependencies

- [`clap`](https://crates.io/crates/clap) 4.x — argument parsing

No other external crates. Math functions (lgamma, regularized incomplete beta, Nelder-Mead) are implemented from scratch.
