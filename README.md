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
  -v genotypes.vcf.gz \
  -b phenotypes.bed.gz \
  -c covariates.txt \
  -r 22:17000000-18000000 \
  -w 1000000 \
  --threshold 0.05 \
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
  -v genotypes.vcf.gz \
  -b phenotypes.bed.gz \
  -c covariates.txt \
  -r 22:17000000-18000000 \
  -p 100 1000 \
  --seed 42 \
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
| Nominal p-value | Exact F(1, df) CDF via Rmath `pf()` | Normal-tail approximation `2·Φ(−|t|)` — anti-conservative for small df |
| Permutation shuffle | Shuffles raw (pre-residualization) phenotype, then re-residualizes each permutation | Shuffles already-residualized phenotype; no re-residualization per permutation |
| Covariate residualization | Eigen `ColPivHouseholderQR` (QR decomposition) | Normal equations (X'X)⁻¹X'y via Gaussian elimination — less numerically stable for ill-conditioned covariate matrices |
| PRNG | `std::random_shuffle` (stdlib, system-seeded) | XorShift64 — seeds are incompatible |
| `ref_factor` output | Reported (REF/ALT orientation) | Not reported |
| Conditional mapping | Implemented (`--mapping`) | Not implemented |
| Interaction testing | Implemented (`--interaction`) | Not implemented |
| Group-aware permutation | Implemented (`--grp`) | Not implemented |
| Chunk/parallelisation | Via Python wrapper | Not implemented (run per-region) |

---

## Dependencies

- [`clap`](https://crates.io/crates/clap) 4.x — argument parsing

No other external crates. Math functions (lgamma, regularized incomplete beta, Nelder-Mead) are implemented from scratch.
