# rust-fastqtl

A Rust reimplementation of the core FastQTL cis-QTL mapping logic.

**Original C++ source:** https://github.com/francois-a/fastqtl

**Reimplemented by:** GPT-5.3-codex and Claude-sonnet-4.6

> **Warning:** The numerical outputs of `rust-fastqtl` have not been formally evaluated against the original FastQTL. Results may differ due to implementation differences described in the [Known differences](#known-differences) section below. Use with caution and validate before drawing conclusions.

---

## Features

- Parse phenotype BED (`.bed` / `.bed.gz`) — loaded once, filtered per region
- Parse genotype VCF (`.vcf` / `.vcf.gz`) within the cis-window (`--window`)
- Stream gzip input (no full-file decompression into memory)
- Use `tabix -h` region retrieval automatically for `.vcf.gz` with `.tbi/.csi` index
- **Multi-region parallel processing**: pass multiple `-r` regions; all regions run concurrently with `rayon` and results are merged in order into a single output file; thread count controlled by `--threads`
- Parse methylation matrix BED (`--bedmethyl`) — rows = CpG sites, columns 4+ = per-sample integer dosage (0/1/2); methylation used as **predictor** (genotype-like), phenotype BED as **outcome** (e.g. splicing → methyl-sQTL; expression → methyl-eQTL)
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
Usage: rust_fastqtl [OPTIONS] --bed <BED> --out <OUT> -r <REGION>...
                    (--vcf <VCF> | --bedmethyl <BEDMETHYL>)

Options:
  -v, --vcf <VCF>                            Input VCF/BCF file (may be gzip-compressed)
  -m, --bedmethyl <BEDMETHYL>                Input methylation matrix BED file (rows=methylation sites, cols 4+=per-sample integer dosage 0/1/2)
  -b, --bed <BED>                            Input BED phenotype file (may be gzip-compressed)
  -o, --out <OUT>                            Merged output file (all regions combined)
      --out-dir <DIR>                        Output directory for per-region files (optional)
  -c, --cov <COV>                            Covariate file (optional)
  -r, --region <REGION>...                   Region(s) to analyse (chr:start-end or chr); repeat for multiple
  -t, --threads <THREADS>                    Parallel threads [default: all CPUs]
  -w, --window <WINDOW>                      Cis-window size in bp [default: 1000000]
      --threshold <THRESHOLD>                P-value threshold for nominal-pass output [default: 1]
      --maf-threshold <MAF_THRESHOLD>        Minor allele frequency filter [default: 0]
      --ma-sample-threshold <N>              Minimum samples carrying the minor allele [default: 0]
  -p, --permute <PERMUTE>...                 Permutation counts (1–3 integers, see below)
      --seed <SEED>                          Random seed [default: 12345]
  -n, --normal                               Apply rank-normal transformation to phenotypes
  -h, --help                                 Print help
```

Exactly one of `--vcf` or `--bedmethyl` must be provided.

### VCF I/O behavior

- For compressed VCF input (`.vcf.gz`), parsing is now streamed.
- If a tabix index (`.tbi` or `.csi`) is present, `rust_fastqtl` automatically uses `tabix -h <vcf> <region>` to fetch only the requested cis-window region.
- If tabix is unavailable (or no index is found), it falls back to sequential scan.

### `--permute` modes

| Invocation | Behaviour |
|---|---|
| `--permute N` | Fixed N permutations |
| `--permute H M` | Adaptive: stop when `n_better ≥ H` or `n_perms ≥ M` |
| `--permute N H M` | Adaptive with guard: stop when `n_perms ≥ N` and (`n_better ≥ H` or `n_perms ≥ M`) |

### Nominal mode

Single region:
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

Multiple regions (all chromosomes, parallel) — merged output only:
```bash
rust_fastqtl \
  -v example/genotypes.vcf.gz \
  -b example/phenotypes.bed.gz \
  -c example/covariates.txt.gz \
  -r chr1 chr2 chr3 chr4 chr5 \
  -t 5 \
  --maf-threshold 0.01 \
  --ma-sample-threshold 5 \
  -o nominal_all.txt
```

Multiple regions — merged output **and** per-region files:
```bash
rust_fastqtl \
  -v example/genotypes.vcf.gz \
  -b example/phenotypes.bed.gz \
  -c example/covariates.txt.gz \
  -r chr1 chr2 chr3 chr4 chr5 \
  -t 5 \
  --maf-threshold 0.01 \
  --ma-sample-threshold 5 \
  -o nominal_all.txt \
  --out-dir results/
# produces: nominal_all.txt  +  results/chr1.txt  results/chr2.txt  ...
```

> `--region` / `-r` accepts any mix of whole-chromosome (`chr1`) and sub-region (`chr1:1000000-2000000`) strings. Each region is processed in a separate rayon thread; phenotype BED is loaded once and filtered per region; genotype VCF is queried independently per region via tabix. `--out` is always written (merged); `--out-dir` is optional and produces one file per region named `<region>.txt` (`:` replaced by `_` in filenames).

With bedMethyl (methyl-sQTL / methyl-eQTL — methylation as predictor):
```bash
rust_fastqtl \
  -m example/methylation.bed.gz \
  -b example/phenotypes.bed.gz \
  -c example/covariates.txt.gz \
  -r 22:17000000-18000000 \
  -w 1000000 \
  --threshold 1.0 \
  -o nominal.txt
```

Output columns:
```
phenotype_id  variant_id  distance  distance_to_body  ma_samples  ma_count  maf  pval  beta  beta_se
```

| Column | Description |
|---|---|
| `distance` | `variant_pos − phenotype_start` (signed; TSS-anchored, matches original FastQTL) |
| `distance_to_body` | Signed distance anchored at the nearer of `phenotype_start` or `phenotype_end`. The closer boundary is chosen as anchor; the value is `variant_pos − anchor` (negative = upstream of that boundary, positive = downstream). More informative than `distance` for long phenotypes such as introns, where the variant may be far from the TSS but close to the splice site. |

### Permutation mode

With VCF genotypes:
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

With bedMethyl (methyl-sQTL / methyl-eQTL):
```bash
rust_fastqtl \
  -m example/methylation.bed.gz \
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
variant_id  distance  distance_to_body  ma_samples  ma_count  maf  p_nom  beta  beta_se
p_empirical  p_beta_adjusted
```

`distance` and `distance_to_body` follow the same definitions as in the nominal pass (see above).

`p_beta_adjusted` is the key column for downstream FDR control (e.g. Storey's q-value or Benjamini-Hochberg).

---

## bedMethyl format (`--bedmethyl`)

A tab-separated methylation matrix file where rows are CpG (or other modified base) sites and columns 4+ are per-sample integer dosage values. Methylation is used as the **predictor** (analogous to SNP genotype), and the phenotype BED (`--bed`) is the **outcome**. Depending on the phenotype this enables:

- **methyl-sQTL**: methylation → splicing
- **methyl-eQTL**: methylation → gene expression
- other phenotype types as appropriate

> Note: standard **mQTL** (genetic variant → methylation) uses `--vcf` with a methylation phenotype BED, not `--bedmethyl`.

Only **integer dosage encoding (0/1/2)** is supported:

```
#chr    start   end     site_id         SAMPLE1  SAMPLE2  SAMPLE3  ...
chr22   100000  100001  CpG_100001      2        1        2        ...
chr22   200000  200001  CpG_200001      0        0        1        ...
```

Each value is an integer genotype-like call:

| Value | Meaning |
|-------|---------|
| `0`   | Low methylation |
| `1`   | Hemi-methylation |
| `2`   | High methylation |

MAF and `ma_count` / `ma_samples` are computed identically to SNV dosage: `maf = min(ref_alleles, alt_alleles) / (2 × n_samples)`, where `ref_alleles = 2×n0 + n1` and `alt_alleles = n1 + 2×n2`. `--maf-threshold` and `--ma-sample-threshold` work the same as for VCF variants.

### Notes

- `start` is 0-based (BED convention); the site position used for windowing is `start + 1`
- Missing values can be encoded as `NA` or `.` — they are imputed by per-site mean
- Sample IDs in the header must match those in the phenotype BED header
- The file may be plain text or bgzipped (`.bed.gz`)

`--bedmethyl` and `--vcf` are mutually exclusive; exactly one must be provided.

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
| Chunk/parallelisation | Via Python wrapper | Native: pass multiple `-r` regions; processed in parallel with rayon |

---

## Dependencies

- [`clap`](https://crates.io/crates/clap) 4.x — argument parsing
- [`rayon`](https://crates.io/crates/rayon) 1.x — data-parallel region processing

Math functions (lgamma, regularized incomplete beta, Nelder-Mead) are implemented from scratch.
