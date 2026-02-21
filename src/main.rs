use std::cmp::Ordering;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Cursor, Write};
use std::process::Command;

use clap::{ArgAction, Parser};

#[derive(Debug, Parser)]
#[command(name = "rust-fastqtl", about = "Fast cis-QTL mapper (Rust port of FastQTL)")]
struct Args {
    /// Input VCF/BCF file (may be gzip-compressed)
    #[arg(short, long)]
    vcf: String,

    /// Input BED phenotype file (may be gzip-compressed)
    #[arg(short, long)]
    bed: String,

    /// Output file path
    #[arg(short, long)]
    out: String,

    /// Covariate file (optional)
    #[arg(short, long)]
    cov: Option<String>,

    /// Genomic region to analyse (chr:start-end or chr)
    #[arg(short, long)]
    region: String,

    /// Cis-window size in bp
    #[arg(short, long, default_value_t = 1_000_000)]
    window: i32,

    /// P-value threshold for nominal-pass output
    #[arg(long, default_value_t = 1.0)]
    threshold: f64,

    /// Minor allele frequency filter
    #[arg(long, default_value_t = 0.0)]
    maf_threshold: f64,

    /// Minimum number of samples carrying the minor allele
    #[arg(long, default_value_t = 0)]
    ma_sample_threshold: usize,

    /// Permutation counts: 1 value = fixed N; 2 values = adaptive (min_hits, max_perms);
    /// 3 values = adaptive with min-perms guard (min_perms, min_hits, max_perms)
    #[arg(short, long, num_args = 1..=3)]
    permute: Option<Vec<usize>>,

    /// Random seed
    #[arg(long, default_value_t = 12345)]
    seed: u64,

    /// Apply rank-normal transformation to phenotypes
    #[arg(short, long, action = ArgAction::SetTrue)]
    normal: bool,
}

#[derive(Clone, Debug)]
struct Region {
    chr: String,
    start: i32,
    end: i32,
}

#[derive(Clone)]
struct Phenotype {
    id: String,
    start: i32,
    values: Vec<f64>,
    sd: f64,
}

#[derive(Clone)]
struct Genotype {
    id: String,
    pos: i32,
    values: Vec<f64>,
    sd: f64,
    maf: f64,
    ma_count: i32,
    ma_samples: usize,
}

struct Residualizer {
    x: Vec<Vec<f64>>, // n_samples x p
    xtx_inv: Option<Vec<Vec<f64>>>,
}

struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    fn new(seed: u64) -> Self {
        let state = if seed == 0 { 0x9e3779b97f4a7c15 } else { seed };
        Self { state }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn gen_usize(&mut self, upper: usize) -> usize {
        if upper <= 1 {
            0
        } else {
            (self.next_u64() % (upper as u64)) as usize
        }
    }

    fn shuffle<T>(&mut self, data: &mut [T]) {
        if data.len() <= 1 {
            return;
        }
        let mut i = data.len() - 1;
        loop {
            let j = self.gen_usize(i + 1);
            data.swap(i, j);
            if i == 0 {
                break;
            }
            i -= 1;
        }
    }
}


fn parse_region(s: &str) -> Result<Region, Box<dyn Error>> {
    if let Some((chr, rest)) = s.split_once(':') {
        let (start_s, end_s) = rest
            .split_once('-')
            .ok_or("region must be chr:start-end")?;
        let start: i32 = start_s.parse()?;
        let end: i32 = end_s.parse()?;
        Ok(Region {
            chr: chr.to_string(),
            start,
            end,
        })
    } else {
        Ok(Region {
            chr: s.to_string(),
            start: 0,
            end: 1_000_000_000,
        })
    }
}

fn open_bufread(path: &str) -> Result<Box<dyn BufRead>, Box<dyn Error>> {
    if path.ends_with(".gz") {
        let out = Command::new("gzip").arg("-cd").arg(path).output()?;
        if !out.status.success() {
            return Err(format!("failed to read gzip file: {}", path).into());
        }
        Ok(Box::new(BufReader::new(Cursor::new(out.stdout))))
    } else {
        let file = File::open(path)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

fn parse_bed(
    path: &str,
    region: &Region,
) -> Result<(Vec<String>, Vec<Phenotype>, HashMap<String, usize>), Box<dyn Error>> {
    let mut reader = open_bufread(path)?;
    let mut line = String::new();

    if reader.read_line(&mut line)? == 0 {
        return Err("empty BED".into());
    }
    let header: Vec<&str> = line.trim_end().split('\t').collect();
    if header.len() < 5 {
        return Err("BED header must have at least 5 columns".into());
    }
    let samples = header[4..].iter().map(|s| (*s).to_string()).collect::<Vec<_>>();

    let mut sample_index = HashMap::new();
    for (i, s) in samples.iter().enumerate() {
        sample_index.insert(s.clone(), i);
    }

    let mut phenotypes = Vec::new();
    line.clear();
    while reader.read_line(&mut line)? > 0 {
        if line.trim().is_empty() {
            line.clear();
            continue;
        }
        let cols: Vec<&str> = line.trim_end().split('\t').collect();
        if cols.len() < 4 + samples.len() {
            line.clear();
            continue;
        }
        let chr = cols[0];
        if chr != region.chr {
            line.clear();
            continue;
        }
        let start0: i32 = cols[1].parse()?;
        let start1 = start0 + 1;
        if start1 < region.start || start1 > region.end {
            line.clear();
            continue;
        }

        let id = cols[3].to_string();
        let mut values = Vec::with_capacity(samples.len());
        for v in &cols[4..] {
            if *v == "NA" || *v == "." {
                values.push(f64::NAN);
            } else {
                values.push(v.parse::<f64>()?);
            }
        }

        phenotypes.push(Phenotype {
            id,
            start: start1,
            values,
            sd: 0.0,
        });
        line.clear();
    }

    Ok((samples, phenotypes, sample_index))
}

fn parse_covariates(
    path: &str,
    sample_index: &HashMap<String, usize>,
    n_samples: usize,
) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
    let mut reader = open_bufread(path)?;
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        return Err("empty covariate file".into());
    }

    let header: Vec<&str> = line.trim_end().split('\t').collect();
    if header.len() < 2 {
        return Err("covariate header must have at least 2 columns".into());
    }

    let mut mapping: Vec<Option<usize>> = Vec::with_capacity(header.len() - 1);
    for s in &header[1..] {
        mapping.push(sample_index.get(*s).copied());
    }

    let mut covariates = Vec::new();
    line.clear();
    while reader.read_line(&mut line)? > 0 {
        if line.trim().is_empty() {
            line.clear();
            continue;
        }
        let cols: Vec<&str> = line.trim_end().split('\t').collect();
        if cols.len() != mapping.len() + 1 {
            return Err(format!("invalid covariate row: {}", line.trim_end()).into());
        }

        let mut row = vec![f64::NAN; n_samples];
        for (i, map_idx) in mapping.iter().enumerate() {
            if let Some(sidx) = map_idx {
                let raw = cols[i + 1];
                row[*sidx] = if raw == "NA" || raw == "." {
                    f64::NAN
                } else {
                    raw.parse::<f64>()?
                };
            }
        }
        impute_nan(&mut row);
        covariates.push(row);
        line.clear();
    }

    Ok(covariates)
}

fn parse_gt_or_ds(sample_field: &str, _fmt_fields: &[&str], idx: usize, gt_mode: bool) -> Option<f64> {
    let fields: Vec<&str> = sample_field.split(':').collect();
    if idx >= fields.len() {
        return None;
    }
    let raw = fields[idx];
    if raw == "." || raw == "./." || raw == ".|." {
        return None;
    }

    if !gt_mode {
        raw.parse::<f64>().ok()
    } else {
        let mut it = raw.split(['/', '|']);
        let a = it.next()?.parse::<i32>().ok()?;
        let b = it.next()?.parse::<i32>().ok()?;
        Some((a + b) as f64)
    }
}

fn maf_stats(values: &[Option<f64>]) -> Option<(f64, i32, usize)> {
    let mut c0 = 0_i32;
    let mut c1 = 0_i32;
    let mut c2 = 0_i32;

    for v in values {
        if let Some(x) = v {
            match x.round() as i32 {
                0 => c0 += 1,
                1 => c1 += 1,
                2 => c2 += 1,
                _ => return None,
            }
        }
    }

    let n = c0 + c1 + c2;
    if n == 0 {
        return None;
    }

    let ref_alleles = 2 * c0 + c1;
    let alt_alleles = c1 + 2 * c2;

    if ref_alleles >= alt_alleles {
        Some((
            alt_alleles as f64 / (2 * n) as f64,
            alt_alleles,
            (c1 + c2) as usize,
        ))
    } else {
        Some((
            ref_alleles as f64 / (2 * n) as f64,
            ref_alleles,
            (c0 + c1) as usize,
        ))
    }
}

fn parse_vcf(
    path: &str,
    region: &Region,
    window: i32,
    sample_index: &HashMap<String, usize>,
    n_samples: usize,
    maf_threshold: f64,
    ma_sample_threshold: usize,
) -> Result<Vec<Genotype>, Box<dyn Error>> {
    let mut reader = open_bufread(path)?;
    let mut line = String::new();

    let mut mapping: Vec<Option<usize>> = Vec::new();
    let mut genotypes = Vec::new();

    let gstart = (region.start - window).max(0);
    let gend = region.end + window;

    while reader.read_line(&mut line)? > 0 {
        if line.starts_with("##") {
            line.clear();
            continue;
        }
        if line.starts_with("#CHROM") {
            let hdr: Vec<&str> = line.trim_end().split('\t').collect();
            if hdr.len() < 10 {
                return Err("VCF header missing sample columns".into());
            }
            mapping.clear();
            for s in &hdr[9..] {
                mapping.push(sample_index.get(*s).copied());
            }
            line.clear();
            continue;
        }
        if line.starts_with('#') {
            line.clear();
            continue;
        }

        let cols: Vec<&str> = line.trim_end().split('\t').collect();
        if cols.len() < 10 {
            line.clear();
            continue;
        }

        let chr = cols[0];
        let pos: i32 = cols[1].parse()?;
        if chr != region.chr || pos < gstart || pos > gend {
            line.clear();
            continue;
        }

        let fmt_fields: Vec<&str> = cols[8].split(':').collect();
        let ds_idx = fmt_fields.iter().position(|f| *f == "DS");
        let gt_idx = fmt_fields.iter().position(|f| *f == "GT");
        let (idx, gt_mode) = if let Some(i) = ds_idx {
            (i, false)
        } else if let Some(i) = gt_idx {
            (i, true)
        } else {
            line.clear();
            continue;
        };

        let mut vals_opt = vec![None; n_samples];
        for (i, sample_field) in cols[9..].iter().enumerate() {
            if let Some(sidx) = mapping.get(i).copied().flatten() {
                vals_opt[sidx] = parse_gt_or_ds(sample_field, &fmt_fields, idx, gt_mode);
            }
        }

        if let Some((maf, ma_count, ma_samples)) = maf_stats(&vals_opt) {
            if maf >= maf_threshold && ma_samples >= ma_sample_threshold {
                let mut values = vals_opt
                    .iter()
                    .map(|v| v.unwrap_or(f64::NAN))
                    .collect::<Vec<f64>>();
                impute_nan(&mut values);
                let id = if cols[2] == "." {
                    format!("{}_{}", chr, pos)
                } else {
                    cols[2].to_string()
                };
                genotypes.push(Genotype {
                    id,
                    pos,
                    values,
                    sd: 0.0,
                    maf,
                    ma_count,
                    ma_samples,
                });
            }
        }

        line.clear();
    }

    Ok(genotypes)
}

fn impute_nan(v: &mut [f64]) {
    let mut sum = 0.0;
    let mut n = 0usize;
    for x in v.iter().copied() {
        if !x.is_nan() {
            sum += x;
            n += 1;
        }
    }
    let mean = if n > 0 { sum / n as f64 } else { 0.0 };
    for x in v.iter_mut() {
        if x.is_nan() {
            *x = mean;
        }
    }
}

fn stddev(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    let ss = v.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>();
    (ss / (v.len() as f64 - 1.0)).sqrt()
}

fn normalize(v: &mut [f64]) {
    if v.is_empty() {
        return;
    }
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    for x in v.iter_mut() {
        *x -= mean;
    }
    let l2 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    let d = if l2 == 0.0 { 1.0 } else { l2 };
    for x in v.iter_mut() {
        *x /= d;
    }
}

fn corr(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y.iter()).map(|(a, b)| a * b).sum()
}

fn slope(c: f64, psd: f64, gsd: f64) -> f64 {
    if psd < 1e-16 || gsd < 1e-16 {
        0.0
    } else {
        c * psd / gsd
    }
}

fn tstat2(c: f64, df: f64) -> f64 {
    let den = 1.0 - c * c;
    if den <= 0.0 {
        f64::INFINITY
    } else {
        df * c * c / den
    }
}

fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * ax);
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let y = 1.0 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * (-ax * ax).exp();
    sign * y
}

fn normal_cdf(x: f64) -> f64 {
    0.5 * (1.0 + erf(x / 2.0_f64.sqrt()))
}

fn inverse_normal_cdf(p: f64) -> f64 {
    // Acklam approximation
    let p = p.clamp(1e-12, 1.0 - 1e-12);
    let a: [f64; 6] = [
        -39.6968302866538,
        220.946098424521,
        -275.928510446969,
        138.357751867269,
        -30.6647980661472,
        2.50662827745924,
    ];
    let b: [f64; 5] = [
        -54.4760987982241,
        161.585836858041,
        -155.698979859887,
        66.8013118877197,
        -13.2806815528857,
    ];
    let c: [f64; 6] = [
        -0.00778489400243029,
        -0.322396458041136,
        -2.40075827716184,
        -2.54973253934373,
        4.37466414146497,
        2.93816398269878,
    ];
    let d: [f64; 4] = [
        0.00778469570904146,
        0.32246712907004,
        2.445134137143,
        3.75440866190742,
    ];

    let plow = 0.02425;
    let phigh = 1.0 - plow;

    if p < plow {
        let q = (-2.0 * p.ln()).sqrt();
        -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    } else if p > phigh {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    } else {
        let q = p - 0.5;
        let r = q * q;
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    }
}

fn pvalue_from_tstat2(t2: f64, _df: f64) -> f64 {
    // FastQTL uses F(1,df); here we use a normal tail approximation for portability.
    if !t2.is_finite() {
        return 0.0;
    }
    let t = t2.max(0.0).sqrt();
    let p = 2.0 * (1.0 - normal_cdf(t));
    p.clamp(0.0, 1.0)
}

// --- Regularized incomplete beta and helpers (for exact F-dist p-values and pbeta) ---

fn lgamma(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    const G: f64 = 7.0;
    if x < 0.5 {
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().ln() - lgamma(1.0 - x)
    } else {
        let xm1 = x - 1.0;
        let mut a = C[0];
        for i in 1..9usize {
            a += C[i] / (xm1 + i as f64);
        }
        let t = xm1 + G + 0.5;
        0.5 * (2.0 * std::f64::consts::PI).ln() + (xm1 + 0.5) * t.ln() - t + a.ln()
    }
}

fn betacf(a: f64, b: f64, x: f64) -> f64 {
    // Continued fraction for regularized incomplete beta (Lentz's method)
    const MAXIT: usize = 200;
    const EPS: f64 = 3.0e-7;
    const FPMIN: f64 = 1.0e-30;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0_f64;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FPMIN { d = FPMIN; }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=MAXIT {
        let mf = m as f64;
        let m2 = 2.0 * mf;
        let aa = mf * (b - mf) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN { d = FPMIN; }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN { c = FPMIN; }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + mf) * (qab + mf) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN { d = FPMIN; }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN { c = FPMIN; }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS { break; }
    }
    h
}

fn pbeta(x: f64, a: f64, b: f64) -> f64 {
    // Regularized incomplete beta I_x(a, b), lower tail
    if x <= 0.0 { return 0.0; }
    if x >= 1.0 { return 1.0; }
    let lbeta_ab = lgamma(a) + lgamma(b) - lgamma(a + b);
    if x < (a + 1.0) / (a + b + 2.0) {
        let front = (a * x.ln() + b * (1.0 - x).ln() - lbeta_ab).exp() / a;
        (front * betacf(a, b, x)).clamp(0.0, 1.0)
    } else {
        let front = (b * (1.0 - x).ln() + a * x.ln() - lbeta_ab).exp() / b;
        (1.0 - front * betacf(b, a, 1.0 - x)).clamp(0.0, 1.0)
    }
}

fn pvalue_from_corr_f(c: f64, df: f64) -> f64 {
    // Exact F(1, df) p-value from Pearson correlation c:
    //   P(F(1,df) > df*c²/(1-c²)) = I_{1-c²}(df/2, 1/2)
    let c2 = c * c;
    if c2 >= 1.0 { return 0.0; }
    pbeta(1.0 - c2, df / 2.0, 0.5)
}

// --- Nelder-Mead minimiser (2-D, derivative-free) ---

fn nelder_mead_2d<F: Fn(f64, f64) -> f64>(
    f: F,
    x0: f64, y0: f64,
    sx: f64, sy: f64,
    tol: f64, max_iter: usize,
) -> (f64, f64) {
    let mut pts = [(x0, y0), (x0 + sx, y0), (x0, y0 + sy)];
    let mut fv = [f(pts[0].0, pts[0].1), f(pts[1].0, pts[1].1), f(pts[2].0, pts[2].1)];
    for _ in 0..max_iter {
        // Sort: pts[0] = best (lowest f), pts[2] = worst
        if fv[0] > fv[1] { pts.swap(0, 1); fv.swap(0, 1); }
        if fv[0] > fv[2] { pts.swap(0, 2); fv.swap(0, 2); }
        if fv[1] > fv[2] { pts.swap(1, 2); fv.swap(1, 2); }
        let size = ((pts[2].0 - pts[0].0).powi(2) + (pts[2].1 - pts[0].1).powi(2)).sqrt();
        if size < tol { break; }
        let cx = (pts[0].0 + pts[1].0) / 2.0;
        let cy = (pts[0].1 + pts[1].1) / 2.0;
        // Reflect worst through centroid of best two
        let rx = 2.0 * cx - pts[2].0;
        let ry = 2.0 * cy - pts[2].1;
        let fr = f(rx, ry);
        if fr < fv[0] {
            // Expansion
            let ex = 3.0 * cx - 2.0 * pts[2].0;
            let ey = 3.0 * cy - 2.0 * pts[2].1;
            let fe = f(ex, ey);
            if fe < fr { pts[2] = (ex, ey); fv[2] = fe; }
            else        { pts[2] = (rx, ry); fv[2] = fr; }
        } else if fr < fv[1] {
            pts[2] = (rx, ry); fv[2] = fr;
        } else {
            // Inside contraction
            let kx = (cx + pts[2].0) / 2.0;
            let ky = (cy + pts[2].1) / 2.0;
            let fk = f(kx, ky);
            if fk <= fv[2] {
                pts[2] = (kx, ky); fv[2] = fk;
            } else {
                // Shrink toward best
                pts[1].0 = (pts[1].0 + pts[0].0) / 2.0;
                pts[1].1 = (pts[1].1 + pts[0].1) / 2.0;
                fv[1] = f(pts[1].0, pts[1].1);
                pts[2].0 = (pts[2].0 + pts[0].0) / 2.0;
                pts[2].1 = (pts[2].1 + pts[0].1) / 2.0;
                fv[2] = f(pts[2].0, pts[2].1);
            }
        }
    }
    if fv[0] <= fv[1] && fv[0] <= fv[2] { pts[0] }
    else if fv[1] <= fv[2] { pts[1] }
    else { pts[2] }
}

// --- Beta distribution MLE via Nelder-Mead ---

fn mle_beta(pvals: &[f64], shape1_init: f64, shape2_init: f64) -> (f64, f64) {
    const S1_MIN: f64 = 0.1;
    const S1_MAX: f64 = 10.0;
    const S2_MIN: f64 = 1.0;
    const S2_MAX: f64 = 1_000_000.0;
    let n = pvals.len() as f64;
    let mut slp = 0.0_f64;
    let mut sl1mp = 0.0_f64;
    for &p in pvals {
        let pp = p.clamp(1e-15, 1.0 - 1e-15);
        slp   += pp.ln();
        sl1mp += (1.0 - pp).ln();
    }
    let neg_ll = |a: f64, b: f64| -> f64 {
        if a < S1_MIN || a > S1_MAX || b < S2_MIN || b > S2_MAX {
            return 1e30;
        }
        let lb = lgamma(a) + lgamma(b) - lgamma(a + b);
        -((a - 1.0) * slp + (b - 1.0) * sl1mp - n * lb)
    };
    let s1 = shape1_init.clamp(S1_MIN, S1_MAX);
    let s2 = shape2_init.clamp(S2_MIN, S2_MAX);
    let (a, b) = nelder_mead_2d(neg_ll, s1, s2, (s1 / 10.0).max(1e-4), (s2 / 10.0).max(1e-4), 0.01, 1000);
    (a.clamp(S1_MIN, S1_MAX), b.clamp(S2_MIN, S2_MAX))
}

// Moment-matching + optional MLE for Beta(shape1, shape2) from p-values
fn fit_beta(pvals: &[f64], n_variants: usize) -> (f64, f64) {
    if pvals.len() < 2 { return (1.0, 1.0); }
    let n = pvals.len() as f64;
    let mean = pvals.iter().sum::<f64>() / n;
    if mean <= 0.0 || mean >= 1.0 { return (1.0, 1.0); }
    let var = pvals.iter().map(|p| (p - mean).powi(2)).sum::<f64>() / (n - 1.0);
    if var <= 0.0 { return (1.0, 1.0); }
    let shape1 = mean * (mean * (1.0 - mean) / var - 1.0);
    let shape2 = shape1 * (1.0 / mean - 1.0);
    if shape1 <= 0.0 || shape2 <= 0.0 { return (1.0, 1.0); }
    if n_variants > 10 {
        mle_beta(pvals, shape1, shape2)
    } else {
        (shape1, shape2)
    }
}

// Golden-section search for effective degrees of freedom (matches C++ learnDF)
fn learn_df(corrs: &[f64], df_init: f64) -> f64 {
    if corrs.len() < 2 { return df_init; }
    let n = corrs.len() as f64;
    let objective = |df: f64| -> f64 {
        let mean = corrs.iter().map(|&c| pvalue_from_corr_f(c, df)).sum::<f64>() / n;
        let var = corrs.iter()
            .map(|&c| { let p = pvalue_from_corr_f(c, df); (p - mean).powi(2) })
            .sum::<f64>() / (n - 1.0);
        if var <= 0.0 { return 1e10; }
        let shape1 = mean * (mean * (1.0 - mean) / var - 1.0);
        (shape1 - 1.0).abs()
    };
    // Golden-section search over [1, max(3*df_init, 100)]
    let phi = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut lo = 1.0_f64;
    let mut hi = (3.0 * df_init).max(100.0);
    let mut x1 = hi - phi * (hi - lo);
    let mut x2 = lo + phi * (hi - lo);
    let mut f1 = objective(x1);
    let mut f2 = objective(x2);
    for _ in 0..50 {
        if hi - lo < 0.01 { break; }
        if f1 < f2 {
            hi = x2; x2 = x1; f2 = f1;
            x1 = hi - phi * (hi - lo);
            f1 = objective(x1);
        } else {
            lo = x1; x1 = x2; f1 = f2;
            x2 = lo + phi * (hi - lo);
            f2 = objective(x2);
        }
    }
    (lo + hi) / 2.0
}

fn rank_normalize(v: &mut [f64]) {
    let mut indexed = v
        .iter()
        .copied()
        .enumerate()
        .collect::<Vec<(usize, f64)>>();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

    let n = v.len();
    let mut ranks = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && indexed[j].1 == indexed[i].1 {
            j += 1;
        }
        let r = (i + j + 1) as f64 / 2.0;
        for k in i..j {
            ranks[indexed[k].0] = r;
        }
        i = j;
    }

    for (i, r) in ranks.iter().enumerate() {
        let p = (r - 0.5) / n as f64;
        v[i] = inverse_normal_cdf(p);
    }
}

fn invert_matrix(mut a: Vec<Vec<f64>>) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    if n == 0 || a[0].len() != n {
        return None;
    }

    let mut inv = vec![vec![0.0; n]; n];
    for (i, row) in inv.iter_mut().enumerate().take(n) {
        row[i] = 1.0;
    }

    for i in 0..n {
        let mut pivot = i;
        let mut best = a[i][i].abs();
        for (r, row) in a.iter().enumerate().take(n).skip(i + 1) {
            if row[i].abs() > best {
                best = row[i].abs();
                pivot = r;
            }
        }
        if best < 1e-12 {
            return None;
        }
        if pivot != i {
            a.swap(i, pivot);
            inv.swap(i, pivot);
        }

        let d = a[i][i];
        for c in 0..n {
            a[i][c] /= d;
            inv[i][c] /= d;
        }

        for r in 0..n {
            if r == i {
                continue;
            }
            let f = a[r][i];
            if f == 0.0 {
                continue;
            }
            for c in 0..n {
                a[r][c] -= f * a[i][c];
                inv[r][c] -= f * inv[i][c];
            }
        }
    }

    Some(inv)
}

impl Residualizer {
    fn new(covariates: &[Vec<f64>], n_samples: usize) -> Self {
        if covariates.is_empty() {
            return Self {
                x: Vec::new(),
                xtx_inv: None,
            };
        }

        let p = covariates.len() + 1;
        let mut x = vec![vec![0.0; p]; n_samples];
        for i in 0..n_samples {
            x[i][0] = 1.0;
            for (j, cov) in covariates.iter().enumerate() {
                x[i][j + 1] = cov[i];
            }
        }

        let mut xtx = vec![vec![0.0; p]; p];
        for i in 0..p {
            for j in 0..p {
                let mut s = 0.0;
                for row in x.iter().take(n_samples) {
                    s += row[i] * row[j];
                }
                xtx[i][j] = s;
            }
        }

        let xtx_inv = invert_matrix(xtx);
        Self { x, xtx_inv }
    }

    fn residualize(&self, y: &mut [f64]) {
        let Some(inv) = &self.xtx_inv else {
            return;
        };

        let n = y.len();
        let p = inv.len();
        let mut xty = vec![0.0; p];
        for (j, xj) in xty.iter_mut().enumerate().take(p) {
            let mut s = 0.0;
            for (i, yi) in y.iter().enumerate().take(n) {
                s += self.x[i][j] * yi;
            }
            *xj = s;
        }

        let mut beta = vec![0.0; p];
        for i in 0..p {
            for (j, xty_j) in xty.iter().enumerate().take(p) {
                beta[i] += inv[i][j] * xty_j;
            }
        }

        for (i, yi) in y.iter_mut().enumerate().take(n) {
            let mut fit = 0.0;
            for (j, bj) in beta.iter().enumerate().take(p) {
                fit += self.x[i][j] * bj;
            }
            *yi -= fit;
        }
    }
}

fn run_nominal(
    out_path: &str,
    phenotypes: &[Phenotype],
    genotypes: &[Genotype],
    cis_window: i32,
    threshold: f64,
    n_cov: usize,
) -> Result<(), Box<dyn Error>> {
    let mut out = BufWriter::new(File::create(out_path)?);

    for p in phenotypes {
        for g in genotypes {
            let dist = g.pos - p.start;
            if dist.abs() > cis_window {
                continue;
            }
            let c = corr(&g.values, &p.values);
            let df = (p.values.len() as i32 - 2 - n_cov as i32) as f64;
            if df <= 0.0 {
                continue;
            }
            let t2 = tstat2(c, df);
            let pval = pvalue_from_tstat2(t2, df);
            if pval <= threshold {
                let b = slope(c, p.sd, g.sd);
                let bse = if t2.is_finite() && t2 > 0.0 {
                    b.abs() / t2.sqrt()
                } else {
                    f64::NAN
                };
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    p.id, g.id, dist, g.ma_samples, g.ma_count, g.maf, pval, b, bse
                )?;
            }
        }
    }

    Ok(())
}

fn run_permutation(
    out_path: &str,
    phenotypes: &[Phenotype],
    genotypes: &[Genotype],
    cis_window: i32,
    n_cov: usize,
    permute: &[usize],
    seed: u64,
) -> Result<(), Box<dyn Error>> {
    let mut out = BufWriter::new(File::create(out_path)?);
    let mut rng = XorShift64::new(seed);

    for (pi, p) in phenotypes.iter().enumerate() {
        let target = genotypes
            .iter()
            .enumerate()
            .filter_map(|(gi, g)| {
                let d = g.pos - p.start;
                if d.abs() <= cis_window {
                    Some((gi, d))
                } else {
                    None
                }
            })
            .collect::<Vec<(usize, i32)>>();

        if target.is_empty() {
            writeln!(out, "{}\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA", p.id)?;
            continue;
        }

        let mut best_corr: f64 = 0.0;
        let mut best_g = target[0].0;
        let mut best_d = target[0].1;
        for (gi, d) in &target {
            let c = corr(&genotypes[*gi].values, &p.values);
            if c.abs() > best_corr.abs() || (c.abs() == best_corr.abs() && d.abs() < best_d.abs()) {
                best_corr = c;
                best_g = *gi;
                best_d = *d;
            }
        }

        let mut n_perm = 0usize;
        let mut n_better = 0usize;
        let mut yperm = p.values.clone();
        let mut perm_corrs: Vec<f64> = Vec::new();

        loop {
            yperm.copy_from_slice(&p.values);
            rng.shuffle(&mut yperm);
            normalize(&mut yperm);

            let mut best_pcorr: f64 = 0.0;
            for (gi, _) in &target {
                let c = corr(&genotypes[*gi].values, &yperm);
                if c.abs() > best_pcorr.abs() {
                    best_pcorr = c;
                }
            }

            if best_pcorr.abs() >= best_corr.abs() {
                n_better += 1;
            }
            perm_corrs.push(best_pcorr);
            n_perm += 1;

            let done = match permute.len() {
                1 => n_perm >= permute[0],
                2 => n_better >= permute[0] || n_perm >= permute[1],
                3 => n_perm >= permute[0] && (n_better >= permute[1] || n_perm >= permute[2]),
                _ => return Err("--permute expects 1..3 integers".into()),
            };
            if done {
                break;
            }
        }

        // Degrees of freedom: learn effective df from permutation correlations
        let true_df = (p.values.len() as i32 - 2 - n_cov as i32) as f64;
        let perm_var = {
            let m = perm_corrs.iter().sum::<f64>() / perm_corrs.len() as f64;
            perm_corrs.iter().map(|c| (c - m).powi(2)).sum::<f64>()
                / (perm_corrs.len() as f64 - 1.0)
        };
        let eff_df = if perm_var > 0.0 {
            learn_df(&perm_corrs, true_df)
        } else {
            true_df
        };

        // Compute permutation p-values with effective df, then fit Beta(shape1, shape2)
        let perm_pvals: Vec<f64> = perm_corrs.iter()
            .map(|&c| pvalue_from_corr_f(c, eff_df.max(1.0)))
            .collect();
        let (beta_shape1, beta_shape2) = fit_beta(&perm_pvals, target.len());

        let t2 = tstat2(best_corr, true_df.max(1.0));
        let p_nom = pvalue_from_tstat2(t2, true_df.max(1.0));
        // Nominal p-value with effective df (used for beta-adjusted p-value)
        let p_nom_eff = pvalue_from_corr_f(best_corr, eff_df.max(1.0));
        let b = slope(best_corr, p.sd, genotypes[best_g].sd);
        let bse = if t2.is_finite() && t2 > 0.0 {
            b.abs() / t2.sqrt()
        } else {
            f64::NAN
        };
        let p_emp = (n_better as f64 + 1.0) / (n_perm as f64 + 1.0);
        // Beta-adjusted p-value: pbeta(p_nom_eff, shape1, shape2)
        let p_beta = pbeta(p_nom_eff, beta_shape1, beta_shape2);

        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            p.id,
            target.len(),
            beta_shape1,
            beta_shape2,
            eff_df,
            p_nom_eff,
            genotypes[best_g].id,
            best_d,
            genotypes[best_g].ma_samples,
            genotypes[best_g].ma_count,
            genotypes[best_g].maf,
            p_nom,
            b,
            bse,
            p_emp,
            p_beta
        )?;

        if (pi + 1) % 100 == 0 {
            eprintln!("processed {} phenotypes", pi + 1);
        }
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    if args.window <= 0 {
        return Err("--window must be positive".into());
    }
    if !(0.0 < args.threshold && args.threshold <= 1.0) {
        return Err("--threshold must satisfy 0 < x <= 1".into());
    }
    if !(0.0..0.5).contains(&args.maf_threshold) {
        return Err("--maf-threshold must satisfy 0 <= x < 0.5".into());
    }

    let region = parse_region(&args.region)?;

    let (samples, mut phenotypes, sample_index) = parse_bed(&args.bed, &region)?;
    if phenotypes.is_empty() {
        return Err("no phenotypes found in selected region".into());
    }

    for p in &mut phenotypes {
        impute_nan(&mut p.values);
        if args.normal {
            rank_normalize(&mut p.values);
        }
    }

    let covariates = if let Some(cov_path) = args.cov.as_ref() {
        parse_covariates(cov_path, &sample_index, samples.len())?
    } else {
        Vec::new()
    };

    let mut genotypes = parse_vcf(
        &args.vcf,
        &region,
        args.window,
        &sample_index,
        samples.len(),
        args.maf_threshold,
        args.ma_sample_threshold,
    )?;
    if genotypes.is_empty() {
        return Err("no genotypes found in selected cis window after filters".into());
    }

    let residualizer = Residualizer::new(&covariates, samples.len());

    for g in &mut genotypes {
        residualizer.residualize(&mut g.values);
        g.sd = stddev(&g.values);
        normalize(&mut g.values);
    }
    for p in &mut phenotypes {
        residualizer.residualize(&mut p.values);
        p.sd = stddev(&p.values);
        normalize(&mut p.values);
    }

    if let Some(permute) = args.permute.as_ref() {
        run_permutation(
            &args.out,
            &phenotypes,
            &genotypes,
            args.window,
            covariates.len(),
            permute,
            args.seed,
        )?;
    } else {
        run_nominal(
            &args.out,
            &phenotypes,
            &genotypes,
            args.window,
            args.threshold,
            covariates.len(),
        )?;
    }

    eprintln!(
        "done: phenotypes={}, genotypes={}, covariates={}",
        phenotypes.len(),
        genotypes.len(),
        covariates.len()
    );

    Ok(())
}
