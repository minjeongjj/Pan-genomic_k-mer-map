# kmer_vcf_gwas

Genome-wide association test of every k-mer in a k-mer presence/absence VCF, for one or more traits.
Referred to https://github.com/voichek/kmersGWAS

## Requirements

- Python 3.8 or later with NumPy and SciPy
- [GEMMA](https://github.com/genetics-statistics/GEMMA) 0.98 (tested with 0.98.1)
- pigz (optional, faster compression of the output)
- bash, awk and gzip for the quick test


## Input

**VCF** (tab-separated, compressed with `bgzip`):

```text
#CHROM  POS  ID                               REF  ALT  QUAL  FILTER  INFO  FORMAT  S1   S2   S3
0       0    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  1    .    .     .       PR    GT      0/0  0/0  0/0
0       0    ACTAATTCATTAGTGGAATAAATTATAGTAA  1    .    .     .       PR    GT      0/0  1/1  1/1
```

- One record per k-mer, with the k-mer sequence in ID.
- GT 0/0 = k-mer present, 1/1 = absent. Missing genotypes are not allowed.
- All records have the same length, as in a VCF written by PLINK 1.9 with CHROM and POS set to 0.
- BGZF compression (`bgzip`, `bcftools -Oz`) is required because the file is read in parallel from block
  boundaries.

If 1/1 means present in your VCF, the P values do not change; `allele1` and `af` then refer to absence.

**Phenotype table** (tab-separated):

```text
sample  Fruit_weight  Plant_height
S1      12.4          85
S2      NA            92
S3      10.1          NA
```

## Usage

```bash
python kmer_vcf_gwas.py --vcf kmers.vcf.gz --pheno traits.tsv --out results --threads 56
```

| Option | Default | Description |
|---|---|---|
| `--vcf` | required | k-mer presence/absence VCF |
| `--pheno` | required | phenotype table |
| `--out` | required | output directory |
| `--threads` | 8 | parallel processes |
| `--maf` | 0.05 | minimum minor-allele frequency of tested k-mers (0: `--mac` only) |
| `--mac` | 5 | minimum minor-allele count of tested k-mers |
| `--gemma` | `gemma` on PATH | GEMMA executable |

### Examples

```bash
python kmer_vcf_gwas.py --vcf kmers.vcf.gz --pheno traits.tsv --out results --threads 56 --gemma /path/to/gemma
```


MAF and MAC can be adjusted in the option.

```bash
python kmer_vcf_gwas.py --vcf kmers.vcf.gz --pheno traits.tsv --out results_maf0.01 --maf 0.01 --mac 3 --threads 56
python kmer_vcf_gwas.py --vcf kmers.vcf.gz --pheno traits.tsv --out results_mac3 --maf 0 --mac 3 --threads 56
```


## Output

```text
results/
└── <trait>/
    └── assoc.txt.gz
```

`assoc.txt.gz` is the GEMMA output for the k-mers with approximate P ≤ 10⁻³, sorted by approximate P. The P value
is `p_lrt` (column 10). `rs` is `<k-mer>_<rank>`; `chr` and `ps` are 0. `allele1` is `P` (present), so `af` is the
fraction of samples that carry the k-mer. The approximate test is slightly conservative, so k-mers with a GEMMA P
close to 10⁻³ can be missing from the file. `l_mle` (column 9) is GEMMA's maximum-likelihood ratio of genetic to
residual variance; 1e5 is GEMMA's upper bound. With a centered kinship matrix the maximum-likelihood fit of a
highly heritable trait can sit at this bound (most `Fruit_weight` rows of the example).

The log has one line per trait. These numbers are not written to a file, so keep the log:

```text
[2026-11-27 16:43:05] Fruit_weight: n 194, 50,783 k-mers tested, Bonferroni -log10 P 6.0067, 39 above, max 14.1745, lambda 0.9110, null pve 0.838158
```

### Reading the results

k-mers at or above the Bonferroni threshold of one trait, as a table and as FASTA

## Restarting

Intermediate files are written to `<out>/work` and removed after every trait has finished. If a run stops, the
same command continues from the last finished step. A trait that already has `<out>/<trait>/assoc.txt.gz` is
skipped, so a trait column added later is run on its own; delete the file to run the trait again. The program
stops if an unfinished `<out>/work` was made from other phenotype values, another VCF (recognized by its file size
and sample names) or another `--maf` or `--mac`.

## Run time

300 pepper genome assemblies (31-mers, 627 GB VCF), 56 threads, measured with the pipeline script this program
was made from (identical results):

- one trait, 280 samples, `--maf 0.01` (7.6 billion k-mers tested): 1 h 21 min, 38 GB peak memory, 13 GB in
  `<out>/work`;
- 120 traits on 47 sample sets, `--maf 0.05`: kinship 26 min, scan 7.8 h, GEMMA 7.1 h. The GEMMA step holds the
  scan results of all traits in memory (about 75 GB here).
