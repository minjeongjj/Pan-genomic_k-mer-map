# 3. Targeted k-mer screening (TKS)

This module screens k-mers derived from targeted genomic sequences for trait associations. Targeted sequence counts define presence/absence, whereas whole-genome counts from [Module 1](../1.Pangenomic_kmer_map/) provide copy-number adjustment. The original screening engine is retained unchanged.

## Requirements

Linux, Bash, Python 3.8+ (recommended), GNU sort, pigz, and the dependencies of Module 1 (including Jellyfish). Install the Python statistical dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

Run commands from `3.TKS/`. The sibling Module 1 directory must be present for Step 01; use `--module1-dir` or `MODULE1_DIR` for another location.

## Inputs

1. `targeted_genome_list.tsv`: two tab-separated columns (`Accession`, `Genome`), listing targeted FASTA/FASTA.gz paths for the entire screening panel, not just the target accessions. See `examples/targeted_genome_list.tsv` (placeholder paths). Include zero-target accessions as valid inputs representing no target sequence, if supported by your Module 1 installation; empty FASTA handling has not been verified.
2. Whole-genome count matrix (`.tsv.gz` or `.tsv`) previously constructed using Module 1, with the same k-mer length, canonical counting and accession assemblies. This is required in `methods` mode, not in `none` mode. Targeted regions should be derived from the corresponding assemblies. Counts should not be binary-converted.
3. `phenotype.tsv`: two whitespace-separated columns, accession ID and numeric phenotype; optional header on first line. Missing phenotypes are omitted from the Welch test.
4. `target_accessions.txt`: one accession ID per line, no header. These define the specificity group. Other PAV accessions form the comparison group for specificity, even if phenotype values are missing.

## Workflow and individual commands

### 01. Build the targeted k-mer map

```bash
 time python3 01_build_targeted_kmer_map.py \
    --input targeted_genome_list.tsv --results results/01_targeted_map \
    --prefix targeted --kmer-size 31 \
    --catalog-threads 30 --catalog-hash 3G \
    --individual-threads 10 --individual-hash 1G \
    --count-jobs 1 --profile-jobs 1 --compression-threads 4 --merge-threads 8
```

One line:

```bash
time python3 01_build_targeted_kmer_map.py --input targeted_genome_list.tsv --results results/01_targeted_map --prefix targeted --kmer-size 31 --catalog-threads 30 --catalog-hash 3G --individual-threads 10 --individual-hash 1G --count-jobs 1 --profile-jobs 1 --compression-threads 4 --merge-threads 8
```

Output: `results/01_targeted_map/05_map/targeted.kmer_merged.table.tsv.gz`. This wrapper calls Module 1 without duplicating its code. Choose hash sizes, threads and jobs for your hardware; these example values are not memory limits. If a targeted matrix already exists, skip Step 01 and use its path in Step 02.

### 02. Prepare the PAV table

```bash
time python3 02_prepare_pav.py \
    --map results/01_targeted_map/05_map/targeted.kmer_merged.table.tsv.gz \
    --phenotype phenotype.tsv --target target_accessions.txt \
    --out-dir results/02_pav --kmer-size 31 --sort-memory 2G --sort-threads 4
```

```bash
time python3 02_prepare_pav.py --map results/01_targeted_map/05_map/targeted.kmer_merged.table.tsv.gz --phenotype phenotype.tsv --target target_accessions.txt --out-dir results/02_pav --kmer-size 31 --sort-memory 2G --sort-threads 4
```

Validates accession IDs, k-mer lengths, and nonnegative integer counts; separates the header, sorts by k-mer under `LC_ALL=C`, and collapses identical duplicates. Conflicting duplicates cause an error. Original targeted counts are retained: the statistical engine interprets counts >=1 as presence and 0 as absence.

Outputs: `PAV_table.tsv` (headerless), `PAV_header.tsv`, `PAV_summary.tsv` under `results/02_pav/`.

### 03. Extract whole-genome CNV counts

```bash
time python3 03_extract_cnv.py \
    --pav results/02_pav/PAV_table.tsv --pav-header results/02_pav/PAV_header.tsv \
    --genome-map /path/to/genome_map.tsv.gz --out-dir results/03_cnv \
    --kmer-size 31 --sort-memory 2G --sort-threads 4 --decompression-threads 2
```

```bash
time python3 03_extract_cnv.py --pav results/02_pav/PAV_table.tsv --pav-header results/02_pav/PAV_header.tsv --genome-map /path/to/genome_map.tsv.gz --out-dir results/03_cnv --kmer-size 31 --sort-memory 2G --sort-threads 4 --decompression-threads 2
```

Loads only targeted k-mer keys into a set, streams the genome-wide matrix once, retains exact matches, and selects/reorders accession columns to match the PAV header. Only retained rows are sorted. Missing PAV accession IDs cause an error. Missing k-mers are reported and cause an error by default; `--allow-missing` permits the original engine to skip unmatched k-mers rather than inventing zero CNV counts. No new Jellyfish counting is performed.

Outputs: `CNV_table.tsv` (headerless), `CNV_header.tsv`, `CNV_summary.tsv`, `missing_kmers.txt` under `results/03_cnv/`.

### 04. Run TKS

```bash
time python3 04_targeted_kmer_screening.py \
    --phenotype phenotype.tsv --resistant target_accessions.txt \
    --pav results/02_pav/PAV_table.tsv --pav-header results/02_pav/PAV_header.tsv \
    --cnv results/03_cnv/CNV_table.tsv --cnv-header results/03_cnv/CNV_header.tsv \
    --copy-adjust methods --weight-scope all --min-resistant 3 --max-other 1 \
    --alpha 0.05 --already-sorted --out results/04_screening/TKS.tsv
```

```bash
time python3 04_targeted_kmer_screening.py --phenotype phenotype.tsv --resistant target_accessions.txt --pav results/02_pav/PAV_table.tsv --pav-header results/02_pav/PAV_header.tsv --cnv results/03_cnv/CNV_table.tsv --cnv-header results/03_cnv/CNV_header.tsv --copy-adjust methods --weight-scope all --min-resistant 3 --max-other 1 --alpha 0.05 --already-sorted --out results/04_screening/TKS.tsv
```

The original engine filters k-mers by target/other presence counts, compares phenotyped carriers and noncarriers with a two-sided Welch t-test, optionally adjusts carrier mean and variance using genomic copy numbers, and applies Benjamini–Hochberg correction to valid tested k-mers. In `methods` mode the weight is the number of phenotyped carriers in the selected scope divided by the sum of their genome-wide counts. `--weight-scope target` restricts weighting to target carriers, not the Welch sample groups. `--copy-adjust none` disables adjustment and does not require CNV files. `--already-sorted` avoids redundant sorting after Steps 02–03. `--output-format csv` selects CSV output. `--resistant` is the original engine's name for the target accession list. CNV counts smaller than targeted counts are logged by the original engine, but not invariably excluded; inspect validation output.

Outputs: `TKS.tsv` (BH-significant), `TKS.all.tsv` (all tested), `TKS.validation.tsv` (CNV consistency findings).

## Run all steps

```bash
time bash Run_pipeline.sh targeted_genome_list.tsv /path/to/genome_map.tsv.gz phenotype.tsv target_accessions.txt results > TKS.log 2>&1
```

Environment variables: `KMER_SIZE` (31), `CATALOG_THREADS` (30), `CATALOG_HASH` (3G), `INDIVIDUAL_THREADS` (10), `INDIVIDUAL_HASH` (1G), `COUNT_JOBS` (1), `PROFILE_JOBS` (1), `COMPRESSION_THREADS` (4), `MERGE_THREADS` (8), `SORT_MEMORY` (2G), `SORT_THREADS` (4), `SORT_TMPDIR` (optional), `DECOMPRESSION_THREADS` (2), `COPY_ADJUST` (methods), `WEIGHT_SCOPE` (all), `MIN_TARGET` (3), `MAX_OTHER` (1), `ALPHA` (0.05), `OUTPUT_FORMAT` (tsv), `TKS_START_STEP` (1), `TKS_STOP_STEP` (4). For `COPY_ADJUST=none`, pass `-` instead of the genome-wide map; Step 03 is skipped.

To reuse an existing targeted map:

```bash
time env TKS_START_STEP=2 TARGETED_MAP=/path/to/targeted_map.tsv.gz bash Run_pipeline.sh - /path/to/genome_map.tsv.gz phenotype.tsv target_accessions.txt results
```

Use a new results directory when regenerating preprocessing outputs; Steps 02–03 refuse to overwrite existing named outputs. The original Step 04 engine may overwrite its output paths. Keep the same input assembly panel and counting conditions throughout.

## Small test (no Jellyfish or genome FASTA required)

```bash
env TKS_START_STEP=2 KMER_SIZE=3 TARGETED_MAP=examples/targeted_map.tsv bash Run_pipeline.sh - examples/genome_map.tsv examples/phenotype.tsv examples/target_accessions.txt test_results
```

Expected: three PAV and CNV rows, two tested k-mers, no BH-significant k-mers at alpha 0.05. The example genome map deliberately has reordered columns and an extra accession. These are synthetic data, not biological results.

## Resource and reproducibility notes

The genome-wide matrix is scanned once but can still take substantial time to decompress and read. Target key storage scales with the number of targeted k-mers. GNU sort requires temporary disk space for retained matrices. Intermediate PAV/CNV tables are uncompressed to match the original statistical engine; the large input maps can remain gzip-compressed. Do not change the original statistical engine without revalidating results. The source engine is preserved byte-for-byte; the synthetic uncompressed Steps 02–04 were tested. Real Jellyfish/pigz execution and large-scale performance are not verified here.
