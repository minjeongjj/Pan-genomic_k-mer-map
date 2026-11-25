# 1. Pan-genomic k-mer map construction

This directory contains the first module of the repository, which constructs an accession-resolved pan-genomic *k*-mer map from genome assemblies. The resulting count matrix provides a genomic variation resource for downstream analyses, including k-mer-based genome-wide association studies (k-mer-GWAS) and TKS.

The module identifies unique canonical *k*-mers across the input assemblies, quantifies their occurrence in each accession, and integrates accession-specific counts into a compressed matrix. K-mers absent from an individual assembly are assigned a count of zero, preserving the order of the pan-genomic catalog across all accessions.

## Scope

This README covers Module 1 only, from genome FASTA input to construction and validation of the pan-genomic *k*-mer count matrix. The downstream k-mer-GWAS and TKS modules are documented separately in their respective directories. The five steps below are internal processing steps of Module 1, not separate repository modules.

Run the commands below from this module's directory.

## Module 1 workflow


| Step | Script | Description |
| --- | --- | --- |
| 01 | `01_build_kmer_catalog.py` | Count canonical k-mers across all assemblies and generate the ordered catalog. |
| 02 | `02_count_individual_kmers.py` | Count canonical k-mers independently for each accession. |
| 03 | `03_build_kmer_profiles.py` | Match individual counts against the catalog and assign zero to absent k-mers. |
| 04 | `04_validate_kmer_profiles.py` | Validate accession-specific profiles. |
| 05 | `05_merge_kmer_profiles.py` | Merge profiles in manifest order, validate the matrix, and compress the output. |

## Requirements

- Linux and Bash
- Python 3.8 or later (recommended baseline; no third-party Python packages required)
- Jellyfish 2
- pigz
- GNU coreutils (`tee`, `paste`, `wc`) and `awk`
- Optional: `cc` or `gcc` for the native streaming matrix validator (awk fallback when no compiler is found)

Install the main dependencies using Conda, for example:

```bash
conda create -n pan_kmer -c conda-forge -c bioconda python=3.10 jellyfish pigz
conda activate pan_kmer
```

Check that `jellyfish --version`, `pigz --version` and `python3 --version` work in your environment. Python 3.6 is not a verified compatibility target for this release.

## Input data

Prepare a tab-separated manifest named `genome_list.tsv` with exactly two columns:

```tsv
Accession	Genome
Acc001	/path/to/Acc001.fasta.gz
Acc002	/path/to/Acc002.fasta
Acc003	/path/to/Acc003.fa.gz
```

- Plain and gzip-compressed FASTA files can be mixed. Gzip is detected from the file signature.
- Genome paths can be absolute or relative to the manifest directory.
- Accession IDs must be unique, start with a letter, digit or underscore, and contain only ASCII letters, digits, underscores, periods or hyphens.
- The manifest order determines the accession column order in the final matrix.
- Use the same genome set, accession IDs, k-mer length, prefix and results directory throughout. Do not modify inputs or intermediate results during a run.

## Quick start

Run all five internal steps from this module directory:

```bash
time bash Run_pipeline.sh genome_list.tsv results pepper > Run_pipeline.log 2>&1
```

To specify resource allocations and compression settings:

```bash
time env \
    KMER_SIZE=31 \
    CATALOG_THREADS=100 \
    CATALOG_HASH=15G \
    CATALOG_COMPRESSION_THREADS=8 \
    INDIVIDUAL_THREADS=30 \
    INDIVIDUAL_HASH=3G \
    COUNT_JOBS=4 \
    PROFILE_JOBS=1 \
    VALIDATION_JOBS=1 \
    COMPRESSION_THREADS=4 \
    COMPRESSION_LEVEL=1 \
    DECOMPRESSION_THREADS=2 \
    CATALOG_SOURCE=header \
    MERGE_THREADS=150 \
    MERGE_DECOMPRESSION_THREADS=1 \
    VALIDATION_MODE=basic \
    bash Run_pipeline.sh genome_list.tsv results pepper > Run_pipeline.log 2>&1
```

Single-line equivalent:

```bash
time env KMER_SIZE=31 CATALOG_THREADS=100 CATALOG_HASH=15G CATALOG_COMPRESSION_THREADS=8 INDIVIDUAL_THREADS=30 INDIVIDUAL_HASH=3G COUNT_JOBS=4 PROFILE_JOBS=1 VALIDATION_JOBS=1 COMPRESSION_THREADS=4 COMPRESSION_LEVEL=1 DECOMPRESSION_THREADS=2 CATALOG_SOURCE=header MERGE_THREADS=150 MERGE_DECOMPRESSION_THREADS=1 VALIDATION_MODE=basic bash Run_pipeline.sh genome_list.tsv results pepper > Run_pipeline.log 2>&1
```

These resource allocations are examples, not automatic hardware limits. Adjust them to available RAM, CPUs and storage throughput before running.

## Step-by-step execution

All scripts provide `--help`. The common options are `--input`, `--results`, `--prefix`, `--compression-level` (default: 1), and `--decompression-threads` (default: 2).

### Step 01. Build the pan-genomic k-mer catalog

Genome sequences are streamed into Jellyfish using canonical counting (`-C`). The dump, ordered k-mer column and unique k-mer count are produced from a single dump stream to avoid repeated decompression and file traversal.

```bash
time python3 01_build_kmer_catalog.py \
    --input genome_list.tsv \
    --results results \
    --prefix pepper \
    --kmer-size 31 \
    --threads 100 \
    --hash-size 15G \
    --compression-threads 8
```

```bash
time python3 01_build_kmer_catalog.py --input genome_list.tsv --results results --prefix pepper --kmer-size 31 --threads 100 --hash-size 15G --compression-threads 8
```

Outputs in `results/01_catalog/`:

- `pepper.cat.jellyfish.dump.gz`: unique canonical k-mers and pooled counts (no header).
- `pepper.cat.jellyfish.header.gz`: ordered k-mer sequences, with `Kmer` as the first line.
- `pepper.cat.jellyfish.dump.line`: number of unique k-mers (excluding the matrix header), followed by the dump path.

Use `--keep-jf` to retain the Jellyfish binary database. Two compressors operate simultaneously in this step; `--compression-threads` applies to each compressor.

### Step 02. Count individual k-mers

Each accession is counted independently using the same k-mer length and canonical counting settings as Step 01. Dump files are compressed directly without retaining uncompressed text dumps.

```bash
time python3 02_count_individual_kmers.py \
    --input genome_list.tsv \
    --results results \
    --prefix pepper \
    --kmer-size 31 \
    --threads 30 \
    --hash-size 3G \
    --jobs 4 \
    --compression-threads 4
```

```bash
time python3 02_count_individual_kmers.py --input genome_list.tsv --results results --prefix pepper --kmer-size 31 --threads 30 --hash-size 3G --jobs 4 --compression-threads 4
```

`--threads` is the number of Jellyfish threads per accession; `--jobs` is the maximum number of concurrent accession jobs. With the example above, up to four count jobs use 120 configured counting threads in total, plus auxiliary threads where applicable. `--hash-size` specifies hash capacity, not a RAM limit in bytes. Use `--keep-jf` if the individual databases are needed later.

Outputs: `results/02_individual_kmers/<Accession>.jellyfish.dump.gz` and `<Accession>.log`.

### Step 03. Build accession-specific k-mer profiles

The original dictionary-based matching algorithm is retained. Individual k-mer counts are loaded into a Python dictionary, then the pan-genomic catalog is streamed in its original order. The count is retrieved for every catalog k-mer; absent k-mers receive zero. Byte-string keys and batched output reduce conversion and write overhead. Compression/decompression subprocesses are started before the dictionary is populated to reduce fork-related memory allocation failures.

```bash
time python3 03_build_kmer_profiles.py \
    --input genome_list.tsv \
    --results results \
    --prefix pepper \
    --jobs 1 \
    --compression-threads 4 \
    --catalog-source header
```

```bash
time python3 03_build_kmer_profiles.py --input genome_list.tsv --results results --prefix pepper --jobs 1 --compression-threads 4 --catalog-source header
```

The default `--catalog-source header` reads the ordered k-mer-only column; use `--catalog-source dump` to read the original two-column catalog instead. Both must originate from the same catalog and preserve identical order.

For example, a catalog containing `AAA`, `AAC`, `AAG` and `ACC`, with individual counts `AAA 3` and `AAG 2`, produces:

```text
Acc001
3
0
2
0
```

Each profile has exactly one accession header followed by one count for every catalog k-mer. Thus a catalog with N unique k-mers requires N + 1 profile lines. Outputs: `results/03_profiles/<Accession>.row.list.gz`, its `.log` and `.meta.json` files.

**Memory warning:** the individual dictionary is held in RAM. Gzip reduces storage requirements but does not reduce dictionary size. Start with `--jobs 1` and measure peak memory before increasing concurrency. The subprocess ordering fix does not eliminate the possibility of genuine out-of-memory errors.

### Step 04. Validate profiles

```bash
time python3 04_validate_kmer_profiles.py \
    --input genome_list.tsv \
    --results results \
    --prefix pepper \
    --jobs 1
```

```bash
time python3 04_validate_kmer_profiles.py --input genome_list.tsv --results results --prefix pepper --jobs 1
```

By default, matching creation metadata is used to check accession identity, row count and file/catalog stamps without rereading the profile. Legacy profiles or profiles without matching metadata are rescanned. To force a complete readback of all profiles, use `--rescan`:

```bash
python3 04_validate_kmer_profiles.py --input genome_list.tsv --results results --prefix pepper --jobs 1 --rescan
```

A rescan checks gzip decoding, the accession header, nonnegative integer counts and the expected line count. Metadata-based validation is not a checksum or proof that every value is correct.

Outputs: `results/04_validation/profile_validation.tsv` and `profile_stamps.json`.

### Step 05. Merge profiles and validate the pan-genomic k-mer map

```bash
time python3 05_merge_kmer_profiles.py \
    --input genome_list.tsv \
    --results results \
    --prefix pepper \
    --threads 150 \
    --validation basic
```

```bash
time python3 05_merge_kmer_profiles.py --input genome_list.tsv --results results --prefix pepper --threads 150 --validation basic
```

The ordered k-mer column and all individual profiles are decompressed as streams, combined using GNU `paste`, validated, and compressed with pigz. An uncompressed merged matrix is not retained. Step 04 must pass before merging.

| Validation | Checks |
| --- | --- |
| `basic` (default) | Exact header and accession order, row and column counts, empty fields, subprocess exit statuses. |
| `full` | All basic checks plus exact k-mer order, nonnegative integer counts and the sum of accession counts against the pooled catalog count, using a further readback. |
| `--verify-gzip` | Additional gzip integrity readback in basic mode. |

For a small test dataset, use `--validation full`. This check can detect order or total-count discrepancies but cannot detect every compensating error between accessions.

Outputs: `results/05_map/pepper.kmer_merged.table.tsv.gz` and `map_validation.tsv`.

## Output format and downstream use

The final compressed TSV contains one row per unique canonical k-mer and one column per accession, in manifest order:

```tsv
Kmer	Acc001	Acc002	Acc003
AAA	3	1	0
AAC	0	2	1
AAG	2	0	4
ACC	0	1	0
```

The example uses short illustrative k-mers; the default analysis uses 31-mers. Counts represent occurrences in the input genome assemblies, not sequencing-read depth. Zero indicates absence from the corresponding assembly's k-mer dump, not necessarily biological absence from the organism.

```text
results/
├── 01_catalog/
│   ├── pepper.cat.jellyfish.dump.gz
│   ├── pepper.cat.jellyfish.header.gz
│   └── pepper.cat.jellyfish.dump.line
├── 02_individual_kmers/
│   ├── Acc001.jellyfish.dump.gz
│   └── Acc001.log
├── 03_profiles/
│   ├── Acc001.row.list.gz
│   ├── Acc001.row.list.gz.log
│   └── Acc001.row.list.gz.meta.json
├── 04_validation/
│   ├── profile_validation.tsv
│   └── profile_stamps.json
└── 05_map/
    ├── pepper.kmer_merged.table.tsv.gz
    └── map_validation.tsv
```

The final `.tsv.gz` matrix is the primary output of Module 1. Retain the accession identifiers and catalog order when preparing inputs for downstream modules. Any module-specific filtering or format conversion is described in the corresponding downstream documentation.

## Resource settings

| Environment variable (runner) | Default | Meaning |
| --- | --- | --- |
| `START_STEP` / `STOP_STEP` | 1 / 5 | Inclusive range of steps to execute |
| `KMER_SIZE` | 31 | k-mer length |
| `CATALOG_THREADS` | 100 | Pooled Jellyfish counting threads |
| `CATALOG_HASH` | 15G | Pooled hash size |
| `CATALOG_COMPRESSION_THREADS` | 8 | Threads per pooled compressor |
| `INDIVIDUAL_THREADS` | 30 | Jellyfish threads per accession |
| `INDIVIDUAL_HASH` | 3G | Hash size per accession |
| `COUNT_JOBS` | 1 | Concurrent accession counting jobs |
| `PROFILE_JOBS` | 1 | Concurrent memory-intensive matching jobs |
| `VALIDATION_JOBS` | 1 | Concurrent profile validation jobs |
| `COMPRESSION_THREADS` | 4 | Compression threads per individual job |
| `COMPRESSION_LEVEL` | 1 | pigz compression level, 1–9 |
| `DECOMPRESSION_THREADS` | 2 | pigz decompression/helper threads |
| `CATALOG_SOURCE` | header | `header` or `dump` for matching |
| `MERGE_THREADS` | 150 | Final matrix compression threads |
| `MERGE_DECOMPRESSION_THREADS` | 1 | Decompression threads per merge input |
| `VALIDATION_MODE` | basic | `basic` or `full` |

Compression level 1 favors speed, but may consume more disk space than level 6 or 9. Standard gzip inflation is predominantly single-threaded per stream, so increasing decompression threads does not provide linear speedup. Step 05 opens one stream per accession plus the catalog; check `ulimit -n` and the process limits before merging hundreds of columns. If permitted, a higher file descriptor limit (e.g. `ulimit -n 4096`) may be useful.

## Restarting and reusing intermediate results

To run Steps 02–05 using a completed catalog:

```bash
time env START_STEP=2 STOP_STEP=5 bash Run_pipeline.sh genome_list.tsv results pepper > Run_pipeline.log 2>&1
```

All three catalog outputs must be present and originate from the same genome set and counting settings. Jellyfish dump ordering may differ between independent runs: never mix profiles generated against different catalogs, even if their k-mer sets are identical.

Existing primary outputs are not automatically overwritten or skipped. After a failed run, inspect logs and temporary files; remove or move aside only files confirmed to be incomplete. If some accessions have completed, Steps 02 or 03 can be rerun with a manifest containing only the outstanding accessions, provided they belong to the original pooled catalog. Use the complete original manifest again for Steps 04 and 05. Keep the original prefix and results directory. Use a new results directory when changing the genome set or counting parameters.

## Implementation and reproducibility

The pipeline uses Jellyfish canonical counting, exact dictionary lookup, zero filling and ordered column-wise merging. Large text outputs are gzip-compressed; temporary `.jf` databases are removed after successful dump generation unless retained explicitly. The current matching engine is Python; no C++ matching engine is included. The optional C program is a streaming matrix validator only.

The previous development checks used small mock datasets and gzip-compatible tools. Real Jellyfish/pigz execution, population-scale performance, and equivalence against independently generated production matrices should be verified in the target environment before publication. No quantitative speedup or memory reduction is claimed.
