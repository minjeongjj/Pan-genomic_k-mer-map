# Pan-genomic k-mer analysis

This repository contains three analysis modules for constructing pan-genomic *k*-mer maps and investigating trait-associated genomic variation using genome assemblies.

## Analysis modules

| Module | Description |
| --- | --- |
| [1. Pan-genomic k-mer map construction](1.Pangenomic_kmer_map/) | Construct an accession-resolved count matrix of canonical k-mers across genome assemblies. |
| [2. k-mer-based genome-wide association studies (k-mer-GWAS)](2.Kmer_GWAS/) | Genome-wide association analysis using k-mer variation and phenotype data. Implementation and usage documentation are pending. |
| [3. Targeted k-mer screening (TKS)](3.TKS/) | Screen k-mers derived from specified genomic sequences for trait associations, with optional whole-genome copy-number adjustment. |



## Getting started

Run the commands from the respective module directories. See each module's README for input formats, dependencies, execution commands and output descriptions. Module 3 reuses the scripts in Module 1; retain the sibling directory structure. For copy-number-adjusted TKS, prepare the whole-genome k-mer map using Module 1 before running Module 3.

Module 1 and Module 3 include example input formats. Example genome paths are placeholders and must be replaced with actual paths. Set thread counts and parallel job limits according to available CPU, RAM and storage resources.

## Implementation status

Modules 1 and 3 contain analysis scripts. Module 2 currently contains an overview only. End-to-end performance and reproducibility on the full dataset should be verified in the intended computing environment before release.
