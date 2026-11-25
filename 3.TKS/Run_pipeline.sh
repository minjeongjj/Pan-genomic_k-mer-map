#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
input=${1:?Usage: bash Run_pipeline.sh targeted_genome_list.tsv genome_map.tsv.gz phenotype.tsv target_accessions.txt [results]}
genome=${2:?Genome-wide map required (use - for none mode)}
phenotype=${3:?Phenotype required}
target=${4:?Target list required}
results=${5:-results}
start=${TKS_START_STEP:-1}; stop=${TKS_STOP_STEP:-4}
mode=${COPY_ADJUST:-methods}
[[ $start =~ ^[1-4]$ && $stop =~ ^[1-4]$ ]] && (( start <= stop )) || { echo 'Invalid TKS step range' >&2; exit 2; }
[[ $mode == methods || $mode == none ]] || { echo 'Invalid COPY_ADJUST' >&2; exit 2; }
if [[ $mode == methods && ! -f $genome ]]; then echo 'Genome-wide map not found' >&2; exit 2; fi
[[ -f $phenotype && -f $target ]] || { echo 'Phenotype or target list not found' >&2; exit 2; }
k=${KMER_SIZE:-31}
pavmap=${TARGETED_MAP:-$results/01_targeted_map/05_map/targeted.kmer_merged.table.tsv.gz}
opts=(--kmer-size "$k" --sort-memory "${SORT_MEMORY:-2G}" --sort-threads "${SORT_THREADS:-4}" --decompression-threads "${DECOMPRESSION_THREADS:-2}")
if [[ -n ${SORT_TMPDIR:-} ]]; then opts+=(--sort-tmpdir "$SORT_TMPDIR"); fi
if (( start <= 1 && stop >= 1 )); then
 time python3 "$here/01_build_targeted_kmer_map.py" --input "$input" --results "$results/01_targeted_map" --kmer-size "$k" --catalog-threads "${CATALOG_THREADS:-30}" --catalog-hash "${CATALOG_HASH:-3G}" --individual-threads "${INDIVIDUAL_THREADS:-10}" --individual-hash "${INDIVIDUAL_HASH:-1G}" --count-jobs "${COUNT_JOBS:-1}" --profile-jobs "${PROFILE_JOBS:-1}" --compression-threads "${COMPRESSION_THREADS:-4}" --merge-threads "${MERGE_THREADS:-8}"
fi
if (( start <= 2 && stop >= 2 )); then
 time python3 "$here/02_prepare_pav.py" --map "$pavmap" --phenotype "$phenotype" --target "$target" --out-dir "$results/02_pav" "${opts[@]}"
fi
cnv=()
if [[ $mode == methods ]]; then
 if (( start <= 3 && stop >= 3 )); then
  time python3 "$here/03_extract_cnv.py" --pav "$results/02_pav/PAV_table.tsv" --pav-header "$results/02_pav/PAV_header.tsv" --genome-map "$genome" --out-dir "$results/03_cnv" "${opts[@]}"
 fi
 cnv=(--cnv "$results/03_cnv/CNV_table.tsv" --cnv-header "$results/03_cnv/CNV_header.tsv")
fi
if (( start <= 4 && stop >= 4 )); then
 time python3 "$here/04_targeted_kmer_screening.py" --phenotype "$phenotype" --resistant "$target" --pav "$results/02_pav/PAV_table.tsv" --pav-header "$results/02_pav/PAV_header.tsv" "${cnv[@]}" --copy-adjust "$mode" --weight-scope "${WEIGHT_SCOPE:-all}" --min-resistant "${MIN_TARGET:-3}" --max-other "${MAX_OTHER:-1}" --alpha "${ALPHA:-0.05}" --already-sorted --out "$results/04_screening/TKS.${OUTPUT_FORMAT:-tsv}"
fi
