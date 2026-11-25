#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
input=${1:?Usage: bash Run_pipeline.sh genome_list.tsv [results] [prefix]}
results=${2:-results}
prefix=${3:-pan_genome}
args=(--input "$input" --results "$results" --prefix "$prefix" --compression-level "${COMPRESSION_LEVEL:-1}" --decompression-threads "${DECOMPRESSION_THREADS:-2}")
k=${KMER_SIZE:-31}
start=${START_STEP:-1}
stop=${STOP_STEP:-5}
if ! [[ $start =~ ^[1-5]$ && $stop =~ ^[1-5]$ ]] || (( start > stop )); then echo "Invalid START_STEP/STOP_STEP" >&2; exit 2; fi
if (( start <= 1 && stop >= 1 )); then
    time python3 "$script_dir/01_build_kmer_catalog.py" "${args[@]}" --kmer-size "$k" --threads "${CATALOG_THREADS:-100}" --hash-size "${CATALOG_HASH:-15G}" --compression-threads "${CATALOG_COMPRESSION_THREADS:-8}"
fi
if (( start <= 2 && stop >= 2 )); then
    time python3 "$script_dir/02_count_individual_kmers.py" "${args[@]}" --kmer-size "$k" --threads "${INDIVIDUAL_THREADS:-30}" --hash-size "${INDIVIDUAL_HASH:-3G}" --jobs "${COUNT_JOBS:-1}" --compression-threads "${COMPRESSION_THREADS:-4}"
fi
if (( start <= 3 && stop >= 3 )); then
    time python3 "$script_dir/03_build_kmer_profiles.py" "${args[@]}" --jobs "${PROFILE_JOBS:-1}" --compression-threads "${COMPRESSION_THREADS:-4}" --catalog-source "${CATALOG_SOURCE:-header}"
fi
if (( start <= 4 && stop >= 4 )); then
    time python3 "$script_dir/04_validate_kmer_profiles.py" "${args[@]}" --jobs "${VALIDATION_JOBS:-1}"
fi
if (( start <= 5 && stop >= 5 )); then
    time python3 "$script_dir/05_merge_kmer_profiles.py" "${args[@]}" --threads "${MERGE_THREADS:-150}" --validation "${VALIDATION_MODE:-basic}" --merge-decompression-threads "${MERGE_DECOMPRESSION_THREADS:-1}"
fi
