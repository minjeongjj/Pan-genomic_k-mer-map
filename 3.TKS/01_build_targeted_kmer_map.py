#!/usr/bin/env python3
"""Build a targeted count matrix using the existing Module 1 pipeline."""
import argparse
import os
import subprocess
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', required=True)
    p.add_argument('--results', default='results/01_targeted_map')
    p.add_argument('--prefix', default='targeted')
    p.add_argument('--module1-dir', default=os.environ.get('MODULE1_DIR', str(Path(__file__).resolve().parent.parent / '1.Pangenomic_kmer_map')))
    p.add_argument('--kmer-size', type=int, default=31)
    p.add_argument('--catalog-threads', type=int, default=30)
    p.add_argument('--catalog-hash', default='3G')
    p.add_argument('--individual-threads', type=int, default=10)
    p.add_argument('--individual-hash', default='1G')
    p.add_argument('--count-jobs', type=int, default=1)
    p.add_argument('--profile-jobs', type=int, default=1)
    p.add_argument('--compression-threads', type=int, default=4)
    p.add_argument('--merge-threads', type=int, default=8)
    a = p.parse_args()
    if min(a.kmer_size,a.catalog_threads,a.individual_threads,a.count_jobs,a.profile_jobs,a.compression_threads,a.merge_threads)<1: p.error('Positive sizes and threads required')
    script = Path(a.module1_dir).resolve() / 'Run_pipeline.sh'
    if not script.is_file(): p.error('Module 1 not found: '+str(script))
    env = dict(os.environ, KMER_SIZE=str(a.kmer_size), CATALOG_THREADS=str(a.catalog_threads), CATALOG_HASH=a.catalog_hash, INDIVIDUAL_THREADS=str(a.individual_threads), INDIVIDUAL_HASH=a.individual_hash, COUNT_JOBS=str(a.count_jobs), PROFILE_JOBS=str(a.profile_jobs), CATALOG_COMPRESSION_THREADS=str(a.compression_threads), COMPRESSION_THREADS=str(a.compression_threads), MERGE_THREADS=str(a.merge_threads))
    subprocess.check_call(['bash',str(script),a.input,a.results,a.prefix],env=env)

if __name__=='__main__': main()
