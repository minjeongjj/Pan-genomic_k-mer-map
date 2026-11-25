#!/usr/bin/env python3
"""Build accession-specific count profiles in parallel."""
from common import *
import sys


def main():
    p = parser(__doc__)
    p.add_argument('--jobs', type=int, default=1, help='Concurrent memory-intensive jobs')
    p.add_argument('--compression-threads', type=int, default=4)
    p.add_argument("--catalog-source",choices=["header","dump"],default="header",help="Use the smaller ordered k-mer column, or the original dump")
    a = p.parse_args()
    rows = genomes(a.input)
    cat = catalog(a, 'dump.gz')
    if not cat.is_file(): raise FileNotFoundError(cat)
    source = catalog(a,'header.gz') if a.catalog_source=='header' else cat
    if not source.is_file(): raise FileNotFoundError(source)
    cat_stamp, source_stamp = stamp(cat), stamp(source)
    (a.results / '03_profiles').mkdir(parents=True, exist_ok=True)
    def worker(row):
        acc, _ = row
        src = a.results / '02_individual_kmers' / (acc + '.jellyfish.dump.gz')
        dest = profile(a, acc)
        tmp = Path(str(dest)+'.tmp')
        meta = profile_meta(dest)
        mt = Path(str(meta)+'.tmp')
        fresh(dest, tmp, meta, mt)
        with Path(str(dest)+'.log').open('w') as log:
            subprocess.run([sys.executable, str(Path(__file__).with_name('build_profile.py')), str(source), str(src), acc, str(tmp), "--compression-threads", str(a.compression_threads), "--decompression-threads", str(a.decompression_threads), "--compression-level", str(a.compression_level), "--row-count-file", str(mt)] + (["--catalog-header"] if a.catalog_source=="header" else []), stdout=log, stderr=subprocess.STDOUT, check=True)
        n = int(mt.read_text())
        if n != expected(a): raise ValueError(f"Wrong row count: {acc}: {n}")
        if stamp(cat)!=cat_stamp or stamp(source)!=source_stamp: raise ValueError("Catalog changed during matching")
        tmp.rename(dest)
        mt.write_text(json.dumps({"accession":acc,"rows":n,"stamp":stamp(dest),"catalog_stamp":cat_stamp,"source":a.catalog_source,"source_stamp":source_stamp}))
        mt.rename(meta)
        print(f'Completed {acc}', flush=True)
    parallel(worker, rows, a.jobs)

if __name__ == '__main__': main()
