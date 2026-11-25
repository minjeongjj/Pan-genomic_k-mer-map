#!/usr/bin/env python3
"""Match individual counts in catalog order using the original dictionary algorithm."""
import argparse
import time
from common import decoded, encoded, positive


def build(catalog, individual, accession, output, compression_threads=4, decompression_threads=2, compression_level=1, catalog_header=False):
    start = time.monotonic()
    # Start all external processes before allocating the large dictionary.
    with decoded(catalog,decompression_threads) as src, encoded(output,compression_threads,compression_level) as dst:
        kmer_map = {}
        # Bytes avoid ASCII decode/encode on billions of lines. Values remain exact.
        with decoded(individual,decompression_threads) as f:
            for i,line in enumerate(f,1):
                parts = line.split()
                if parts:
                    kmer_map[parts[0]] = parts[1]
                if i % 50000000 == 0: print(f'Loaded {i:,} lines',flush=True)
        print(f'Loaded {len(kmer_map):,} unique k-mers in {time.monotonic()-start:.2f} s',flush=True)
        start = time.monotonic()
        n = 0
        get = kmer_map.get
        batch = bytearray()
        extend = batch.extend
        dst.write(accession.encode('ascii')+b'\n')
        if catalog_header and src.readline()!=b'Kmer\n': raise ValueError('Invalid catalog header')
        for line in src:
            if catalog_header:
                kmer=line.rstrip(b'\r\n')
                if not kmer: raise ValueError('Empty k-mer in header')
            else:
                parts=line.split(None,1)
                if not parts: continue
                kmer=parts[0]
            extend(get(kmer,b'0')); extend(b'\n')
            n += 1
            if len(batch)>=1024*1024:
                dst.write(batch); batch.clear()
            if n % 100000000 == 0: print(f'Matched {n:,} rows',flush=True)
        if batch: dst.write(batch)
    print(f'Matched {n:,} rows in {time.monotonic()-start:.2f} s',flush=True)
    return n+1

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('catalog','individual','accession','output'): p.add_argument(arg)
    p.add_argument('--compression-threads',type=positive,default=4)
    p.add_argument('--decompression-threads',type=positive,default=2)
    p.add_argument('--compression-level',type=int,choices=range(1,10),default=1)
    p.add_argument("--row-count-file")
    p.add_argument("--catalog-header",action="store_true")
    a=p.parse_args(); n=build(a.catalog,a.individual,a.accession,a.output,a.compression_threads,a.decompression_threads,a.compression_level,a.catalog_header)
    if a.row_count_file:
        from pathlib import Path
        Path(a.row_count_file).write_text(str(n))
