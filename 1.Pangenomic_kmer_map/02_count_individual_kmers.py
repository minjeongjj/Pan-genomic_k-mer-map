#!/usr/bin/env python3
"""Count canonical k-mers and write compressed individual dumps."""
from common import *


def main():
    p=parser(__doc__)
    p.add_argument('--kmer-size',type=int,default=31)
    p.add_argument('--threads',type=int,default=30)
    p.add_argument('--hash-size',default='3G')
    p.add_argument('--jobs',type=int,default=1)
    p.add_argument('--compression-threads',type=int,default=4)
    p.add_argument('--keep-jf',action='store_true')
    a=p.parse_args()
    if min(a.kmer_size,a.threads,a.jobs,a.compression_threads)<1: p.error('Positive numeric options required')
    out=a.results/'02_individual_kmers'; out.mkdir(parents=True,exist_ok=True)
    def worker(row):
        acc,_=row
        jf=out/(acc+'.jellyfish.jf'); dump=out/(acc+'.jellyfish.dump.gz')
        tj,td=Path(str(jf)+'.tmp'),Path(str(dump)+'.tmp')
        fresh(jf,dump,tj,td)
        with (out/(acc+'.log')).open('w') as log:
            start = time.monotonic()
            count_genomes([row],tj,a.kmer_size,a.threads,a.hash_size,log,a.decompression_threads)
            print(f"Count: {time.monotonic()-start:.2f} s",file=log,flush=True)
            start = time.monotonic()
            dump_jf(tj,td,a.compression_threads,log,a.compression_level)
            print(f"Dump: {time.monotonic()-start:.2f} s",file=log,flush=True)
        td.rename(dump)
        if a.keep_jf: tj.rename(jf)
        else: tj.unlink()
        print(f'Completed {acc}',flush=True)
    parallel(worker,genomes(a.input),a.jobs)

if __name__=='__main__': main()
