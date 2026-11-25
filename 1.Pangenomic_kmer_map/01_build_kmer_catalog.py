#!/usr/bin/env python3
"""Build the compressed catalog, ordered header and row count in one dump pass."""
from common import *


def main():
    p = parser(__doc__)
    p.add_argument('--output-dir', type=Path)
    p.add_argument('--kmer-size', type=positive, default=31)
    p.add_argument('--threads', type=positive, default=100)
    p.add_argument('--hash-size', default='15G')
    p.add_argument('--compression-threads', type=positive, default=8, help='Workers per compressor; two run concurrently')
    p.add_argument('--keep-jf', action='store_true')
    a = p.parse_args()
    rows = genomes(a.input)
    out = a.output_dir or a.results/'01_catalog'; out.mkdir(parents=True, exist_ok=True)
    base = out/(a.prefix+'.cat.jellyfish.')
    jf,dump,header,lines = [Path(str(base)+s) for s in ('jf','dump.gz','header.gz','dump.line')]
    tj,td,th,tl = [Path(str(x)+'.tmp') for x in (jf,dump,header,lines)]
    fresh(jf,dump,header,lines,tj,td,th,tl)
    start = time.monotonic()
    count_genomes(rows,tj,a.kmer_size,a.threads,a.hash_size,decompression_threads=a.decompression_threads)
    timed('Jellyfish count',start)
    start = time.monotonic()
    # tee sends the original dump to pigz and the same ordered stream to awk.
    # awk writes the header and counts rows without reading compressed output again.
    with encoded(td,a.compression_threads,a.compression_level) as d, encoded(th,a.compression_threads,a.compression_level) as h:
        h.write(b'Kmer\n'); h.flush()
        procs = []
        try:
            j = subprocess.Popen(['jellyfish','dump','-c',str(tj)],stdout=subprocess.PIPE); procs.append(j)
            t = subprocess.Popen(['tee','/dev/fd/'+str(d.fileno())],stdin=j.stdout,stdout=subprocess.PIPE,pass_fds=(d.fileno(),)); procs.append(t); j.stdout.close()
            with tl.open('wb') as countfile:
                w = subprocess.Popen(['awk','{print $1} END {printf "%.0f\\n", NR > "/dev/stderr"}'],stdin=t.stdout,stdout=h,stderr=countfile,env={**os.environ,'LC_ALL':'C'}); procs.append(w); t.stdout.close()
                codes=[x.wait() for x in reversed(procs)]
                if any(codes): raise RuntimeError(f'Dump pipeline failed: {codes}')
        finally:
            for x in procs:
                if x.poll() is None: x.terminate()
            for x in procs: x.wait()
    n = int(tl.read_text().strip())
    if n < 1: raise ValueError('Empty catalog')
    tl.write_text(f'{n}\t{dump}\n')
    for x,y in [(td,dump),(th,header),(tl,lines)]: x.rename(y)
    if a.keep_jf: tj.rename(jf)
    else: tj.unlink()
    timed('Dump, header and row count',start)
    print(f'Completed: {n:,} unique k-mers',flush=True)

if __name__=='__main__': main()
