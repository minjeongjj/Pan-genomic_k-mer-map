import argparse
import os
import time
import json
import csv
import gzip
import re
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager


def genomes(path):
    path = path.resolve()
    rows, seen = [], set()
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames != ['Accession', 'Genome']:
            raise ValueError('Expected Accession<TAB>Genome header')
        for row in reader:
            acc, name = row['Accession'], row['Genome']
            if None in row or not acc or not name or not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', acc) or acc in seen:
                raise ValueError(f'Invalid manifest row {reader.line_num}')
            genome = (path.parent / Path(name).expanduser()).resolve()
            if not genome.is_file(): raise FileNotFoundError(genome)
            rows.append((acc, genome)); seen.add(acc)
    if not rows: raise ValueError('Empty genome list')
    return rows


def parser(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--input', required=True, type=Path)
    p.add_argument('--results', type=Path, default=Path('results'))
    p.add_argument('--prefix', default='pan_genome')
    p.add_argument('--decompression-threads', type=positive, default=2)
    p.add_argument('--compression-level', type=int, choices=range(1,10), default=1)
    return p


def catalog(a, suffix):
    return a.results / '01_catalog' / (a.prefix + '.cat.jellyfish.' + suffix)


def parallel(fn, rows, jobs):
    if jobs < 1: raise ValueError('Jobs must be positive')
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        return list(pool.map(fn, rows))


def fresh(*paths):
    for p in paths:
        if p.exists(): raise FileExistsError(p)


def compressed(path):
    with Path(path).open('rb') as f: return f.read(2) == b'\x1f\x8b'


@contextmanager
def decoded(path, threads=2):
    # Native decompression, with exit status checked after complete consumption.
    if not compressed(path):
        with Path(path).open('rb') as f: yield f
        return
    with subprocess.Popen(['pigz','-dc','-p',str(threads),'--',str(path)], stdout=subprocess.PIPE) as p:
        try: yield p.stdout
        except BaseException:
            p.stdout.close(); p.terminate(); raise
        finally: p.stdout.close()
        if p.wait(): raise RuntimeError(f'Decompression failed: {path}')


@contextmanager
def encoded(path, threads, level=1):
    if threads < 1: raise ValueError('Compression threads must be positive')
    with Path(path).open('xb') as f, subprocess.Popen(['pigz','-p',str(threads),f'-{level}','-c'], stdin=subprocess.PIPE, stdout=f) as p:
        try: yield p.stdin
        except BaseException:
            p.stdin.close(); p.terminate(); raise
        finally: p.stdin.close()
        if p.wait(): raise RuntimeError(f'Compression failed: {path}')


def line_count(path, threads=2):
    with decoded(path, threads) as f:
        return int(subprocess.check_output(['wc','-l'], stdin=f))


def expected(a):
    return int(catalog(a, 'dump.line').read_text().split()[0]) + 1


def profile(a, acc):
    return a.results / '03_profiles' / (acc + '.row.list.gz')


def feed_genomes(rows, out, decompression_threads=2):
    # Preserve FASTA record boundaries without retaining uncompressed genomes.
    for acc, path in rows:
        with decoded(path, decompression_threads) as src:
            chunk = src.read(8*1024*1024)
            if not chunk.startswith(b'>'): raise ValueError(f'Expected FASTA: {path}')
            while chunk:
                out.write(chunk)
                last = chunk[-1:]
                chunk = src.read(8*1024*1024)
            if last != b'\n': out.write(b'\n')


def count_genomes(rows, jf, k, threads, size, log=None, decompression_threads=2):
    cmd = ['jellyfish','count','-m',str(k),'-t',str(threads),'-C','-s',size,'-o',str(jf),'/dev/stdin']
    if len(rows) == 1 and not compressed(rows[0][1]):
        cmd[-1] = str(rows[0][1])
        subprocess.run(cmd, stdout=log, stderr=log, check=True)
        return
    with subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=log, stderr=log) as p:
        try: feed_genomes(rows,p.stdin,decompression_threads)
        except BaseException:
            p.terminate(); raise
        finally: p.stdin.close()
        if p.wait(): raise RuntimeError('Jellyfish count failed')


def dump_jf(jf, dest, threads, log=None, level=1):
    with encoded(dest,threads,level) as out:
        subprocess.run(['jellyfish','dump','-c',str(jf)],stdout=out,stderr=log,check=True)


def positive(value):
    n = int(value)
    if n < 1: raise argparse.ArgumentTypeError('Must be positive')
    return n


def stamp(path):
    s = Path(path).stat()
    return [s.st_size, s.st_mtime_ns, s.st_ctime_ns, s.st_ino, s.st_dev]


def profile_meta(path):
    return Path(str(path) + '.meta.json')


def timed(label, start):
    print(f'{label}: {time.monotonic()-start:.2f} s', flush=True)


# All streams are checked, including compressor and decompressor exit status.
def checked_pipeline(commands, src=None, dst=None):
    procs = []
    try:
        prev = src
        for i, cmd in enumerate(commands):
            p = subprocess.Popen(cmd, stdin=prev, stdout=dst if i == len(commands)-1 else subprocess.PIPE, env={**os.environ, 'LC_ALL':'C'})
            if procs: procs[-1].stdout.close()
            procs.append(p)
            prev = p.stdout
        codes = [p.wait() for p in reversed(procs)]
        if any(codes): raise RuntimeError(f'Pipeline failed: {codes}')
    finally:
        for p in procs:
            if p.poll() is None: p.terminate()
        for p in procs: p.wait()
