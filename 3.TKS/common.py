"""Streaming matrix I/O and GNU sort helpers for TKS preprocessing."""
import contextlib
import os
import subprocess
import tempfile
from pathlib import Path

@contextlib.contextmanager
def reader(path, threads=2):
    with open(path, 'rb') as f:
        compressed = f.read(2) == b'\x1f\x8b'
        f.seek(0)
        if not compressed:
            yield f
            return
        p = subprocess.Popen(['pigz', '-dc', '-p', str(threads)], stdin=f, stdout=subprocess.PIPE)
        try:
            yield p.stdout
            p.stdout.close()
            if p.wait():
                raise RuntimeError('Decompression failed: {}'.format(path))
        finally:
            p.stdout.close()
            if p.poll() is None:
                p.terminate()
                p.wait()


def header(f):
    fields = f.readline().split()
    if len(fields) < 2 or fields[0].lower() not in (b'kmer', b'#kmer'):
        raise ValueError('Expected Kmer followed by accession IDs in matrix header')
    if len(set(fields[1:])) != len(fields)-1:
        raise ValueError('Duplicate accession IDs')
    return fields[1:]


def sort_body(src, dst, memory, threads, scratch):
    env = dict(os.environ, LC_ALL='C')
    with open(dst, 'wb') as out:
        subprocess.check_call(['sort', '-b', '-k1,1', '-S', memory, '--parallel='+str(threads), '-T', str(scratch), str(src)], stdout=out, env=env)


def new_dir(path, names):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    for name in names:
        if (path / name).exists():
            raise FileExistsError(str(path / name))
    return path


def key_check(key, k):
    if len(key) != k or key.translate(None, b'ACGT'):
        raise ValueError('Invalid k-mer: {!r}'.format(key))
