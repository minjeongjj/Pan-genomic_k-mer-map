#!/usr/bin/env python3
"""Prepare a headerless sorted targeted count matrix; retain original counts."""
import argparse
import tempfile
import time
from pathlib import Path
from common import reader, header, sort_body, new_dir, key_check
from importlib import import_module


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--map', required=True)
    p.add_argument('--phenotype', required=True)
    p.add_argument('--target', required=True)
    p.add_argument('--out-dir', default='results/02_pav')
    p.add_argument('--kmer-size', type=int, default=31)
    p.add_argument('--sort-memory', default='2G')
    p.add_argument('--sort-threads', type=int, default=4)
    p.add_argument('--decompression-threads', type=int, default=2)
    p.add_argument('--sort-tmpdir')
    a = p.parse_args()
    if min(a.kmer_size,a.sort_threads,a.decompression_threads)<1: p.error('Positive sizes/threads required')
    start = time.time()
    # Reuse the original engine's phenotype and target parsers.
    engine = import_module('04_targeted_kmer_screening')
    phen = engine.read_phenotype(a.phenotype)
    target = engine.read_names(a.target)
    out = new_dir(a.out_dir, ['PAV_table.tsv','PAV_header.tsv','PAV_summary.tsv'])
    scratch = Path(a.sort_tmpdir or out); scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(scratch)) as td:
        raw = Path(td)/'raw'; sorted_path = Path(td)/'sorted'
        with reader(a.map,a.decompression_threads) as f, open(raw,'wb') as w:
            acc = header(f); ids = {x.decode('utf-8') for x in acc}
            missing = set(target)-ids
            if missing: raise ValueError('Target IDs absent from PAV: '+','.join(sorted(missing)))
            matched = len(set(phen)&ids)
            if matched<4: raise ValueError('Fewer than four phenotyped PAV accessions')
            n=0
            for line in f:
                fields=line.split()
                if len(fields)!=len(acc)+1: raise ValueError('Invalid PAV row {}'.format(n+1))
                key_check(fields[0],a.kmer_size)
                if any(not x.isdigit() for x in fields[1:]): raise ValueError('Invalid PAV count at '+str(fields[0]))
                w.write(b'\t'.join(fields)+b'\n'); n+=1
            if not n: raise ValueError('Empty targeted map')
        sort_body(raw,sorted_path,a.sort_memory,a.sort_threads,scratch)
        # Reject conflicting duplicates, collapse identical ones as in the original engine.
        unique=0; prev=None
        with open(sorted_path,'rb') as f, open(out/'PAV_table.tsv.tmp','wb') as w:
            for line in f:
                key=line.partition(b'\t')[0]
                if prev is not None and key==prev.partition(b'\t')[0]:
                    if line!=prev: raise ValueError('Conflicting PAV duplicate: '+str(key))
                    continue
                w.write(line); prev=line; unique+=1
        (out/'PAV_table.tsv.tmp').rename(out/'PAV_table.tsv')
        (out/'PAV_header.tsv').write_bytes(b'Kmer\t'+b'\t'.join(acc)+b'\n')
        (out/'PAV_summary.tsv').write_text('metric\tvalue\ninput_rows\t{}\nunique_kmers\t{}\naccessions\t{}\nphenotyped_accessions\t{}\ntarget_accessions\t{}\nphenotype_ids_absent\t{}\n'.format(n,unique,len(acc),matched,len(target),len(set(phen)-ids)))
    print('Prepared {:,} PAV k-mers in {:.1f}s'.format(unique,time.time()-start),flush=True)

if __name__=='__main__': main()
