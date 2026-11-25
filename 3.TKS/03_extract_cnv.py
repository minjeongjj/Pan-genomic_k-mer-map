#!/usr/bin/env python3
"""Extract exact targeted k-mer matches from a whole-genome map in one scan."""
import argparse
import tempfile
import time
from pathlib import Path
from common import reader, header, sort_body, new_dir, key_check


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pav',required=True)
    p.add_argument('--pav-header',required=True)
    p.add_argument('--genome-map',required=True)
    p.add_argument('--out-dir',default='results/03_cnv')
    p.add_argument('--kmer-size',type=int,default=31)
    p.add_argument('--sort-memory',default='2G')
    p.add_argument('--sort-threads',type=int,default=4)
    p.add_argument('--decompression-threads',type=int,default=2)
    p.add_argument('--sort-tmpdir')
    p.add_argument('--allow-missing',action='store_true',help='Report missing keys and allow original engine to skip them')
    a=p.parse_args()
    if min(a.kmer_size,a.sort_threads,a.decompression_threads)<1: p.error('Positive sizes/threads required')
    start=time.time()
    out=new_dir(a.out_dir,['CNV_table.tsv','CNV_header.tsv','CNV_summary.tsv','missing_kmers.txt'])
    scratch=Path(a.sort_tmpdir or out); scratch.mkdir(parents=True,exist_ok=True)
    with open(a.pav_header,'rb') as f: pav_acc=header(f)
    # Start decompressor before allocating the target key set (avoid large-process fork).
    with tempfile.TemporaryDirectory(dir=str(scratch)) as td:
        raw=Path(td)/'raw'; sorted_path=Path(td)/'sorted'
        with reader(a.genome_map,a.decompression_threads) as f, open(raw,'wb') as w:
            genome_acc=header(f); idx={x:i+1 for i,x in enumerate(genome_acc)}
            missing_ids=set(pav_acc)-set(genome_acc)
            if missing_ids: raise ValueError('PAV accessions absent from genome-wide map: '+str(sorted(missing_ids)))
            # Preserve every PAV accession for the original CNV consistency checks.
            cols=[idx[x] for x in pav_acc]
            same=genome_acc==pav_acc
            wanted=set()
            with open(a.pav,'rb') as pf:
                for line in pf:
                    key=line.split(None,1)[0]; key_check(key,a.kmer_size); wanted.add(key)
            if not wanted: raise ValueError('Empty PAV table')
            print('Loaded {:,} target keys; scanning genome-wide map'.format(len(wanted)),flush=True)
            scanned=hits=0
            for line in f:
                pos=line.find(b'\t'); key=line[:pos]
                if pos<0: raise ValueError('Genome map must be tab-separated')
                scanned+=1
                if key in wanted:
                    fields=line.split()
                    if len(fields)!=len(genome_acc)+1: raise ValueError('Invalid genome-wide row for '+str(key))
                    selected=fields[1:] if same else [fields[i] for i in cols]
                    if any(not x.isdigit() for x in selected): raise ValueError('Invalid CNV count for '+str(key))
                    w.write(key+b'\t'+b'\t'.join(selected)+b'\n'); hits+=1
                if scanned%100000000==0: print('[Progress] {:,} scanned; {:,} hits; {:.1f}s'.format(scanned,hits,time.time()-start),flush=True)
        del wanted
        sort_body(raw,sorted_path,a.sort_memory,a.sort_threads,scratch)
        # Merge-compare sorted keys, detect duplicates and enumerate missing keys without another large set.
        unique=missing=0; prev=None
        with open(sorted_path,'rb') as cf, open(a.pav,'rb') as pf, open(out/'CNV_table.tsv.tmp','wb') as w, open(out/'missing_kmers.txt','wb') as mf:
            pl=next(pf,None)
            for line in cf:
                key=line.partition(b'\t')[0]
                if prev is not None and key==prev.partition(b'\t')[0]:
                    if line!=prev: raise ValueError('Conflicting CNV duplicate: '+str(key))
                    continue
                while pl is not None and pl.partition(b'\t')[0]<key:
                    mf.write(pl.partition(b'\t')[0]+b'\n'); missing+=1; pl=next(pf,None)
                if pl is None or pl.partition(b'\t')[0]!=key: raise ValueError('PAV must be sorted and unique')
                pl=next(pf,None); w.write(line); prev=line; unique+=1
            while pl is not None:
                mf.write(pl.partition(b'\t')[0]+b'\n'); missing+=1; pl=next(pf,None)
        (out/'CNV_summary.tsv').write_text('metric\tvalue\ngenome_rows_scanned\t{}\nmatched_unique_kmers\t{}\nmissing_kmers\t{}\naccessions\t{}\n'.format(scanned,unique,missing,len(pav_acc)))
        if missing and not a.allow_missing:
            raise ValueError('{} targeted k-mers absent from genome map; inspect missing_kmers.txt; use a new output directory to retry'.format(missing))
        (out/'CNV_table.tsv.tmp').rename(out/'CNV_table.tsv')
        (out/'CNV_header.tsv').write_bytes(b'Kmer\t'+b'\t'.join(pav_acc)+b'\n')
    print('Extracted {:,} CNV rows; {:,} missing; {:.1f}s'.format(unique,missing,time.time()-start),flush=True)

if __name__=='__main__': main()
