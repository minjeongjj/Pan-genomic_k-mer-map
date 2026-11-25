#!/usr/bin/env python3
"""Merge compressed profiles with paste and validate during compression."""
from common import *
from contextlib import ExitStack
import tempfile
import shutil


def main():
    p=parser(__doc__)
    p.add_argument('--threads',type=positive,default=150)
    p.add_argument('--validation',choices=['basic','full'],default='basic')
    p.add_argument('--merge-decompression-threads',type=positive,default=1,help='Per input stream; hundreds of inputs already run concurrently')
    p.add_argument('--verify-gzip',action='store_true',help='Optional extra pass to verify final gzip CRC after writing')
    a=p.parse_args(); ids=[acc for acc,_ in genomes(a.input)]; n=expected(a)
    report=a.results/'04_validation'/'profile_validation.tsv'
    with report.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    if [r['Accession'] for r in rows]!=ids or any(r['Status']!='PASS' or int(r['Observed_rows'])!=n for r in rows):
        raise ValueError('Run Step 04 successfully before merging')
    stamp_file=report.with_name('profile_stamps.json')
    if stamp_file.exists():
        saved=json.loads(stamp_file.read_text())
        if any(saved.get(acc)!=stamp(profile(a,acc)) for acc in ids): raise ValueError('Profiles changed since validation; rerun Step 04')
    paths=[catalog(a,'header.gz')]+[profile(a,acc) for acc in ids]
    out=a.results/'05_map'; out.mkdir(parents=True,exist_ok=True)
    dest=out/(a.prefix+'.kmer_merged.table.tsv.gz'); tmp=Path(str(dest)+'.tmp'); fresh(dest,tmp)
    start=time.monotonic(); header='\t'.join(['Kmer']+ids)
    with tempfile.TemporaryDirectory() as work, ExitStack() as stack:
        # A native validator is optional. The portable awk fallback also streams.
        cc=shutil.which('cc') or shutil.which('gcc')
        exe=Path(work)/'validate_stream'
        if cc:
            subprocess.run([cc,'-O3','-std=c99',str(Path(__file__).with_name('validate_stream.c')),'-o',str(exe)],check=True)
            validator=[str(exe),str(n),str(len(ids)+1),header]
        else:
            script='BEGIN{FS="\\t"} NR==1 && $0!=header {bad=1} NF!=cols || $0 ~ /^\\t|\\t\\t|\\t$/ {bad=1} {print} END{if(bad || NR!=rows) exit 1}'
            validator=['awk','-v',f'header={header}','-v',f'cols={len(ids)+1}','-v',f'rows={n}',script]
        streams=[stack.enter_context(decoded(path,a.merge_decompression_threads)) for path in paths]
        fds=[f.fileno() for f in streams]
        with encoded(tmp,a.threads,a.compression_level) as dst:
            with subprocess.Popen(['paste',*[f'/dev/fd/{fd}' for fd in fds]],pass_fds=tuple(fds),stdout=subprocess.PIPE) as paste:
                try:
                    subprocess.run(validator,stdin=paste.stdout,stdout=dst,check=True,env={**os.environ,'LC_ALL':'C'})
                finally: paste.stdout.close()
                if paste.wait(): raise RuntimeError('paste failed')
    timed('Merge, inline validation and compression',start)
    if a.validation=='full':
        start=time.monotonic()
        with decoded(tmp,a.decompression_threads) as merged, decoded(catalog(a,'dump.gz'),a.decompression_threads) as cat:
            next(merged)
            for i,line in enumerate(cat,2):
                kmer,total=line.split(); fields=next(merged).rstrip(b'\n').split(b'\t')
                if fields[0]!=kmer or any(not x.isdigit() for x in fields[1:]) or sum(map(int,fields[1:]))!=int(total):
                    raise ValueError(f'Count or order mismatch at row {i}')
            if next(merged,None) is not None: raise ValueError('Extra merged rows')
        timed('Full validation',start)
    elif a.verify_gzip:
        start=time.monotonic(); subprocess.run(['pigz','-t','-p',str(a.decompression_threads),str(tmp)],check=True); timed('Final gzip integrity scan',start)
    tmp.rename(dest)
    integrity='readback' if a.validation=='full' or a.verify_gzip else 'compressor_exit_only'
    (out/'map_validation.tsv').write_text(f'Mode\tRows\tColumns\tGzip_check\tStatus\n{a.validation}\t{n}\t{len(ids)+1}\t{integrity}\tPASS\n')
    print(f'Completed: {dest}',flush=True)

if __name__=='__main__': main()
