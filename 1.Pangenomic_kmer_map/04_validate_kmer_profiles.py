#!/usr/bin/env python3
"""Validate profiles using creation metadata, or rescan compressed data."""
from common import *


def main():
    p=parser(__doc__)
    p.add_argument('--jobs',type=positive,default=1)
    p.add_argument('--rescan',action='store_true',help='Read all profiles even when matching creation metadata exists')
    a=p.parse_args(); n=expected(a); cat=catalog(a,'dump.gz')
    def worker(row):
        acc,_=row; path=profile(a,acc); meta=profile_meta(path)
        if not path.is_file(): return (acc,'NA',n,'missing','FAIL')
        current=stamp(path)
        if meta.is_file() and not a.rescan:
            m=json.loads(meta.read_text())
            if m.get('source') in ('header','dump') and m.get('source_stamp')==stamp(catalog(a,'header.gz' if m['source']=='header' else 'dump.gz')) and m.get('stamp')==current and m.get('catalog_stamp')==stamp(cat) and m.get('accession')==acc and m.get('rows')==n:
                return (acc,n,n,'creation_metadata','PASS')
        # Validate legacy profiles in one pass, including single nonnegative integer field.
        with decoded(path,a.decompression_threads) as src:
            script='NR==1 {if ($0!=acc) bad=1; next} $0 !~ /^[0-9]+$/ {bad=1} END {printf "%.0f\\n", NR; if(bad || NR!=n) exit 1}'
            r=subprocess.run(['awk','-v',f'acc={acc}','-v',f'n={n}',script],stdin=src,stdout=subprocess.PIPE,env={**os.environ,'LC_ALL':'C'})
        count=int(r.stdout.strip()); ok=r.returncode==0 and current==stamp(path)
        return (acc,count,n,'rescan','PASS' if ok else 'FAIL')
    results=parallel(worker,genomes(a.input),a.jobs)
    out=a.results/'04_validation'; out.mkdir(parents=True,exist_ok=True)
    with (out/'profile_validation.tsv').open('w') as f:
        f.write('Accession\tObserved_rows\tExpected_rows\tMethod\tStatus\n')
        for row in results: f.write('\t'.join(map(str,row))+'\n')
    if any(r[-1]!='PASS' for r in results): raise SystemExit('Profile validation failed')
    (out/'profile_stamps.json').write_text(json.dumps({acc:stamp(profile(a,acc)) for acc,_ in genomes(a.input)}))
    print('All profiles passed',flush=True)

if __name__=='__main__': main()
