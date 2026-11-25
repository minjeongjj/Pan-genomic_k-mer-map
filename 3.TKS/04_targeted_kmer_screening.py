#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TKS: NLR PAV specificity -> original DSI Welch test, optionally CNV Methods -> BH.

PAV counts (NLR): >=1 presence; 0 absence. 301 accession columns.
CNV counts (genome-wide): used ONLY for copy-number weighting.
--weight-scope all: use all phenotyped PAV carriers (default).
--weight-scope target: use only phenotyped PAV carriers from the target list.
The target list always defines the PAV specificity groups.

Methods: w = number of weight-scope PAV carriers / sum(their genomic CNV counts).
Adjusted carrier mean = w*mean1+(1-w)*mean0
Adjusted carrier variance = w*var1+(1-w)*var0+w*(1-w)*(mean1-mean0)**2.
Welch: n1/n0 are PAV carrier/noncarrier sample counts, not copies.

BH runs on valid k-mers after specificity/sample/CNV validation filters.
"""
import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

MISSING = {'NA', 'NAN', 'NULL', '-', '.'}
FIELDS = ['kmer','mean_present','mean_absent','present_count','absent_count',
          't_stat','p_value','p_adj','significant','w','m','copy_sum',
          'resistant_present','other_present']
RAW_FIELDS = [x for x in FIELDS if x not in ('p_adj','significant')]
VALID_FIELDS = ['kmer','accession','pav_count','cnv_count','reason']


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--phenotype', required=True, help='Fixed original DSI: accession value')
    ap.add_argument('--pav', required=True, help='NLR k-mer count table')
    ap.add_argument('--pav-header', required=True, help='PAV accession header')
    ap.add_argument('--cnv', help='Genome-wide k-mer count table (required for methods)')
    ap.add_argument('--cnv-header', help='CNV accession header (required for methods)')
    ap.add_argument('--resistant', required=True, help='Target/resistant accessions for specificity only')
    ap.add_argument('--out', required=True, help='Output path (e.g. mj.tsv or mj.csv)')
    ap.add_argument('--output-format', choices=('tsv', 'csv'), default=None,
                    help='Output delimiter/extension; default: infer from --out extension, otherwise tsv')
    ap.add_argument('--copy-adjust', choices=('none','methods'), default='methods')
    ap.add_argument('--weight-scope', choices=('all','target'), default='all',
                    help='methods weighting: all PAV carriers (default) or target PAV carriers only')
    ap.add_argument('--min-resistant', type=int, default=3)
    ap.add_argument('--max-other', type=int, default=1)
    ap.add_argument('--alpha', type=float, default=0.05)
    ap.add_argument('--sort-memory', default=os.environ.get('SORT_MEMORY','2G'))
    ap.add_argument('--sort-tmpdir', default=os.environ.get('SORT_TMPDIR'))
    ap.add_argument('--already-sorted', action='store_true', help='Skip GNU sort; inputs must be sorted by k-mer')
    ap.add_argument('--debug-kmers', help='File with k-mer IDs to debug (one per line)')
    args = ap.parse_args()
    if args.copy_adjust == 'methods' and not (args.cnv and args.cnv_header):
        ap.error('--cnv and --cnv-header are required for methods')
    if not 0 < args.alpha < 1 or args.min_resistant < 0 or args.max_other < 0:
        ap.error('Invalid alpha/specificity threshold')
    # Select an actual delimiter, not merely a filename extension.
    root, ext = os.path.splitext(args.out)
    ext = ext.lower()
    if ext not in ('', '.tsv', '.csv'):
        ap.error('--out must end with .tsv or .csv (or have no extension)')
    args.output_format = args.output_format or (ext[1:] if ext else 'tsv')
    if ext:
        args.out = root + '.' + args.output_format
    else:
        args.out = args.out + '.' + args.output_format
    args.delimiter = '\t' if args.output_format == 'tsv' else ','
    return args


def read_header(path):
    with open(path, encoding='utf-8') as f:
        for line in f:
            p = line.split()
            if not p: continue
            accs = p[1:] if p[0].startswith('#') or p[0].lower() == 'kmer' else p
            if not accs or len(accs) != len(set(accs)):
                raise ValueError('Empty or duplicated accession header: '+path)
            return accs
    raise ValueError('Empty header: '+path)


def read_names(path):
    ids = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            p=line.split()
            if p and not p[0].startswith('#'):
                ids.append(p[0])
    if not ids or len(ids)!=len(set(ids)):
        raise ValueError('Empty or duplicated ID list: '+path)
    return ids


def read_phenotype(path):
    data = {}
    with open(path, encoding='utf-8') as f:
        for ln,line in enumerate(f,1):
            p=line.split()
            if not p or p[0].startswith('#') or len(p)<2: continue
            if p[1].upper() in MISSING: continue
            try: y=float(p[1])
            except ValueError:
                if ln==1: continue
                raise ValueError('Bad DSI at line {}'.format(ln))
            if not math.isfinite(y): continue
            if p[0] in data: raise ValueError('Duplicate phenotype ID: '+p[0])
            data[p[0]]=y
    if not data: raise ValueError('No numeric phenotype values')
    return data


def sort_copy(src, dst, tmpdir, memory):
    env=dict(os.environ); env['LC_ALL']='C'
    with open(dst,'wb') as out:
        subprocess.check_call(['sort','-b','-k1,1','-S',memory,'-T',tmpdir,src],stdout=out,env=env)


def next_row(f):
    for line in f:
        if line.strip() and not line.lstrip().startswith('#'):
            return line.split()
    return None


def unique_rows(f, source, counts):
    """Yield sorted k-mer rows; safely collapse exactly identical duplicates.

    If duplicate k-mers have differing counts, fail instead of silently
    choosing or combining biologically different genotype vectors.
    """
    previous = None
    previous_key = None
    while True:
        row = next_row(f)
        if row is None:
            break
        key = row[0]
        if previous_key is not None:
            if key < previous_key:
                raise ValueError('Unsorted {} key: previous={} current={}'.format(
                    source, previous_key, key))
            if key == previous_key:
                # Whitespace variations are normalized by line.split().
                if row == previous:
                    counts['duplicate_{}_identical'.format(source.lower())] += 1
                    continue
                raise ValueError('Conflicting {} rows for k-mer {}: '
                                 'same k-mer has different accession counts'.format(source, key))
        previous_key = key
        previous = row
        yield row


def parse_count(s):
    if s.upper() in MISSING: return None
    n=int(s)
    if n<0: raise ValueError('Negative count')
    return n


def compute_bh(raw_file, all_file, sig_file, alpha, delimiter):
    pvals=[]
    with open(raw_file, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f, delimiter=delimiter): pvals.append(float(row['p_value']))
    if pvals:
        reject,qvals,_,_=multipletests(np.asarray(pvals,dtype=float), alpha=alpha, method='fdr_bh')
    else:
        reject=[]; qvals=[]
    num=0
    with open(raw_file,newline='',encoding='utf-8') as rf, \
         open(all_file,'w',newline='',encoding='utf-8') as af, \
         open(sig_file,'w',newline='',encoding='utf-8') as sf:
        wa=csv.DictWriter(af,fieldnames=FIELDS,delimiter=delimiter); ws=csv.DictWriter(sf,fieldnames=FIELDS,delimiter=delimiter)
        wa.writeheader(); ws.writeheader()
        for i,row in enumerate(csv.DictReader(rf, delimiter=delimiter)):
            row['p_adj']='{:.17g}'.format(float(qvals[i]))
            row['significant']='1' if reject[i] else '0'
            wa.writerow(row)
            if reject[i]: ws.writerow(row); num+=1
    return num


def main():
    args=parse_args(); start=time.time()
    phenotype=read_phenotype(args.phenotype)
    pav_accs=read_header(args.pav_header)
    pav_idx={a:i+1 for i,a in enumerate(pav_accs)}
    target=read_names(args.resistant); target_set=set(target)
    missing_target=target_set-set(pav_accs)
    if missing_target: raise ValueError('Target IDs missing from PAV header: '+','.join(sorted(missing_target)))
    target_indices=[pav_idx[a] for a in target]
    other_indices=[pav_idx[a] for a in pav_accs if a not in target_set]
    pheno_not_pav=set(phenotype)-set(pav_accs)
    if pheno_not_pav:
        print('WARNING: {} phenotype IDs absent from PAV header'.format(len(pheno_not_pav)),file=sys.stderr)
    samples=[(pav_idx[a],a,y) for a,y in phenotype.items() if a in pav_idx]
    if len(samples)<4: raise ValueError('Too few phenotyped accessions matched to PAV')
    if args.copy_adjust=='methods':
        cnv_accs=read_header(args.cnv_header)
        cnv_idx={a:i+1 for i,a in enumerate(cnv_accs)}
        required_accs={a for _,a,_ in samples} if args.weight_scope=='all' else ({a for _,a,_ in samples} & target_set)
        missing=required_accs-set(cnv_accs)
        if missing:
            raise ValueError('Methods mode with --weight-scope {} requires CNV header columns for each possible PAV carrier in the selected scope. Missing {} IDs; examples: {}'.format(args.weight_scope,len(missing),','.join(sorted(missing)[:12])))
    else: cnv_accs=[]; cnv_idx={}
    debug=set(read_names(args.debug_kmers)) if args.debug_kmers else set()
    counts=dict(total=0,specificity_fail=0,invalid_pav=0,insufficient=0,
                missing_cnv_kmer=0,invalid_cnv=0,cnv_less_than_pav=0,
                cnv_zero_for_carrier=0,cnv_checked=0,invalid_test=0,tested=0,
                duplicate_pav_identical=0,duplicate_cnv_identical=0)
    print('Phenotypes: {} / PAV: {} / target: {} / other: {} / CNV: {}'.format(
        len(samples),len(pav_accs),len(target),len(other_indices),len(cnv_accs)),flush=True)
    print('Mode={} | Weight scope={} | Specificity target >= {}, other <= {} | MAF disabled | BH alpha={} | Output={}'.format(
        args.copy_adjust,args.weight_scope,args.min_resistant,args.max_other,args.alpha,args.output_format),flush=True)
    outdir=os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir,exist_ok=True)
    scratch=args.sort_tmpdir or outdir
    os.makedirs(scratch,exist_ok=True)
    stem = os.path.splitext(args.out)[0]
    all_file = stem + '.all.' + args.output_format
    validation_file = stem + '.validation.' + args.output_format
    with tempfile.TemporaryDirectory(prefix='tks_pav_cnv_',dir=scratch) as td:
        pp=args.pav; cp=args.cnv
        if not args.already_sorted:
            print('Sorting PAV...',flush=True)
            pp=os.path.join(td,'pav.sorted'); sort_copy(args.pav,pp,td,args.sort_memory)
            if args.copy_adjust=='methods':
                print('Sorting CNV...',flush=True)
                cp=os.path.join(td,'cnv.sorted'); sort_copy(args.cnv,cp,td,args.sort_memory)
        raw=os.path.join(td,'raw.' + args.output_format)
        with open(pp,encoding='utf-8') as pf, \
             open(raw,'w',newline='',encoding='utf-8') as rawout, \
             open(validation_file,'w',newline='',encoding='utf-8') as vf:
            wr=csv.DictWriter(rawout,fieldnames=RAW_FIELDS,delimiter=args.delimiter); wr.writeheader()
            vw=csv.DictWriter(vf,fieldnames=VALID_FIELDS,delimiter=args.delimiter); vw.writeheader()
            cf=open(cp,encoding='utf-8') if args.copy_adjust=='methods' else None
            pav_rows=unique_rows(pf,'PAV',counts)
            cnv_rows=unique_rows(cf,'CNV',counts) if cf else None
            crow=next(cnv_rows,None) if cnv_rows else None
            try:
                while True:
                    pav=next(pav_rows,None)
                    if pav is None: break
                    key=pav[0]; counts['total']+=1
                    if counts['total']%1000000==0:
                        print('[Progress] {:,} PAV | {:,} tested | {:.0f}s'.format(counts['total'],counts['tested'],time.time()-start),flush=True)
                    if len(pav)!=len(pav_accs)+1:
                        counts['invalid_pav']+=1; continue
                    try:
                        pv=[parse_count(v) for v in pav[1:]]
                    except ValueError:
                        counts['invalid_pav']+=1; continue
                    if any(v is None for v in pv):
                        counts['invalid_pav']+=1; continue
                    nr=sum(pv[i-1]>=1 for i in target_indices)
                    if nr<args.min_resistant:
                        counts['specificity_fail']+=1; continue
                    no=sum(pv[i-1]>=1 for i in other_indices)
                    if no>args.max_other:
                        counts['specificity_fail']+=1; continue
                    carriers=[(a,y,pv[idx-1]) for idx,a,y in samples if pv[idx-1]>=1]
                    absent=[y for idx,a,y in samples if pv[idx-1]==0]
                    n1=len(carriers); n0=len(absent)
                    if n1<2 or n0<2:
                        counts['insufficient']+=1; continue
                    if args.copy_adjust=='methods':
                        while crow is not None and crow[0]<key:
                            crow=next(cnv_rows,None)
                        if crow is None or crow[0]!=key:
                            counts['missing_cnv_kmer']+=1; continue
                        current=crow; crow=next(cnv_rows,None)
                        if len(current)!=len(cnv_accs)+1:
                            counts['invalid_cnv']+=1; continue
                        try:
                            cnv_values={acc:parse_count(current[cnv_idx[acc]]) for acc in pav_accs if acc in cnv_idx}
                        except (ValueError,IndexError):
                            counts['invalid_cnv']+=1; continue
                        for acc,pavcount in zip(pav_accs,pv):
                            if acc not in cnv_values: continue
                            c=cnv_values[acc]
                            if c is None: continue
                            counts['cnv_checked']+=1
                            if c<pavcount:
                                counts['cnv_less_than_pav']+=1
                                vw.writerow(dict(kmer=key,accession=acc,pav_count=pavcount,
                                                 cnv_count=c,reason='cnv_less_than_pav'))
                                if pavcount>=1 and c==0:
                                    counts['cnv_zero_for_carrier']+=1
                            # PAV=0,CNV>0 is normal; never flag or filter.
                        weight_carriers = carriers if args.weight_scope=='all' else [x for x in carriers if x[0] in target_set]
                        copy_list=[cnv_values.get(a) for a,_,_ in weight_carriers]
                        if not weight_carriers:
                            counts['invalid_cnv']+=1
                            vw.writerow(dict(kmer=key,accession='',pav_count='',cnv_count='',reason='no_weight_scope_carriers'))
                            continue
                        if any(c is None or c<1 for c in copy_list):
                            counts['invalid_cnv']+=1
                            for (a,_,pavcount),c in zip(weight_carriers,copy_list):
                                if c is None or c<1:
                                    vw.writerow(dict(kmer=key,accession=a,pav_count=pavcount,
                                                     cnv_count='' if c is None else c,
                                                     reason='missing_or_zero_carrier_cnv'))
                            continue
                        m=len(weight_carriers); copy_sum=sum(copy_list); w=float(m)/float(copy_sum)
                    else:
                        m=n1; copy_sum=n1; w=1.0
                    raw1=np.array([y for _,y,_ in carriers],dtype=float)
                    raw0=np.array(absent,dtype=float)
                    mean1=float(raw1.mean()); mean0=float(raw0.mean())
                    var1=float(raw1.var(ddof=1)); var0=float(raw0.var(ddof=1))
                    adjmean=w*mean1+(1-w)*mean0
                    adjvar=w*var1+(1-w)*var0+w*(1-w)*(mean1-mean0)**2
                    result=stats.ttest_ind_from_stats(
                        mean1=adjmean,std1=math.sqrt(max(0.,adjvar)),nobs1=n1,
                        mean2=mean0,std2=math.sqrt(max(0.,var0)),nobs2=n0,equal_var=False)
                    t=float(result.statistic); p=float(result.pvalue)
                    if not math.isfinite(t) or not math.isfinite(p):
                        counts['invalid_test']+=1; continue
                    wr.writerow(dict(kmer=key,mean_present=adjmean,mean_absent=mean0,
                                     present_count=n1,absent_count=n0,t_stat=t,p_value=p,
                                     w=w,m=m,copy_sum=copy_sum,resistant_present=nr,other_present=no))
                    counts['tested']+=1
                    if key in debug:
                        print('[DEBUG] {} | n1={} n0={} | target_present={} other_present={} | '
                              'weight_scope={} CNV_carriers={} m={} copy_sum={} w={} | mean1={} mean0={} var1={} var0={} '
                              '| adjmean={} adjvar={} t={} p={}'.format(
                                  key,n1,n0,nr,no,args.weight_scope,copy_list if args.copy_adjust=='methods' else 'unused',
                                  m,copy_sum,w,mean1,mean0,var1,var0,adjmean,adjvar,t,p),
                              file=sys.stderr,flush=True)
            finally:
                if cf: cf.close()
        n_sig=compute_bh(raw,all_file,args.out,args.alpha,args.delimiter)
    print('\n========== Summary ==========',flush=True)
    for key,val in counts.items():
        print('{:<26}: {:,}'.format(key,val),flush=True)
    print('BH significant             : {:,}'.format(n_sig),flush=True)
    print('Significant {}            : {}'.format(args.output_format.upper(),args.out),flush=True)
    print('All tested {}             : {}'.format(args.output_format.upper(),all_file),flush=True)
    print('Validation {}             : {}'.format(args.output_format.upper(),validation_file),flush=True)
    print('Elapsed seconds            : {:.1f}'.format(time.time()-start),flush=True)


if __name__=='__main__':
    main()
