#!/usr/bin/env python3
"""k-mer GWAS on a k-mer presence/absence VCF.

usage: kmer_vcf_gwas.py --vcf kmers.vcf.gz --pheno traits.tsv --out results --threads 56
"""
import argparse
import gzip
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import get_context

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
import numpy as np  # noqa: E402
from scipy import special  # noqa: E402

KEEP_LOGP = 3.0         # approximate -log10 P needed for the GEMMA pass
KINSHIP_MAF = 0.05      # MAF and minor-allele count of the k-mers used for the kinship matrix,
KINSHIP_MAC = 5         # whatever --maf and --mac are
MAX_RANGES = 2000       # byte ranges the VCF is cut into
MIN_RANGE = 1 << 20     # smallest byte range, so small VCFs get fewer ranges
KINSHIP_RANGES = 61     # evenly spaced byte ranges used for the kinship matrix
TEXT_LINES = 100_000    # records parsed at a time
MISSING = {"", "NA", "na", "NaN", "nan", ".", "-"}
GEMMA_COLUMNS = "chr\trs\tps\tn_miss\tallele1\tallele0\taf\tlogl_H1\tl_mle\tp_lrt\n"
# F histogram for lambda: 0.001 steps up to 4, 0.1 steps up to 60
EDGES = np.concatenate([np.arange(0, 4, 0.001), np.arange(4, 60.0001, 0.1)])
NB = len(EDGES) + 1
BGZF_MAGIC = b"\x1f\x8b\x08\x04"
BGZF_XTRA = b"\x06\x00\x42\x43\x02\x00"
BGZF_HDR = 18
_G = {}
_GEMMA_PROCS = []


def log(msg):
    print(f"[{time.strftime('%F %T')}] {msg}", flush=True)


def die(msg):
    for p in _GEMMA_PROCS:
        if p.poll() is None:
            p.kill()
    sys.exit(f"ERROR: {msg}")


def write_json(path, obj):
    with open(path + ".tmp", "w") as fh:
        json.dump(obj, fh, indent=1)
    os.replace(path + ".tmp", path)


def min_count(n, maf, mac):
    return max(mac, int(np.ceil(maf * n)))


def maf_ok(count, n, maf, mac):
    lo = min_count(n, maf, mac)
    return (count >= lo) & (count <= n - lo)


def read_vcf_layout(path):
    with open(path, "rb") as fh:
        if not block_size_at(fh.read(BGZF_HDR), 0):
            die(f"{path} is not compressed with bgzip")
    samples = first = None
    with gzip.open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b"##"):
                continue
            if line.startswith(b"#CHROM"):
                samples = line.rstrip(b"\r\n").decode().split("\t")[9:]
                continue
            first = line
            break
    if not samples or first is None:
        die(f"{path}: no #CHROM line or no records")
    if len(set(samples)) != len(samples):
        die(f"{path}: duplicated sample names")
    f = first.rstrip(b"\n").split(b"\t")
    if len(f) != 9 + len(samples) or f[8] != b"GT":
        die(f"{path}: records must have FORMAT GT and one genotype per sample")
    prefix = sum(len(c) + 1 for c in f[:9])
    lay = dict(prefix=prefix, id_start=len(f[0]) + len(f[1]) + 2, k=len(f[2]), width=prefix + 4 * len(samples))
    if len(first) != lay["width"]:
        die(f"{path}: genotypes must be 0/0 or 1/1")
    return samples, lay


def read_pheno(path):
    with open(path) as fh:
        lines = [ln.rstrip("\r\n") for ln in fh if ln.strip()]
    if len(lines) < 2:
        die(f"{path}: needs a header line and data lines")
    sep = "\t" if "\t" in lines[0] else None
    head = [h.strip() for h in lines[0].split(sep)]
    traits = head[1:]
    if not traits:
        die(f"{path}: no trait columns")
    for t in traits:
        if not t or t.startswith(".") or "/" in t or t in ("work", "summary.tsv"):
            die(f"{path}: trait name {t!r} cannot be used as an output directory name")
    if len(set(traits)) != len(traits):
        die(f"{path}: duplicated trait names")
    values = {t: {} for t in traits}
    seen = set()
    for ln in lines[1:]:
        f = [x.strip() for x in ln.split(sep)]
        if len(f) != len(head):
            die(f"{path}: {len(f)} fields instead of {len(head)} in line '{ln[:60]}'")
        if f[0] in seen:
            die(f"{path}: sample {f[0]} appears twice")
        seen.add(f[0])
        for t, v in zip(traits, f[1:]):
            if v in MISSING:
                continue
            try:
                x = float(v)
            except ValueError:
                die(f"{path}: value {v!r} of {f[0]} ({t}) is not a number; write missing values as NA")
            if x == -9 or not np.isfinite(x):
                die(f"{path}: value {v!r} of {f[0]} ({t}); write missing values as NA")
            values[t][f[0]] = x
    return traits, values


def make_runs(samples, traits, values, base):
    pos = {s: i for i, s in enumerate(samples)}
    runs = []
    for t in traits:
        ph = values[t]
        idx = sorted(pos[s] for s in ph if s in pos)
        if len(idx) < 20:
            die(f"trait {t}: only {len(idx)} phenotyped samples are in the VCF")
        if 2 * min_count(len(idx), base["maf"], base["mac"]) > len(idx):
            die(f"trait {t}: with {len(idx)} samples no k-mer can pass --maf {base['maf']} and --mac {base['mac']}")
        y = [ph[samples[i]] for i in idx]
        key = json.dumps([base["vcf_size"], base["maf"], base["mac"], [samples[i] for i in idx], y])
        runs.append(dict(trait=t, idx=idx, y=y, key=hashlib.sha256(key.encode()).hexdigest()[:16],
                         not_in_vcf=sorted(s for s in ph if s not in pos)))
    return runs


def prepare_work(work, base, todo, vcf):
    path = os.path.join(work, "config.json")
    if os.path.exists(path):
        cfg = json.load(open(path))
        keys = {r["trait"]: r["key"] for r in cfg["runs"]}
        if any(cfg.get(k) != v for k, v in base.items()) or any(keys.get(r["trait"]) != r["key"] for r in todo):
            die(f"{work} was made from other inputs; remove it to start again")
        log(f"continuing in {work}")
    else:
        sets, runs = {}, []
        for r in todo:
            sid = sets.setdefault(tuple(r["idx"]), len(sets))
            runs.append(dict(trait=r["trait"], key=r["key"], set=sid, y=r["y"]))
        cfg = dict(base, sets=[dict(id=i, idx=list(ix)) for ix, i in sets.items()], runs=runs)
        os.makedirs(work, exist_ok=True)
        write_json(path, cfg)
    cfg["vcf"] = vcf
    return cfg


def block_size_at(buf, i):
    if i + BGZF_HDR > len(buf) or buf[i:i + 4] != BGZF_MAGIC or buf[i + 10:i + 16] != BGZF_XTRA:
        return 0
    size = struct.unpack("<H", buf[i + 16:i + 18])[0] + 1
    return size if size > BGZF_HDR else 0


def find_block_start(fh, offset, probe=1 << 20, chain=3):
    pos = offset
    fsize = os.fstat(fh.fileno()).st_size
    while pos < fsize:
        fh.seek(pos)
        buf = fh.read(probe)
        if not buf:
            return -1
        i = 0
        while True:
            i = buf.find(BGZF_MAGIC, i)
            if i < 0 or i + BGZF_HDR > len(buf):
                break
            ok, j, n = True, i, 0
            while n < chain:
                bs = block_size_at(buf, j)
                if bs == 0:
                    ok = j + BGZF_HDR > len(buf)
                    break
                j += bs
                n += 1
            if ok and block_size_at(buf, i):
                return pos + i
            i += 1
        pos += max(probe - BGZF_HDR, 1)
    return -1


def iter_blocks(fh, pos, read_size):
    fh.seek(pos)
    buf, base, i = b"", pos, 0
    while True:
        if len(buf) - i < BGZF_HDR:
            more = fh.read(read_size)
            if not more:
                return
            buf, base, i = buf[i:] + more, base + i, 0
            if len(buf) < BGZF_HDR:
                return
        bsize = block_size_at(buf, i)
        if bsize == 0:
            raise RuntimeError(f"no BGZF block header at byte {base + i}")
        while len(buf) - i < bsize:
            more = fh.read(read_size)
            if not more:
                raise RuntimeError(f"truncated BGZF block at byte {base + i}")
            buf, base, i = buf[i:] + more, base + i, 0
        yield base + i, zlib.decompress(buf[i + BGZF_HDR:i + bsize - 8], -15)
        i += bsize


def byte_ranges(cfg):
    n = cfg["n_ranges"]
    e = np.linspace(0, cfg["vcf_size"], n + 1).astype(np.int64)
    return [(i, int(e[i]), int(e[i + 1]), i == 0, i == n - 1) for i in range(n)]


def stream_records(vcf, task, lay, nsamp, read_size=1 << 24):
    # a range starts after its first newline and ends at the first newline past its end
    rid, start, end, is_first, is_last = task
    width, prefix = lay["width"], lay["prefix"]
    k0, k1 = lay["id_start"], lay["id_start"] + lay["k"]
    gt_cols = prefix + 4 * np.arange(nsamp)

    def parse(text):
        n = len(text) // width
        if n == 0:
            return 0, None
        arr = np.frombuffer(bytes(text[:n * width]), dtype=np.uint8).reshape(n, width)
        if not ((arr[:, width - 1] == 10).all() and (arr[:, prefix - 1] == 9).all()
                and (arr[:, k0 - 1] == 9).all() and (arr[:, k1] == 9).all()):
            raise RuntimeError("VCF records differ in length")
        gt = arr[:, gt_cols]
        if (((gt != 48) & (gt != 49)) | (arr[:, gt_cols + 2] != gt)).any():
            raise RuntimeError("genotype other than 0/0 or 1/1")
        return n * width, (arr[:, k0:k1].copy(), (gt == 48).astype(np.uint8))

    with open(vcf, "rb") as fh:
        pos = 0 if is_first else find_block_start(fh, start)
        if pos < 0 or (not is_last and pos >= end):
            return
        text = bytearray()
        head_done, tail_from = False, None
        for bstart, data in iter_blocks(fh, pos, read_size):
            if not is_last and bstart >= end and tail_from is None:
                tail_from = len(text)
            text += data
            if not head_done:
                if is_first:
                    k = text.find(b"#CHROM")
                    nl = text.find(b"\n", k) if k >= 0 else -1
                else:
                    nl = text.find(b"\n")
                if nl < 0:
                    continue
                del text[:nl + 1]
                if tail_from is not None:
                    tail_from = max(tail_from - (nl + 1), 0)
                head_done = True
                if not is_first and text[:1] == b"#":
                    raise RuntimeError("a byte range starts inside the VCF header")
            if tail_from is not None:
                nl = text.find(b"\n", tail_from)
                if nl >= 0:
                    del text[nl + 1:]
                    break
                continue
            if len(text) >= TEXT_LINES * width:
                used, blk = parse(text)
                del text[:used]
                if blk is not None:
                    yield blk
        used, blk = parse(text)
        del text[:used]
        if blk is not None:
            yield blk
        if text:
            raise RuntimeError(f"byte range {rid}: {len(text)} trailing bytes")


def _init(d):
    _G.clear()
    _G.update(d)


def kin_file(work, sid, ext):
    return os.path.join(work, "kinship", f"set{sid:03d}{ext}")


def kinship_worker(args):
    vcf, task, lay, nsamp = args
    acc = [[np.zeros((len(ix), len(ix))), 0] for ix in _G["sets"]]
    for _, pres in stream_records(vcf, task, lay, nsamp):
        for j, ix in enumerate(_G["sets"]):
            x = pres[:, ix]
            x = x[maf_ok(x.sum(1), len(ix), KINSHIP_MAF, KINSHIP_MAC)].astype(np.float64)
            if not len(x):
                continue
            x -= x.mean(1, keepdims=True)
            acc[j][0] += x.T @ x
            acc[j][1] += len(x)
    return acc


def run_kinship(cfg, work, threads):
    """centered relatedness matrix K = Xc Xc' / m (GEMMA -gk 1) per sample set"""
    todo = [s for s in cfg["sets"] if not os.path.exists(kin_file(work, s["id"], ".npz"))]
    if not todo:
        return
    os.makedirs(os.path.join(work, "kinship"), exist_ok=True)
    ranges = byte_ranges(cfg)[::cfg["kinship_stride"]]
    log(f"kinship: {len(todo)} sample set(s)")
    total = [[np.zeros((len(s["idx"]), len(s["idx"]))), 0] for s in todo]
    jobs = [(cfg["vcf"], t, cfg["layout"], len(cfg["samples"])) for t in ranges]
    with get_context("fork").Pool(threads, initializer=_init,
                                  initargs=({"sets": [np.array(s["idx"]) for s in todo]},)) as pool:
        for acc in pool.imap(kinship_worker, jobs):
            for j, (m, n) in enumerate(acc):
                total[j][0] += m
                total[j][1] += n
    for s, (m, n) in zip(todo, total):
        if n == 0:
            die(f"no k-mer with MAF >= {KINSHIP_MAF} and MAC >= {KINSHIP_MAC} for the kinship matrix "
                f"of sample set {s['id']}")
        kin = m / n
        p = kin_file(work, s["id"], "")
        np.savetxt(p + ".gemma.txt.tmp", kin, fmt="%.8g", delimiter="\t")
        os.replace(p + ".gemma.txt.tmp", p + ".gemma.txt")
        with open(p + ".npz.tmp", "wb") as fh:
            np.savez(fh, K=kin, n_kmers=n)
        os.replace(p + ".npz.tmp", p + ".npz")


def reml_delta(y, S, U):
    """EMMA REML on the spectral decomposition of K; returns delta = sigma_e^2 / sigma_g^2"""
    n = len(y)
    yr, xr = U.T @ y, U.T @ np.ones(n)

    def nll(logd):
        w = 1.0 / (S + np.exp(logd))
        xwx = (w * xr * xr).sum()
        xwy = (w * xr * yr).sum()
        rss = (w * yr * yr).sum() - xwy * xwy / xwx
        return 0.5 * ((n - 1) * np.log(rss / (n - 1)) + np.log(1 / w).sum() + np.log(xwx))

    grid = np.linspace(-10, 10, 201)
    vals = np.array([nll(g) for g in grid])
    i = int(np.argmin(vals))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    for _ in range(60):  # golden-section refinement
        m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
        if nll(m1) < nll(m2):
            hi = m2
        else:
            lo = m1
    return float(np.exp((lo + hi) / 2))


def build_sets(cfg, work):
    # w = 1/(S + delta), a = U'1, b = U'x; F = (Sxy^2 / x'x) / ((y'y - Sxy^2 / x'x) / (n - 2)) after whitening
    nsamp = len(cfg["samples"])
    sets = []
    for s in cfg["sets"]:
        runs = [j for j, r in enumerate(cfg["runs"]) if r["set"] == s["id"]]
        S, U = np.linalg.eigh(np.load(kin_file(work, s["id"], ".npz"))["K"])
        S = np.clip(S, 0, None)
        n = len(s["idx"])
        a_ = U.T @ np.ones(n)
        V, Wa, W, Saa, yty = [], [], [], [], []
        for j in runs:
            y = np.asarray(cfg["runs"][j]["y"], dtype=np.float64)
            w = 1.0 / (S + reml_delta(y, S, U))
            sw = np.sqrt(w)
            astar, ystar = sw * a_, sw * (U.T @ y)
            ytil = ystar - astar * (astar @ ystar) / (astar @ astar)
            V.append(sw * ytil)
            Wa.append(w * a_)
            W.append(w)
            Saa.append(astar @ astar)
            yty.append(ytil @ ytil)
        # U padded to all VCF samples (zero rows outside the set), so presence rows @ P = x_set @ U
        Pm = np.zeros((nsamp, n), dtype=np.float32)
        Pm[np.array(s["idx"])] = U
        df = n - 2
        sets.append(dict(idx=np.array(s["idx"]), n=n, df=df, P=Pm, runs=runs,
                         VWa=np.ascontiguousarray(np.array(V + Wa).T, dtype=np.float32),
                         W=np.ascontiguousarray(np.array(W).T, dtype=np.float32),
                         Saa=np.array(Saa), yty=np.array(yty),
                         fthr=special.stdtrit(df, 0.5 * 10.0 ** (-KEEP_LOGP)) ** 2 * (1 - 1e-6)))
    return sets


def f_to_logp(f, df):
    p = 2 * special.stdtr(df, -np.sqrt(np.maximum(f, 0)))
    with np.errstate(divide="ignore"):
        return -np.log10(np.maximum(p, 1e-300))


def f_bin(f):
    lo = np.floor(f * 1000.0) + 1
    hi = 4000 + np.floor((f - 4.0) * 10.0) + 1
    return np.minimum(np.where(f < 4.0, lo, hi), NB - 1).astype(np.int64)


def scan_init(d):
    _init(d)
    setmat = np.zeros((d["nsamp"], len(d["sets"])), dtype=np.float32)
    for j, s in enumerate(d["sets"]):
        setmat[s["idx"], j] = 1
    _G["setmat"] = setmat


def range_prefix(work, i):
    return os.path.join(work, "scan", f"range{i:05d}")


def scan_worker(args):
    vcf, task, lay, nsamp, maf, mac, pre = args
    sets, nrun = _G["sets"], _G["nrun"]
    hist = np.zeros((nrun, NB), dtype=np.int64)
    tested = np.zeros(nrun, dtype=np.int64)
    kept = [([], []) for _ in range(nrun)]
    u_kmers, u_bits, n_union, records = [], [], 0, 0
    for kmers, pres in stream_records(vcf, task, lay, nsamp):
        nb = len(kmers)
        records += nb
        hits, any_keep = [], np.zeros(nb, dtype=bool)
        pf = pres.astype(np.float32)
        counts = pf @ _G["setmat"]
        for js, s in enumerate(sets):
            rows = np.flatnonzero(maf_ok(counts[:, js], s["n"], maf, mac))
            if not len(rows):
                continue
            b = pf[rows] @ s["P"]
            T = len(s["runs"])
            P = (b @ s["VWa"]).astype(np.float64)
            Sxx = ((b * b) @ s["W"]).astype(np.float64)
            del b
            Sxy, Sxa = P[:, :T], P[:, T:]
            xx = np.maximum(Sxx - Sxa * Sxa / s["Saa"], 1e-12)
            ssr = Sxy * Sxy / xx
            F = ssr / (np.maximum(s["yty"] - ssr, 1e-300) / s["df"])
            for c, j in enumerate(s["runs"]):
                tested[j] += len(rows)
                hist[j] += np.bincount(f_bin(F[:, c]), minlength=NB)
                sel = np.flatnonzero(F[:, c] >= s["fthr"])
                if not len(sel):
                    continue
                lp = f_to_logp(F[sel, c], s["df"])
                ok = lp >= KEEP_LOGP
                if not ok.any():
                    continue
                sel = sel[ok]
                hits.append((j, rows[sel], lp[ok].astype(np.float32)))
                any_keep[rows[sel]] = True
        if any_keep.any():
            u = np.flatnonzero(any_keep)
            where = np.full(nb, -1, dtype=np.int64)
            where[u] = np.arange(len(u)) + n_union
            u_kmers.append(kmers[u])
            u_bits.append(np.packbits(pres[u], axis=1))
            for j, rec, lp in hits:
                kept[j][0].append(where[rec].astype(np.uint32))
                kept[j][1].append(lp)
            n_union += len(u)
    res = {"hist": hist.astype(np.int32), "tested": tested}
    for j, (ix, lp) in enumerate(kept):
        if ix:
            res[f"r{j}_idx"] = np.concatenate(ix)
            res[f"r{j}_logp"] = np.concatenate(lp)
    np.save(pre + ".kmers.tmp.npy", np.concatenate(u_kmers) if u_kmers else np.zeros((0, lay["k"]), np.uint8))
    np.save(pre + ".bits.tmp.npy", np.concatenate(u_bits) if u_bits else np.zeros((0, (nsamp + 7) // 8), np.uint8))
    np.savez(pre + ".tmp.npz", **res)
    for ext in (".kmers", ".bits"):
        os.replace(pre + ext + ".tmp.npy", pre + ext + ".npy")
    os.replace(pre + ".tmp.npz", pre + ".npz")
    write_json(pre + ".json", dict(records=records))
    return records


def run_scan(cfg, work, threads):
    os.makedirs(os.path.join(work, "scan"), exist_ok=True)
    todo = [t for t in byte_ranges(cfg) if not os.path.exists(range_prefix(work, t[0]) + ".json")]
    if not todo:
        return
    sets = build_sets(cfg, work)
    nsamp = len(cfg["samples"])
    log(f"scan: {len(cfg['runs'])} trait(s)")
    jobs = [(cfg["vcf"], t, cfg["layout"], nsamp, cfg["maf"], cfg["mac"], range_prefix(work, t[0])) for t in todo]
    with get_context("fork").Pool(threads, initializer=scan_init,
                                  initargs=(dict(sets=sets, nsamp=nsamp, nrun=len(cfg["runs"])),)) as pool:
        for _ in pool.imap_unordered(scan_worker, jobs):
            pass


def load_scan(cfg, work, threads):
    def one(i):
        pre = range_prefix(work, i)
        z = np.load(pre + ".npz")
        return dict(arr={k: z[k] for k in z.files}, kmers=np.load(pre + ".kmers.npy"), bits=np.load(pre + ".bits.npy"))

    with ThreadPoolExecutor(threads) as ex:
        return list(ex.map(one, range(cfg["n_ranges"])))


def gather_run(j, scan, lay):
    parts = [d for d in scan if f"r{j}_idx" in d["arr"]]
    if not parts:
        return np.zeros((0, lay["k"]), np.uint8), np.zeros((0, 1), np.uint8)
    lp = np.concatenate([d["arr"][f"r{j}_logp"] for d in parts])
    km = np.concatenate([d["kmers"][d["arr"][f"r{j}_idx"]] for d in parts])
    bits = np.concatenate([d["bits"][d["arr"][f"r{j}_idx"]] for d in parts])
    order = np.argsort(-lp, kind="stable")
    return km[order], bits[order]


def write_bed(prefix, pres, ids):
    """presence -> homozygous allele1 'P' (code 00), absence -> homozygous allele2 'A' (code 11)"""
    nsamp = pres.shape[1]
    codes = np.where(pres == 1, 0, 3).astype(np.uint8)
    pad = (-nsamp) % 4
    if pad:
        codes = np.concatenate([codes, np.zeros((len(codes), pad), np.uint8)], axis=1)
    c = codes.reshape(len(codes), -1, 4)
    packed = (c[:, :, 0] | (c[:, :, 1] << 2) | (c[:, :, 2] << 4) | (c[:, :, 3] << 6)).astype(np.uint8)
    with open(prefix + ".bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        fh.write(packed.tobytes())
    with open(prefix + ".bim", "w") as fh:
        for i in ids:
            fh.write(f"0\t{i}\t0\t0\tP\tA\n")


def start_gemma(cfg, work, j, km, bits, jobs, gemma):
    r = cfg["runs"][j]
    wd = os.path.join(work, "gemma", f"r{j}")
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    idx = np.array(cfg["sets"][r["set"]]["idx"])
    names = [cfg["samples"][i] for i in idx]
    n = len(idx)
    # GEMMA's own MAF filter must not drop k-mers that passed the scan
    gmaf = f"{min(0.01, (min_count(n, cfg['maf'], cfg['mac']) - 0.5) / n):.6g}"
    procs = []
    for c, ix in enumerate(x for x in np.array_split(np.arange(len(km)), jobs) if len(x)):
        pre = os.path.join(wd, f"c{c:03d}")
        pres = np.unpackbits(bits[ix], axis=1, count=len(cfg["samples"]))[:, idx]
        write_bed(pre, pres, [f"{bytes(k).decode()}_{i + 1}" for i, k in zip(ix, km[ix])])
        with open(pre + ".fam", "w") as fh:
            for sm, v in zip(names, r["y"]):
                fh.write(f"{sm}\t{sm}\t0\t0\t0\t{float(v):.10g}\n")
        cmd = [gemma, "-bfile", pre, "-k", kin_file(work, r["set"], ".gemma.txt"), "-lmm", "2",
               "-maf", gmaf, "-miss", "0.5", "-outdir", wd, "-o", f"c{c:03d}"]
        with open(pre + ".log", "w") as lf:
            procs.append((subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=wd), pre))
        _GEMMA_PROCS.append(procs[-1][0])
    return procs, wd


def lambda_from_hist(hist, df):
    tot = hist.sum()
    cum = np.cumsum(hist)
    b = int(np.searchsorted(cum, tot / 2))
    lo = 0.0 if b == 0 else EDGES[b - 1]
    hi = EDGES[b] if b < len(EDGES) else EDGES[-1]
    before = cum[b - 1] if b > 0 else 0
    med = lo + (hi - lo) * (tot / 2 - before) / max(hist[b], 1)
    return float(special.chdtri(1, special.fdtrc(1, df, med)) / special.chdtri(1, 0.5))


def count_assoc(files, tested):
    if not files:
        return dict(kmers_in_assoc=0, n_bonferroni=0, max_logp=None)
    prog = ('FNR > 1 { p = $10 + 0; n++; if (p <= pb) c++; if (n == 1 || p < m) m = p } '
            'END { printf "%d %d %.10e\\n", n, c, m }')
    out = subprocess.run(["awk", "-F", "\t", "-v", f"pb={0.05 / tested!r}", prog, *files],
                         check=True, capture_output=True, text=True).stdout.split()
    m = float(out[2])
    return dict(kmers_in_assoc=int(out[0]), n_bonferroni=int(out[1]),
                max_logp=round(float(-np.log10(m)), 4) if m > 0 else None)


def merge_gz(parts, dest, threads):
    pigz = shutil.which("pigz")
    cmd = [pigz, "-p", str(threads), "-c"] if pigz else ["gzip", "-c"]
    with open(dest + ".tmp", "wb") as fo:
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=fo)
        if not parts:
            p.stdin.write(GEMMA_COLUMNS.encode())
        for i, f in enumerate(parts):
            with open(f, "rb") as fi:
                head = fi.readline()
                if i == 0:
                    p.stdin.write(head)
                shutil.copyfileobj(fi, p.stdin, 1 << 24)
        p.stdin.close()
        if p.wait() != 0:
            die(f"compression of {dest} failed")
    os.replace(dest + ".tmp", dest)


def finish_run(cfg, out, j, procs, wd, tested, hist, pz):
    r = cfg["runs"][j]
    for p, pre in procs:
        if p.wait() != 0:
            die(f"GEMMA failed for {r['trait']}, see {pre}.log")
    parts = [pre for _, pre in procs]
    for pre in parts:
        os.remove(pre + ".bed")
    assoc = [pre + ".assoc.txt" for pre in parts]
    n = len(r["y"])
    cnt = count_assoc(assoc, tested[j])
    pve = None
    if parts and os.path.exists(parts[0] + ".log.txt"):
        for line in open(parts[0] + ".log.txt"):
            if "pve estimate in the null model" in line:
                pve = float(line.split("=")[-1])
    d = os.path.join(out, r["trait"])
    os.makedirs(d, exist_ok=True)
    merge_gz(assoc, os.path.join(d, "assoc.txt.gz"), pz)
    shutil.rmtree(wd, ignore_errors=True)
    log(f"{r['trait']}: n {n}, {int(tested[j]):,} k-mers tested, Bonferroni -log10 P "
        f"{-np.log10(0.05 / tested[j]):.4f}, {cnt['n_bonferroni']:,} above, max {cnt['max_logp']}, "
        f"lambda {lambda_from_hist(hist[j], n - 2):.4f}, null pve {pve}")


def run_exact(cfg, work, out, gemma, threads, todo):
    if any(not os.path.exists(range_prefix(work, i) + ".json") for i in range(cfg["n_ranges"])):
        die("scan is incomplete")
    runs = [j for j, r in enumerate(cfg["runs"]) if r["trait"] in todo]
    log(f"gemma: {len(runs)} trait(s)")
    scan = load_scan(cfg, work, threads)
    tested = {j: sum(int(d["arr"]["tested"][j]) for d in scan) for j in runs}
    hist = {j: sum(d["arr"]["hist"][j].astype(np.int64) for d in scan) for j in runs}
    pz = min(4, threads)
    jobs = max(1, threads - pz - 2)
    pending = None
    for j in runs:
        trait = cfg["runs"][j]["trait"]
        if tested[j] == 0:
            die(f"{trait}: no k-mer passed --maf and --mac")
        km, bits = gather_run(j, scan, cfg["layout"])
        if pending is not None:
            for p, _ in pending[1]:
                p.wait()
        procs, wd = start_gemma(cfg, work, j, km, bits, jobs, gemma)
        del km, bits
        if pending is not None:
            pj, pprocs, pwd = pending
            finish_run(cfg, out, pj, pprocs, pwd, tested, hist, pz)
        pending = (j, procs, wd)
    pj, pprocs, pwd = pending
    finish_run(cfg, out, pj, pprocs, pwd, tested, hist, pz)




def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="k-mer presence/absence VCF compressed with bgzip")
    ap.add_argument("--pheno", required=True, help="phenotype table: sample ID column, then one column per trait")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--threads", type=int, default=8, help="parallel processes (default 8)")
    ap.add_argument("--maf", type=float, default=0.05,
                    help="minimum minor-allele frequency of tested k-mers (default 0.05; 0 = --mac only)")
    ap.add_argument("--mac", type=int, default=5, help="minimum minor-allele count of tested k-mers (default 5)")
    ap.add_argument("--gemma", default="gemma", help="GEMMA executable (default: gemma on PATH)")
    a = ap.parse_args()
    if not 0 <= a.maf < 0.5:
        die("--maf must be at least 0 and below 0.5")
    if a.mac < 1:
        die("--mac must be at least 1")
    if a.threads < 1:
        die("--threads must be at least 1")
    gemma = shutil.which(a.gemma)
    if not gemma:
        die(f"GEMMA not found ({a.gemma}); give its path with --gemma")

    vcf = os.path.abspath(a.vcf)
    samples, lay = read_vcf_layout(vcf)
    traits, values = read_pheno(a.pheno)
    size = os.path.getsize(vcf)
    n_ranges = int(min(MAX_RANGES, max(1, size // MIN_RANGE)))
    base = dict(vcf_size=size, samples=samples, layout=lay, maf=a.maf, mac=a.mac, n_ranges=n_ranges,
                kinship_stride=-(-n_ranges // KINSHIP_RANGES))
    out = os.path.abspath(a.out)
    work = os.path.join(out, "work")
    todo = {}
    for r in make_runs(samples, traits, values, base):
        done = os.path.join(out, r["trait"], "assoc.txt.gz")
        if os.path.exists(done):
            log(f"{r['trait']}: assoc.txt.gz exists, skipped")
            continue
        if r["not_in_vcf"]:
            log(f"{r['trait']}: {len(r['not_in_vcf'])} phenotyped sample(s) not in the VCF: "
                f"{', '.join(r['not_in_vcf'][:5])}{' ...' if len(r['not_in_vcf']) > 5 else ''}")
        todo[r["trait"]] = r

    if todo:
        cfg = prepare_work(work, base, list(todo.values()), vcf)
        run_kinship(cfg, work, a.threads)
        run_scan(cfg, work, a.threads)
        run_exact(cfg, work, out, gemma, a.threads, set(todo))
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
