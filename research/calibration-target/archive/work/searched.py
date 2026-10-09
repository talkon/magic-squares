#!/usr/bin/env python3
"""(P, S) already searched (n = 6): every record with P and S in the
scratchpad's JSONL outputs (except calib-target), in the tar.xz archives of
JSONL there, and in research/wip/data/all_records.tsv.xz. Writes
searched.json: {"P_key": [S, ...]} plus per-P record kinds."""
import gzip, json, lzma, os, re, sys, tarfile, collections
SP = "/tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad"
OUT = os.path.join(SP, "calib-target/work/searched.json")
rP = re.compile(rb'"P":\[([0-9,]*)\]')
rS = re.compile(rb'"S":(\d+)')
rn = re.compile(rb'"n":(\d+)')
rt = re.compile(rb'"type":"(\w+)"')
seen = collections.defaultdict(set)
kinds = collections.Counter()

def feed(line):
    m = rP.search(line)
    if not m:
        return
    s = rS.search(line)
    if not s:
        return
    n = rn.search(line)
    if n and n.group(1) != b"6":
        return
    t = rt.search(line)
    t = t.group(1).decode() if t else "?"
    if t not in ("sum", "dsum", "csum", "dchunk", "square", "dsquare", "csquare"):
        return
    P = [int(x) for x in m.group(1).split(b",") if x]
    while P and P[-1] == 0:
        P.pop()
    seen["_".join(map(str, P))].add(int(s.group(1)))
    kinds[t] += 1

def opener(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    if path.endswith(".xz"):
        return lzma.open(path, "rb")
    return open(path, "rb")

nfiles = 0
for root, dirs, files in os.walk(SP):
    if os.path.abspath(root).startswith(os.path.join(SP, "calib-target")):
        dirs[:] = []
        continue
    for f in files:
        p = os.path.join(root, f)
        try:
            if re.search(r"\.jsonl(\.gz|\.xz)?$", f):
                with opener(p) as fh:
                    for line in fh:
                        feed(line)
                nfiles += 1
            elif f.endswith(".tar.xz") or f.endswith(".tar"):
                with tarfile.open(p) as tf:
                    for mem in tf:
                        if mem.isfile() and re.search(r"\.jsonl(\.gz)?$", mem.name):
                            fh = tf.extractfile(mem)
                            if mem.name.endswith(".gz"):
                                fh = gzip.open(fh)
                            for line in fh:
                                feed(line)
                            nfiles += 1
        except Exception as e:
            print("skip", p, e, file=sys.stderr)
with lzma.open("/home/user/magic-squares/research/wip/data/all_records.tsv.xz", "rb") as fh:
    for line in fh:
        feed(line)
json.dump({k: sorted(v) for k, v in seen.items()}, open(OUT, "w"))
print(nfiles, "files;", len(seen), "P;", sum(len(v) for v in seen.values()), "(P, S);", dict(kinds))
