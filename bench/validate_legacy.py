#!/usr/bin/env python3
"""
Cross-check benchmark instances against the legacy arrangement program.

For every instance in the given file, write the enumeration file with
bin/enumerate, run the legacy bin/arrangement_{n} on it, and compute the
number of squares and the same order-independent hash that bench computes.
Prints the instance file with the legacy values filled in, and reports any
disagreement with expectations already present in the file.

usage: bench/validate_legacy.py bench/quick.txt [--bin bin] [--tmp DIR]
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time

MASK = (1 << 64) - 1


def fnv1a(values):
    h = 1469598103934665603
    for x in values:
        for b in range(8):
            h ^= (x >> (8 * b)) & 0xFF
            h = (h * 1099511628211) & MASK
    return h


def square_hash(rows, cols):
    r = sorted(tuple(sorted(v)) for v in rows)
    c = sorted(tuple(sorted(v)) for v in cols)
    fr = [x for v in r for x in v]
    fc = [x for v in c for x in v]
    first, second = (fr, fc) if fr <= fc else (fc, fr)
    return fnv1a(first + second)


def parse_legacy(path, n):
    """yield (rows, cols) for each solution in arrangement.c output"""
    with open(path) as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        if lines[i] == "solution found":
            rows = [tuple(map(int, lines[i + 1 + k].split()[1:])) for k in range(n)]
            cols = [tuple(map(int, lines[i + 2 + n + k].split()[1:])) for k in range(n)]
            yield rows, cols
            i += 2 + 2 * n
        else:
            i += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("instances")
    ap.add_argument("--bin", default="bin")
    ap.add_argument("--tmp", default=None)
    args = ap.parse_args()

    tmp = args.tmp or tempfile.mkdtemp(prefix="validate_legacy_")
    os.makedirs(tmp, exist_ok=True)
    bad = 0
    for line in open(args.instances):
        raw = line.rstrip("\n")
        if not raw.strip() or raw.strip().startswith("#"):
            print(raw)
            continue
        fields = [f.strip() for f in raw.split("|")]
        n, exps, S = int(fields[0]), fields[1].split(), int(fields[2])
        exp_count = fields[3] if len(fields) > 3 else "-"
        exp_hash = fields[4] if len(fields) > 4 else "-"
        efile = os.path.join(tmp, f"e_{n}_{'_'.join(exps)}_{S}.txt")
        afile = os.path.join(tmp, f"a_{n}_{'_'.join(exps)}_{S}.txt")
        subprocess.run([os.path.join(args.bin, "enumerate"), "--vec-size", str(n),
                        "--min-sum", str(S), "--max-sum", str(S), "--file", efile, *exps],
                       check=True, stderr=subprocess.DEVNULL)
        t = time.time()
        with open(afile, "w") as out:
            subprocess.run([os.path.join(args.bin, f"arrangement_{n}"), "--file", efile,
                            "--sum", str(S)], check=True, stdout=out)
        elapsed = time.time() - t
        count, h = 0, 0
        for rows, cols in parse_legacy(afile, n):
            count += 1
            h = (h + square_hash(rows, cols)) & MASK
        hs = f"{h:016x}"
        mismatch = (exp_count not in ("-", str(count))) or (exp_hash not in ("-", hs))
        bad += mismatch
        print(f"{n} | {' '.join(exps)} | {S} | {count} | {hs}"
              + (f"   # MISMATCH with {exp_count} {exp_hash}" if mismatch else ""))
        print(f"# legacy time {elapsed:.3f}s", file=sys.stderr)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
