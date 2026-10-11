#!/usr/bin/env python3
"""
Extract the tarballs of scripts/laptop/collect.py into one directory, for
decide.py --units-dir DIR and scheduler.py ingest --counts-only DIR/U*.jsonl.
Python 3.8+, standard library only.

usage: unpack.py DIR results-*.tar.xz...

* Each tarball's files are checked against its manifest.json (sha256 of
  every unit and small file; no unit missing or extra). Members must be
  plain files with plain names (no directories, links or "..").
  A tarball that fails is refused, and then nothing is extracted.
* The tarballs are merged in the order they were written (the manifests'
  time and part number): a unit packed twice (run again with run.py
  --rerun-cut) keeps its later version, as do progress.log and units.tsv.
"""
import argparse
import hashlib
import json
import lzma
import os
import re
import shutil
import sys
import tarfile

SAFE = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]*$")
UNIT = re.compile(r"U\d+\.jsonl$")


class Refused(Exception):
    pass


def stage(path, sdir):
    """extract the tarball into the empty directory sdir and check it
    against its manifest; returns the manifest"""
    got = {}
    try:
        with tarfile.open(path, "r|xz") as tf:
            for m in tf:
                if not m.isreg():
                    raise Refused(f"member {m.name!r} is not a plain file")
                if not SAFE.match(m.name) or m.name in got:
                    raise Refused(f"member name {m.name!r} not allowed")
                h = hashlib.sha256()
                src = tf.extractfile(m)
                dst = os.path.join(sdir, m.name)
                with open(dst, "wb") as f:
                    while True:
                        b = src.read(1 << 20)
                        if not b:
                            break
                        h.update(b)
                        f.write(b)
                os.utime(dst, (m.mtime, m.mtime))
                got[m.name] = h.hexdigest()
    except (tarfile.TarError, lzma.LZMAError, EOFError, OSError) as e:
        raise Refused(f"cannot be read ({type(e).__name__}: {e})")
    if "manifest.json" not in got:
        raise Refused("no manifest.json (an older collect.sh's tarball? tar -xJf it)")
    try:
        with open(os.path.join(sdir, "manifest.json")) as f:
            man = json.load(f)
        want = dict(man["units"])
        files = dict(man["files"])
    except (ValueError, KeyError, TypeError) as e:
        raise Refused(f"manifest.json is not valid ({e})")
    want = {n: u["sha256"] for n, u in want.items()}
    want.update(files)
    del got["manifest.json"]
    if set(got) != set(want):
        miss = sorted(set(want) - set(got))[:5]
        extra = sorted(set(got) - set(want))[:5]
        raise Refused(f"files not as in its manifest (missing {miss}, not listed {extra})")
    bad = sorted(n for n in got if got[n] != want[n])
    if bad:
        raise Refused(f"sha256 does not match the manifest: {bad[:5]}")
    os.remove(os.path.join(sdir, "manifest.json"))
    return man


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("dir")
    ap.add_argument("tarballs", nargs="+")
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    tmp = os.path.join(a.dir, f".unpack-{os.getpid()}")
    os.makedirs(tmp)
    try:
        staged = []
        for k, path in enumerate(a.tarballs):
            sdir = os.path.join(tmp, str(k))
            os.makedirs(sdir)
            try:
                man = stage(path, sdir)
            except Refused as e:
                sys.exit(f"refused {path}: {e}")
            staged.append(((float(man.get("created", 0)), int(man.get("part", 0)),
                            os.path.basename(path)), sdir, path, man))
        staged.sort(key=lambda x: x[0])
        seen = {}
        for _, sdir, path, man in staged:
            for name in os.listdir(sdir):
                os.replace(os.path.join(sdir, name), os.path.join(a.dir, name))
                if UNIT.match(name):
                    seen[name] = seen.get(name, 0) + 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    twice = sum(1 for n in seen.values() if n > 1)
    total = sum(1 for n in os.listdir(a.dir) if UNIT.match(n))
    print(f"{len(staged)} tarballs checked: {len(seen)} units"
          + (f" ({twice} packed more than once: the later kept)" if twice else "")
          + f"; {a.dir} has {total} units")


if __name__ == "__main__":
    main()
