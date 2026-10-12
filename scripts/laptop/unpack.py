#!/usr/bin/env python3
"""
Extract the tarballs of scripts/laptop/collect.py into one directory, for
decide.py --units-dir DIR and scheduler.py ingest --counts-only DIR/U*.jsonl.
Python 3.8+, standard library only.

usage: unpack.py DIR results-*.tar.xz...

* Each tarball's files are checked against its manifest.json (sha256 of
  every unit and small file; no unit missing or extra). Members must be
  plain files with plain names (no directories, links, sparse files or
  ".."), each under 1 GB and under 4 GB in all. A tarball that fails is
  refused, and then nothing is extracted.
* The tarballs are merged in the order they were written (the manifests'
  collect number "seq" and part number; an older collect.py's tarball,
  with a time in place of the number, before all the numbered ones): a
  unit packed twice (run again with run.py --rerun-cut) keeps its later
  version, as do progress.log and units.tsv. (The numbers count the
  collects of one output directory: tarballs of two machines are ordered
  by number, then by name.)
  DIR/unpacked.json records the version of each file in DIR, so this
  also holds across calls: an older tarball unpacked later into the same
  DIR does not replace a newer file.
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
RECORD = "unpacked.json"
MEMBER_MAX = 1 << 30    # bytes of one member
TOTAL_MAX = 4 << 30     # bytes of one tarball


class Refused(Exception):
    pass


def stage(path, sdir):
    """extract the tarball into the empty directory sdir and check it
    against its manifest; returns the manifest"""
    got = {}
    total = 0
    try:
        with tarfile.open(path, "r|xz") as tf:
            for m in tf:
                if not m.isreg() or m.type == tarfile.GNUTYPE_SPARSE or m.sparse is not None:
                    raise Refused(f"member {m.name!r} is not a plain file")
                if not SAFE.match(m.name) or m.name in got or m.name == RECORD:
                    raise Refused(f"member name {m.name!r} not allowed")
                if m.size > MEMBER_MAX or total + m.size > TOTAL_MAX:
                    raise Refused(f"member {m.name!r} too large ({m.size} bytes)")
                total += m.size
                h = hashlib.sha256()
                size = 0
                src = tf.extractfile(m)
                dst = os.path.join(sdir, m.name)
                with open(dst, "wb") as f:
                    while True:
                        b = src.read(1 << 20)
                        if not b:
                            break
                        size += len(b)
                        if size > m.size:
                            raise Refused(f"member {m.name!r} longer than its header")
                        h.update(b)
                        f.write(b)
                os.utime(dst, (m.mtime, m.mtime))
                got[m.name] = h.hexdigest()
    except (tarfile.TarError, lzma.LZMAError, EOFError, OSError, ValueError) as e:
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
            # (seq: collect.py's number; created: an older one's time)
            new = "seq" in man
            staged.append(((int(new), float(man.get("seq" if new else "created", 0)),
                            int(man.get("part", 0)), os.path.basename(path)), sdir, path, man))
        staged.sort(key=lambda x: x[0])
        # the version of each file already in DIR (an earlier call)
        rec_path = os.path.join(a.dir, RECORD)
        try:
            with open(rec_path) as f:
                rec = json.load(f)
        except FileNotFoundError:
            rec = {}
        except ValueError:
            sys.exit(f"{rec_path} is damaged: remove it (then the files given "
                     f"later replace those in {a.dir})")
        # (an older unpack.py's version [time, part, name]: before the
        # numbered ones)
        rec = {n: ([0] + v if len(v) == 3 else v) for n, v in rec.items()}
        seen = {}
        older = 0
        for key, sdir, path, man in staged:
            ver = list(key)
            for name in os.listdir(sdir):
                if UNIT.match(name):
                    seen[name] = seen.get(name, 0) + 1
                if name in rec and rec[name][:3] > ver[:3] \
                        and os.path.exists(os.path.join(a.dir, name)):
                    older += UNIT.match(name) is not None
                    continue
                os.replace(os.path.join(sdir, name), os.path.join(a.dir, name))
                rec[name] = ver
        with open(rec_path + ".tmp", "w") as f:
            json.dump(rec, f, indent=0, sort_keys=True)
        os.replace(rec_path + ".tmp", rec_path)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    twice = sum(1 for n in seen.values() if n > 1)
    total = sum(1 for n in os.listdir(a.dir) if UNIT.match(n))
    print(f"{len(staged)} tarballs checked: {len(seen)} units"
          + (f" ({twice} packed more than once: the later kept)" if twice else "")
          + (f"; {older} older than the copy already in {a.dir}: not replaced"
             if older else "")
          + f"; {a.dir} has {total} units")


if __name__ == "__main__":
    main()
