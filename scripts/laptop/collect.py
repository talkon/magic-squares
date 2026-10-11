#!/usr/bin/env python3
"""
Pack the outputs of scripts/laptop/run.py that were not packed before.
Python 3.8+, standard library only.

usage: collect.py [DIR]     (DIR: run.py's --out, default laptop-runs)

* Writes research/laptop/results-<host>-<YYYYmmdd-HHMM>[-partNN].tar.xz
  with the unit files (U<i>.jsonl) not packed before, and those changed
  since they were packed (run.py --rerun-cut replaces a unit). The units
  are split over parts of at most ~45 MB each (GitHub warns at 50 MB and
  refuses files over 100 MB); no part is over 49 MB.
* Every part also carries the small files: meta_*.json, progress.log,
  units.tsv, MAGIC.txt and the build's record (build-laptop/build_info.txt
  and check_quick.txt, as build-laptop_*.txt). And a manifest.json: its
  units and files with their sha256 (scripts/laptop/unpack.py checks them).
* DIR/collected.json records each unit packed (sha256, size, tarball). It
  is written after each part is complete, so an interrupted collect is
  simply run again. With nothing new it writes no tarball.
"""
import argparse
import hashlib
import io
import json
import lzma
import os
import re
import socket
import subprocess
import sys
import tarfile
import time

UNIT = re.compile(r"U\d+\.jsonl$")
SMALL = re.compile(r"(meta_.*\.json|progress\.log|units\.tsv|MAGIC\.txt)$")
BUILD = ("build_info.txt", "check_quick.txt")
RECORD = "collected.json"
# MB: the target size of a part, and the hard limit (a part found over it
# is written again with fewer units)
PART_MB = 45.0
HARD_MB = 49.0
# raw bytes the compressor may hold before its output shows (measured
# under 0.6 MB at presets 6 to 9e; JSONL compresses ~7-8x, counted as 2x)
LAG = 4 << 20
# xz preset: the default, 6. On 255 MB of unit files: 7.5x in 125 s with
# 0.1 GB of memory; 6 | EXTREME 3.5% smaller in 1.5x the time; 9 | EXTREME
# 5.7% smaller in 2.2x the time with 0.7 GB
PRESET = 6


class XzSink:
    """the file-like object tarfile writes to: compresses into f and counts
    the bytes in and out"""

    def __init__(self, f, preset):
        self.f = f
        self.c = lzma.LZMACompressor(format=lzma.FORMAT_XZ, check=lzma.CHECK_CRC64,
                                     preset=preset)
        self.raw = 0
        self.out = 0

    def write(self, b):
        if self.f.closed:   # (tarfile's cleanup after an interrupt)
            return len(b)
        self.raw += len(b)
        x = self.c.compress(b)
        self.out += len(x)
        self.f.write(x)
        return len(b)

    def finish(self):
        x = self.c.flush()
        self.out += len(x)
        self.f.write(x)


def add(tf, name, data, mtime):
    ti = tarfile.TarInfo(name)
    ti.size = len(data)
    ti.mtime = int(mtime)
    ti.mode = 0o644
    tf.addfile(ti, io.BytesIO(data))


def read(path):
    """(bytes, size, mtime_ns) of a file, the stat of the file read"""
    with open(path, "rb") as f:
        st = os.fstat(f.fileno())
        data = f.read()
    return data, st.st_size, st.st_mtime_ns


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def write_part(tmp, small, names, out_dir, limit, cap, info):
    """write the part to tmp: the small files, then units from names while
    the part stays under limit bytes (at least one, at most cap), then
    manifest.json. Returns ([(name, sha256, size, mtime_ns)], file size)"""
    packed = []
    with open(tmp, "wb") as f:
        sink = XzSink(f, info["preset"])
        tf = tarfile.open(fileobj=sink, mode="w|", format=tarfile.PAX_FORMAT)
        for name, (data, mtime) in small.items():
            add(tf, name, data, mtime)
        for name in names:
            if cap is not None and len(packed) >= cap:
                break
            data, size, mtime_ns = read(os.path.join(out_dir, name))
            # (the output not yet shown is at most LAG raw bytes; the unit
            # and the manifest's line at least 2x smaller compressed)
            est = sink.out + (min(sink.raw, LAG) + size + 1024 + 150 * (len(packed) + 1)) / 2
            if packed and est > limit:
                break
            add(tf, name, data, mtime_ns / 1e9)
            packed.append((name, sha256(data), size, mtime_ns))
        man = dict(info["manifest"])
        man["units"] = {n: {"sha256": h, "size": s} for n, h, s, _ in packed}
        man["files"] = {n: sha256(d) for n, (d, _) in small.items()}
        add(tf, "manifest.json", json.dumps(man, indent=0).encode(), time.time())
        tf.close()
        sink.finish()
        f.flush()
        os.fsync(f.fileno())
    return packed, os.path.getsize(tmp)


def save_record(path, rec):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(rec, f, indent=0, sort_keys=True)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def host_name():
    """this machine's name; on macOS its LocalHostName, which does not
    change with the network as gethostname() can"""
    h = ""
    if sys.platform == "darwin":
        try:
            h = subprocess.run(["scutil", "--get", "LocalHostName"], stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL, universal_newlines=True,
                               timeout=10).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            h = ""
    h = h or socket.gethostname().split(".")[0]
    h = re.sub(r"[^A-Za-z0-9-]", "-", h).strip("-")
    return h or "host"


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(os.path.dirname(here))
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("dir", nargs="?", default=os.path.join(root, "laptop-runs"))
    # (for tests)
    ap.add_argument("--dest", default=os.path.join(root, "research", "laptop"),
                    help=argparse.SUPPRESS)
    ap.add_argument("--build", default=os.path.join(root, "build-laptop"),
                    help=argparse.SUPPRESS)
    ap.add_argument("--part-mb", type=float, default=PART_MB, help=argparse.SUPPRESS)
    ap.add_argument("--host", default=None, help=argparse.SUPPRESS)
    ap.add_argument("--preset", type=int, default=PRESET, help=argparse.SUPPRESS)
    a = ap.parse_args()
    out_dir = a.dir
    if not os.path.isdir(out_dir):
        sys.exit(f"{out_dir}: no such directory (run.py's --out)")
    limit = a.part_mb * 1e6
    hard = limit * HARD_MB / PART_MB
    host = a.host or host_name()
    os.makedirs(a.dest, exist_ok=True)

    rec_path = os.path.join(out_dir, RECORD)
    try:
        with open(rec_path) as f:
            rec = json.load(f)
    except FileNotFoundError:
        rec = {"units": {}, "tarballs": []}
    except ValueError:
        sys.exit(f"{rec_path} is damaged: remove it to pack every unit again")

    # the units new or changed since they were packed (size and mtime as
    # recorded: unchanged; else compared by sha256)
    todo = []
    changed = 0
    touched = False
    for name in sorted(os.listdir(out_dir)):
        if not UNIT.match(name):
            continue
        r = rec["units"].get(name)
        if r is None:
            todo.append(name)
            continue
        st = os.stat(os.path.join(out_dir, name))
        if st.st_size == r["size"] and st.st_mtime_ns == r.get("mtime_ns"):
            continue
        data, size, mtime_ns = read(os.path.join(out_dir, name))
        if sha256(data) == r["sha256"]:
            r["size"], r["mtime_ns"] = size, mtime_ns
            touched = True
            continue
        todo.append(name)
        changed += 1
    todo.sort(key=lambda n: int(n[1:-6]))
    if not todo:
        if touched:
            save_record(rec_path, rec)
        print(f"nothing new to collect ({len(rec['units'])} units packed before)")
        return

    # the small files, read once: every part carries the same copies
    small = {}
    for name in sorted(os.listdir(out_dir)):
        if SMALL.match(name):
            data, _, mtime_ns = read(os.path.join(out_dir, name))
            small[name] = (data, mtime_ns / 1e9)
    for name in BUILD:
        for path in (os.path.join(a.build, name),
                     os.path.join(out_dir, "build-laptop_" + name)):
            if os.path.exists(path):
                data, _, mtime_ns = read(path)
                small["build-laptop_" + name] = (data, mtime_ns / 1e9)
                break

    # a name not used yet (a second collect in the same minute)
    stamp = time.strftime("%Y%m%d-%H%M")
    existing = os.listdir(a.dest)
    k = 0
    while True:
        base = f"results-{host}-{stamp}" + (f"-{k}" if k else "")
        if not any(x.startswith(base + ".") or x.startswith(base + "-part") for x in existing):
            break
        k += 1
    for x in existing:
        if x.startswith("results-") and x.endswith(".tar.xz.tmp"):
            os.remove(os.path.join(a.dest, x))

    written = []
    i = 0
    part = 0
    t0 = time.time()
    while i < len(todo):
        part += 1
        tmp = os.path.join(a.dest, f"{base}-part{part:02d}.tar.xz.tmp")
        info = {"preset": a.preset,
                "manifest": {"format": 1, "host": host, "collect": base, "part": part,
                             "created": time.time(),
                             "created_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z")}}
        cap = None
        while True:
            try:
                packed, size = write_part(tmp, small, todo[i:], out_dir, limit, cap, info)
            except KeyboardInterrupt:
                os.remove(tmp)
                done = " ".join(w[0] for w in written) or "none"
                sys.exit(f"interrupted; the parts written and recorded (git add them): {done}; "
                         f"run collect again for the rest")
            if size <= hard or len(packed) == 1:
                break
            cap = max(1, len(packed) * 3 // 4)
        if size > hard:
            os.remove(tmp)
            sys.exit(f"{packed[0][0]} alone is {size / 1e6:.1f} MB compressed, over the "
                     f"{hard / 1e6:.1f} MB limit of a part")
        i += len(packed)
        single = part == 1 and i == len(todo)
        name = base + (".tar.xz" if single else f"-part{part:02d}.tar.xz")
        os.replace(tmp, os.path.join(a.dest, name))
        for n, h, s, m in packed:
            rec["units"][n] = {"sha256": h, "size": s, "mtime_ns": m, "tarball": name}
        rec["tarballs"].append({"name": name, "units": len(packed), "bytes": size,
                                "created": info["manifest"]["created_iso"]})
        save_record(rec_path, rec)
        written.append((name, size, len(packed), sum(p[2] for p in packed)))

    raw = sum(w[3] for w in written)
    comp = sum(w[1] for w in written)
    print(f"packed {len(todo)} unit{'s' if len(todo) > 1 else ''}"
          + (f" ({changed} changed since packed before)" if changed else "")
          + f", {raw / 1e6:.1f} MB, into {len(written)} file{'s' if len(written) > 1 else ''} "
          f"of {comp / 1e6:.1f} MB in {time.time() - t0:.0f} s:")
    rel = os.path.relpath(a.dest, root)
    if rel.startswith(".."):
        rel = a.dest
    for name, size, n, _ in written:
        print(f"  {os.path.join(rel, name)} ({size / 1e6:.1f} MB, {n} units)")
    print("(here the counts go to decide.py and into the scheduler's state; the laptop's "
          "times stay out of its time laws)")
    files = " ".join(os.path.join(rel, w[0]) for w in written)
    print(f"to send them back: git add {files} && git commit -m 'laptop run results' && git push")


if __name__ == "__main__":
    main()
