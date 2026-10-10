#!/usr/bin/env python3
"""Tests of scripts/machine_cal.py that run no search: the ISA and core
parsers on canned /proc/cpuinfo (Sapphire Rapids, Zen 3, Zen 4, Cascade
Lake, a Graviton) and sysctl (an M1) text, the bench parser, the node
checks, and the speed arithmetic."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import machine_cal as mc  # noqa: E402

COMMON = ("fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx "
          "fxsr sse sse2 ht syscall nx pdpe1gb rdtscp lm constant_tsc rep_good nopl xtopology "
          "pni pclmulqdq ssse3 fma cx16 sse4_1 sse4_2 x2apic movbe popcnt aes xsave avx f16c "
          "rdrand hypervisor lahf_lm abm bmi1 avx2 bmi2 erms adx smap clflushopt sha_ni")
SPR = COMMON + (" avx512f avx512dq rdseed avx512ifma avx512cd avx512bw avx512vl xsaveopt "
                "avx_vnni avx512_bf16 avx512vbmi umip pku ospke waitpkg avx512_vbmi2 gfni vaes "
                "vpclmulqdq avx512_vnni avx512_bitalg tme avx512_vpopcntdq la57 rdpid "
                "bus_lock_detect cldemote movdiri movdir64b enqcmd fsrm md_clear serialize "
                "tsxldtrk pconfig arch_lbr amx_bf16 avx512_fp16 amx_tile amx_int8 flush_l1d")
CLX = COMMON + (" hle rtm mpx avx512f avx512dq rdseed avx512cd avx512bw avx512vl xsaveopt "
                "xsavec xgetbv1 xsaves ida arat pku ospke avx512_vnni md_clear flush_l1d "
                "arch_capabilities")
ZEN3 = COMMON + (" mmxext fxsr_opt cr8_legacy sse4a misalignsse 3dnowprefetch osvw topoext "
                 "perfctr_core ssbd ibrs ibpb stibp vmmcall fsgsbase invpcid rdseed clzero "
                 "xsaveerptr rdpru wbnoinvd arat npt nrip_save umip pku ospke vaes "
                 "vpclmulqdq rdpid fsrm")
ZEN4 = ZEN3 + (" avx512f avx512dq avx512ifma avx512cd avx512bw avx512vl avx512_bf16 "
               "avx512vbmi avx512_vbmi2 gfni avx512_vnni avx512_bitalg avx512_vpopcntdq "
               "flush_l1d")


def cpuinfo(vendor, model, flags, cpus):
    """canned /proc/cpuinfo: cpus = [(physical id, core id)], one block each"""
    nper = {}
    for p, c in cpus:
        nper.setdefault(p, set()).add(c)
    sib = {p: sum(1 for q, _ in cpus if q == p) for p in nper}
    out = []
    for i, (p, c) in enumerate(cpus):
        out.append(f"processor\t: {i}\nvendor_id\t: {vendor}\ncpu family\t: 6\nmodel\t\t: 1\n"
                   f"model name\t: {model}\nphysical id\t: {p}\nsiblings\t: {sib[p]}\n"
                   f"core id\t\t: {c}\ncpu cores\t: {len(nper[p])}\nflags\t\t: {flags}\n"
                   f"bogomips\t: 4200.00\n")
    return "\n".join(out) + "\n"


SMT4 = [(0, 0), (0, 1), (0, 0), (0, 1)]    # 2 cores x 2 threads
NOSMT4 = [(0, 0), (0, 1), (0, 2), (0, 3)]

GRAVITON = "".join(
    f"processor\t: {i}\nBogoMIPS\t: 243.75\nFeatures\t: fp asimd evtstrm aes pmull sha1 sha2 "
    f"crc32 atomics fphp asimdhp cpuid asimdrdm lrcpc dcpop asimddp ssbs\nCPU implementer\t: "
    f"0x41\nCPU architecture: 8\nCPU variant\t: 0x3\nCPU part\t: 0xd0c\nCPU revision\t: 1\n\n"
    for i in range(4))

M1_SYSCTL = """hw.ncpu: 8
hw.byteorder: 1234
hw.activecpu: 8
hw.physicalcpu: 8
hw.physicalcpu_max: 8
hw.logicalcpu: 8
hw.logicalcpu_max: 8
hw.perflevel0.physicalcpu: 4
hw.perflevel0.logicalcpu: 4
hw.perflevel0.name: Performance
hw.perflevel1.physicalcpu: 4
hw.perflevel1.name: Efficiency
hw.optional.floatingpoint: 1
hw.optional.neon: 1
hw.optional.arm64: 1
hw.optional.armv8_crc32: 1
hw.optional.arm.FEAT_LSE: 1
machdep.cpu.brand_string: Apple M1
machdep.cpu.core_count: 8
machdep.cpu.thread_count: 8
"""

LSCPU_SMT = """Architecture:                            x86_64
CPU(s):                                  4
Thread(s) per core:                      2
Core(s) per socket:                      2
Socket(s):                               1
"""

FAST = ["avx2", "avx512bitalg", "avx512bw", "avx512f", "avx512vbmi", "avx512vpopcntdq", "bmi2",
        "gfni"]


def test_cpuinfo():
    spr = mc.parse_cpuinfo(cpuinfo("GenuineIntel", "Intel(R) Xeon(R) Platinum 8481C CPU @ 2.70GHz",
                                   SPR, NOSMT4))
    assert spr["isa"] == FAST, spr
    assert spr["logical"] == 4 and spr["physical"] == 4 and spr["arch"] == "x86_64"
    assert mc.isa_class(spr["isa"]).startswith("avx512 full")
    zen4 = mc.parse_cpuinfo(cpuinfo("AuthenticAMD", "AMD EPYC 9R14 96-Core Processor", ZEN4, SMT4))
    assert zen4["isa"] == FAST, zen4
    assert zen4["logical"] == 4 and zen4["physical"] == 2
    assert mc.isa_class(zen4["isa"]).startswith("avx512 full")
    zen3 = mc.parse_cpuinfo(cpuinfo("AuthenticAMD", "AMD EPYC 7R13 Processor", ZEN3, SMT4))
    assert zen3["isa"] == ["avx2", "bmi2"], zen3
    assert zen3["physical"] == 2 and zen3["logical"] == 4
    assert mc.isa_class(zen3["isa"]).startswith("avx2")
    clx = mc.parse_cpuinfo(cpuinfo("GenuineIntel", "Intel(R) Xeon(R) Platinum 8259CL CPU @ 2.50GHz",
                                   CLX, SMT4))
    assert clx["isa"] == ["avx2", "avx512bw", "avx512f", "bmi2"], clx
    cls = mc.isa_class(clx["isa"])
    assert "cascadelake" in cls and "avx512vpopcntdq/avx512vbmi/gfni/avx512bitalg" in cls, cls
    arm = mc.parse_cpuinfo(GRAVITON)
    assert arm["arch"] == "arm64" and arm["isa"] == ["arm64"] and arm["logical"] == 4
    assert arm["physical"] is None
    assert mc.isa_class(arm["isa"]).startswith("arm64")
    # this machine's own text parses (whatever it is)
    if os.path.exists("/proc/cpuinfo"):
        own = mc.parse_cpuinfo(open("/proc/cpuinfo").read())
        assert own["logical"] >= 1
    print("cpuinfo ok")


def test_lscpu_sysctl():
    ls = mc.parse_lscpu(LSCPU_SMT)
    assert ls == {"threads_per_core": 2, "cores_per_socket": 2, "sockets": 1}, ls
    m1 = mc.parse_sysctl(M1_SYSCTL)
    assert m1["arch"] == "arm64" and m1["isa"] == ["arm64"], m1
    assert m1["logical"] == 8 and m1["physical"] == 8 and m1["perf"] == 4
    assert m1["model"] == "Apple M1"
    intel_mac = mc.parse_sysctl("hw.physicalcpu: 8\nhw.logicalcpu: 16\nhw.optional.avx2_0: 1\n"
                                "hw.optional.avx512f: 1\nhw.optional.avx512bw: 1\n"
                                "machdep.cpu.leaf7_features: RDWRFSGS TSC_THREAD_OFFSET BMI1 AVX2 "
                                "BMI2 AVX512F AVX512DQ AVX512CD AVX512BW AVX512VL\n"
                                "machdep.cpu.brand_string: Intel(R) Core(TM) i9-10910 CPU\n")
    assert intel_mac["isa"] == ["avx2", "avx512bw", "avx512f", "bmi2"], intel_mac
    assert intel_mac["physical"] == 8 and intel_mac["logical"] == 16
    print("lscpu / sysctl ok")


def test_bench_and_checks():
    text = """instance                   N labels        nodes   sq check   time(s) Mnodes/s
6x6 P=10.4.3.2 S=327     451    70         8688    1    ok    0.0010     8.76
6x6 P=10.6.3.1.0.1 S=648  1682   120       960236    1    ok    0.1556     6.17
TOTAL                                   1770779         ok    0.2543     6.96
"""
    assert mc.parse_bench(text) == (1770779, True, 2)
    bad = text.replace("1    ok    0.0010", "1  FAIL    0.0010").replace("ok    0.2543", "FAIL    0.2543")
    assert mc.parse_bench(bad)[1] is False
    t1 = mc.TESTS[0]
    rec = {"type": "sum", "nodes": mc.REFS["carry512"]["T1"], "squares": 18}
    assert mc.check_test(t1, rec, "carry512", mc.REFS) is None
    assert "nodes" in mc.check_test(t1, dict(rec, nodes=1), "carry512", mc.REFS)
    assert "squares" in mc.check_test(t1, dict(rec, squares=17), "carry512", mc.REFS)
    assert mc.check_test(t1, dict(rec, nodes=mc.REFS["matrix"]["T1"]), "matrix", mc.REFS) is None
    t3 = mc.TESTS[2]
    assert "pairs" in mc.check_test(t3, {"nodes": mc.REFS["carry512"]["T3"], "pairs": 1},
                                    "carry512", mc.REFS)
    # the bench gate refuses an unknown path
    try:
        mc.exactness("/nonexistent", "carry256", mc.REFS, say=lambda *a: None)
        raise AssertionError("no failure")
    except SystemExit as e:
        assert "no reference" in str(e)
    print("bench / checks ok")


def test_speeds():
    # K copies whose mean wall is the reference CPU / s run at K x s
    # reference CPU-seconds per second
    rows = {t["name"]: {"wall": t["ref_cpu"] / 0.5, "cpu": t["ref_cpu"] / 0.4} for t in mc.TESTS}
    sp, pp = mc.band(rows, ["T1", "T2"], 4)
    assert abs(sp - 2.0) < 1e-12 and abs(pp - 0.4) < 1e-12
    sp, pp = mc.band(rows, ["T3", "T4"], 1)
    assert abs(sp - 0.5) < 1e-12
    print("speeds ok")


if __name__ == "__main__":
    test_cpuinfo()
    test_lscpu_sysctl()
    test_bench_and_checks()
    test_speeds()
    print("ok")
