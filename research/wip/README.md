# Work in progress (saved 2026-10-08 before a usage-limit stop)

Not integrated yet; kept here so nothing is lost if the session's container
is reclaimed. Each directory under patches/ applies with `git am`:

| patches/ | base | what |
|---|---|---|
| cx-profile, cx-decomp, cx-magic, cx-prune | fdb77fc | instrumentation and negative-result prototypes of the complexity study |
| c2-classhall | be625d8 | round 2: the V_d profile's DPROF instrumentation (c2/profile), the class-matching support prototype (c2/vdclass, `-DCLASS_SUP`) and the stage-A counts of its pair rules (research/ideas.md, "Pair rules on top of the class support") |

Notes: complexity-profile.md (where the N^4 comes from), complexity-study.md
(ideas, picks, prototype measurements, verifier reports).
Scheduler v2 has since been merged (55855a9; research/scheduler-v2.md).
cx-wide, cx-pretest and cx-dfirst have since been integrated (branch
cx/integrated, with the reviewers' fixes, rebased onto scheduler v2;
research/ideas.md, "Integration of cx/wide, cx/pretest and cx/dfirst"), so
their patches and the integration's snapshot are gone from here.

data/all_records.tsv.xz: every distinct msearch output record of this
session's runs (172,701 sum records, 11,547 square records, done/dchunk/dsum
records), one per line as `source-file<TAB>json`, deduplicated by line. The
source path tells which run a record came from (e.g. collect/ = the 36 seed
P, forecast/state/units/ = the scheduler's 40 units, existence/ = large-N
and sampled runs; sampled runs only searched every k-th first row).

Round 2 (complexity of the d-first search; patches against be625d8):
c2-classhall (the instrumentation of the pair-rules no-go), c2-profile and
c2-vdclass (instrumentation), c2-partner (partner pruning: negative).
Profile: complexity-dfirst-profile.md. c2/classsup (class-matching support
in V_d, verified: 0.66x CPU at N >= 20k), c2/star (star cover of the d loop,
verified: 1.045-1.065x) and c2/classhall's write-up (pair rules: no-go) have
since been integrated (branch integ/round2: research/ideas.md, "Integration
of the round-2 d-first changes"), so their patches are gone from here.
