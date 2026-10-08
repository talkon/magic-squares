# Work in progress (saved 2026-10-08 before a usage-limit stop)

Not integrated yet; kept here so nothing is lost if the session's container
is reclaimed. Each directory under patches/ applies with `git am`:

| patches/ | base | what |
|---|---|---|
| cx-profile, cx-decomp, cx-magic, cx-prune | fdb77fc | instrumentation and negative-result prototypes of the complexity study |

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
