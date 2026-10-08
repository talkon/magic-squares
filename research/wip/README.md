# Work in progress (saved 2026-10-08 before a usage-limit stop)

Not integrated yet; kept here so nothing is lost if the session's container
is reclaimed. Each directory under patches/ applies with `git am`:

| patches/ | base | what |
|---|---|---|
| cx-wide | fdb77fc | carried bitsets up to 512 labels, per-r1 width (verified: 2.8-4x above 256 labels) |
| cx-pretest | fdb77fc | count-only pretest of deep children (verified: 1.1-1.2x at N >= 15k) |
| cx-dfirst | fdb77fc | msearch --diag-first, r1/d sampling (verified: 0.25-0.35x plain CPU at N 15-23k) |
| cx-integrated | a3812f0 | the three above merged, with the reviewers' fixes (integration in progress) |
| cx-profile, cx-decomp, cx-magic, cx-prune | fdb77fc | instrumentation and negative-result prototypes of the complexity study |
| scheduler-v2 | a3812f0 | scheduler with the analytic E(P, S) model and a wider pool (in progress; uncommitted.diff is the worktree's state on top of the commits) |

Notes: complexity-profile.md (where the N^4 comes from), complexity-study.md
(ideas, picks, prototype measurements, verifier reports),
scheduler-v2-design.md (designs and the merged plan).

data/all_records.tsv.xz: every distinct msearch output record of this
session's runs (172,701 sum records, 11,547 square records, done/dchunk/dsum
records), one per line as `source-file<TAB>json`, deduplicated by line. The
source path tells which run a record came from (e.g. collect/ = the 36 seed
P, forecast/state/units/ = the scheduler's 40 units, existence/ = large-N
and sampled runs; sampled runs only searched every k-th first row).
