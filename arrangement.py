import time
from typing import Callable


Vec = tuple[int, int, int, int, int, int]


def make_vecs(
    S: int, filt_rows: list[Vec], start: float
) -> tuple[dict[int, int], list[Vec], int]:
    eltc = {}
    eltcl = []

    while True:
        if len(filt_rows) <= 11:
            # print(f"[{S}] (Eliminated by reduction)")
            print(
                f"[{S}] (Completed in {time.time()-start} seconds by reduction)",
                flush=True,
            )
            return {}, [], 0
        eltc = {}
        for row in filt_rows:
            for elt in row:
                eltc[elt] = eltc.get(elt, 0) + 1
        eltcl = sorted(eltc.items(), key=lambda x: x[1], reverse=True)
        if eltcl[-1][1] >= 2:
            break
        else:
            filt_rows = [row for row in filt_rows if all(eltc[elt] > 1 for elt in row)]

    # print(f"[{S}] (Setup) reduced from {l1} vecs to {len(filt_rows)} vecs", flush=True)

    rename_map = {}
    inv_map = {}
    for i, elt in enumerate(eltcl):
        rename_map[elt[0]] = i
        inv_map[i] = elt[0]

    vecs = sorted(
        [
            tuple(sorted([rename_map[i] for i in row], reverse=True))
            for row in filt_rows
        ],
        reverse=True,
    )
    # for vec in vecs:
    #   print(vec)
    nvecs = len(vecs)
    return inv_map, vecs, nvecs


def make_intersections(vecs: list[Vec]) -> dict[int, dict[int, int]]:
    d = {}
    for i in range(len(vecs)):
        dd = {}
        for j in range(len(vecs)):
            dd[j] = len(set(vecs[i]).intersection(set(vecs[j])))
        d[i] = dd
    return d


def recorder(
    inv_map: dict[int, int], start: float, vecs: list[Vec], S: int, count: list[int]
) -> Callable[
    [list[Vec], list[Vec], list[int], list[int]],
    tuple[tuple[list[Vec], list[Vec]] | None, tuple[list[Vec], list[Vec]] | None],
]:
    def translate(rows: list[Vec], cols: list[Vec]) -> tuple[list[Vec], list[Vec]]:
        rs = [tuple(inv_map[e] for e in row) for row in rows]
        cs = [tuple(inv_map[e] for e in col) for col in cols]
        return rs, cs

    def record(
        rows: list[Vec], cols: list[Vec], ris: list[int], cis: list[int]
    ) -> tuple[tuple[list[Vec], list[Vec]] | None, tuple[list[Vec], list[Vec]] | None]:
        count[0] += 1
        if count[0] % 100000 == 0:
            print(
                f"[{S}] (Log) {count[0]} {time.time()-start} {len(vecs)} {sorted(ris+cis)}",
                flush=True,
            )

        score = len(rows) * len(cols)
        rs, cs = translate(rows, cols)

        if score == 25:
            missing_elt = 6 * S - sum(set().union(*rs).union(*cs))
            last_row = list(set().union(*rs).difference(set().union(*cs))) + [
                missing_elt
            ]
            pl = 1
            p1 = 1
            for e1, el in zip(rs[0], last_row):
                pl = pl * el
                p1 = p1 * e1
            # print(e, p1, pl, flush=True)
            if p1 == pl:
                print(f"[{S}] (SOLUTION  6x6) {(rs, cs)}", flush=True)
                return ((rs, cs), None)
            if p1 != pl:
                print(f"[{S}] (Near miss 5x5) {(rs, cs)}", flush=True)
                return (None, (rs, cs))

        if score > 25:
            raise Exception(f"{score} should never happen?")
        return (None, None)

    return record


def searcher(
    nvecs: int,
    vecs: list[Vec],
    record: Callable[
        [list[Vec], list[Vec], list[int], list[int]],
        tuple[tuple[list[Vec], list[Vec]] | None, tuple[list[Vec], list[Vec]] | None],
    ],
    sols: list[tuple[list[Vec], list[Vec]]],
    nms: list[tuple[list[Vec], list[Vec]]],
    intersections: dict[int, dict[int, int]],
) -> Callable[[int, list[Vec], list[Vec], list[int], list[int], set[int]], None]:
    def search_aux(
        i: int,
        rows: list[Vec],
        cols: list[Vec],
        ris: list[int],
        cis: list[int],
        unmatched: set[int],
    ):
        nonlocal sols
        nonlocal nms
        # print(i, rows, cols)
        if i < nvecs:
            mr = vecs[i][0]
            if unmatched and max(unmatched) > mr:
                return
            # union_rows = set().union(*rows)
            # union_cols = set().union(*cols)
            # if not union_cols.issuperset(sorted(e for e in union_rows if e > mr)):
            #   return
            # if not union_rows.issuperset(sorted(e for e in union_cols if e > mr)):
            #   return

        nm, sol = record(rows, cols, ris, cis)
        if nm:
            nms += [nm]
        if sol:
            sols += [sol]

        for j in range(i, nvecs):
            vec = vecs[j]

            is_valid_row = all(intersections[j][ri] == 0 for ri in ris) and all(
                intersections[j][ci] == 1 for ci in cis
            )
            if is_valid_row:
                um = unmatched.symmetric_difference(vec)
                search_aux(j + 1, rows + [vec], cols, ris + [j], cis, um)

            is_valid_col = (
                len(rows) > 0
                and all(intersections[j][ri] == 1 for ri in ris)
                and all(intersections[j][ci] == 0 for ci in cis)
            )
            if is_valid_col:
                um = unmatched.symmetric_difference(vec)
                search_aux(j + 1, rows, cols + [vec], ris, cis + [j], um)

    return search_aux


def process(
    input: tuple[int, list[Vec]]
) -> tuple[int, list[tuple[list[Vec], list[Vec]]], list[tuple[list[Vec], list[Vec]]]]:
    S, filt_rows = input
    start = time.time()

    inv_map, vecs, nvecs = make_vecs(S, filt_rows, start)
    if not inv_map:
        return 0, [], []

    # print(f"[{S}] (Starting)", flush=True)
    # print(">>> finished in", time.time()-start, "seconds")

    # print(">>> Precomputing intersections")
    # start = time.time()
    intersections = make_intersections(vecs)
    # print(">>> finished in", time.time()-start, "seconds")

    count = [0]
    search_start = time.time()
    sols = []
    nms = []

    record = recorder(inv_map, start, vecs, S, count)
    search_aux = searcher(nvecs, vecs, record, sols, nms, intersections)

    # print(f"[{S}] (Searching)", flush=True)
    search_aux(0, [], [], [], [], set())
    # print(f"[{S}] (Log) {count[0]} {time.time()-search_start} {len(vecs)}")
    print(
        f"[{S}] (Completed in {time.time()-start} seconds: explored {count[0]} configs, {len(nms)} near misses, {len(sols)} solutions)",
        flush=True,
    )

    return count[0], sols, nms
