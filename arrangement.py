import time


def process(input):
    S, filt_rows = input
    start = time.time()

    # print(f"[{S}] (Starting)", flush=True)

    # filt_rows = row_dict[S]
    # sorted([row for row in all_rows if sum(row) == S and len(set(row)) == len(row)], reverse=True)
    l1 = len(filt_rows)
    # print(filt_rows)

    eltc = {}
    eltcl = []

    while True:
        if len(filt_rows) <= 11:
            # print(f"[{S}] (Eliminated by reduction)")
            print(
                f"[{S}] (Completed in {time.time()-start} seconds by reduction)",
                flush=True,
            )
            return 0, [], []
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

    # print(">>> finished in", time.time()-start, "seconds")

    def make_intersections(vecs):
        d = {}
        for i in range(len(vecs)):
            dd = {}
            for j in range(len(vecs)):
                dd[j] = len(set(vecs[i]).intersection(set(vecs[j])))
            d[i] = dd
        return d

    # print(">>> Precomputing intersections")
    # start = time.time()
    intersections = make_intersections(vecs)
    # print(">>> finished in", time.time()-start, "seconds")

    count = [0]
    search_start = time.time()
    sols = []
    nms = []

    def record(rows, cols, ris, cis, count, sols, nms):
        count[0] += 1
        if count[0] % 100000 == 0:
            rs = [tuple(inv_map[e] for e in row) for row in rows]
            cs = [tuple(inv_map[e] for e in col) for col in cols]
            print(
                f"[{S}] (Log) {count[0]} {time.time()-start} {len(vecs)} {sorted(ris+cis)}",
                flush=True,
            )

        score = len(rows) * len(cols)

        if score == 25:
            rs = [tuple(inv_map[e] for e in row) for row in rows]
            cs = [tuple(inv_map[e] for e in col) for col in cols]
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
            # if p1 == pl:
            #   print(f"[{S}] (SOLUTION  6x6) {(rs, cs)}", flush=True)
            #   sols += [(rs, cs)]
            if p1 != pl:
                print(f"[{S}] (Near miss 5x5) {(rs, cs)}", flush=True)
                nms += [(rs, cs)]

        if score == 36:
            rs = [tuple(inv_map[e] for e in row) for row in rows]
            cs = [tuple(inv_map[e] for e in col) for col in cols]
            print(f"[{S}] (SOLUTION  6x6) {(rs, cs)}", flush=True)
            sols += [(rs, cs)]

    def search_aux(i, rows, cols, ris, cis, unmatched):
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

        record(rows, cols, ris, cis, count, sols, nms)

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

    # print(f"[{S}] (Searching)", flush=True)
    search_aux(0, [], [], [], [], set())
    # print(f"[{S}] (Log) {count[0]} {time.time()-search_start} {len(vecs)}")
    print(
        f"[{S}] (Completed in {time.time()-start} seconds: explored {count[0]} configs, {len(nms)} near misses, {len(sols)} solutions)",
        flush=True,
    )

    return count[0], sols, nms
