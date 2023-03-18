tup_bank = {}
primes = (2, 3, 5, 7, 11, 13, 17, 19, 23)


def gen_tup(s, n):
    if (s, n) in tup_bank.keys():
        return tup_bank[(s, n)]

    if s == 0:
        ret = [(0,) * n]
        tup_bank[(s, n)] = ret
        return ret

    elif n == 1:
        ret = [(s,)]
        tup_bank[(s, n)] = ret
        return ret

    out = []
    for i in range(0, s + 1):
        prev = gen_tup(s - i, n - 1)
        out += [(i,) + p for p in prev]
    tup_bank[(s, n)] = out
    return out


def gen_rows(pows, n):
    rows = {(1,) * n}
    for p, pow in zip(primes, pows):
        ptups = gen_tup(pow, n)
        new_rows = set()
        for row in rows:
            for tup in ptups:
                new_rows.add(
                    tuple(
                        sorted(
                            [elt * (p**r) for elt, r in zip(row, tup)], reverse=True
                        )
                    )
                )
        rows = new_rows
    return rows


def split_by_sum(all_rows):
    counts = {}
    unique_counts = {}
    for row in all_rows:
        s = sum(row)
        counts[s] = counts.get(s, 0) + 1
        if len(row) == len(set(row)):
            unique_counts[s] = unique_counts.get(s, 0) + 1

    ucs = sorted(unique_counts.items(), key=lambda x: x[1], reverse=True)
    return ucs


def rows_to_rowdict(all_rows):
    ucs = split_by_sum(all_rows)
    S_list = [S for S, c in ucs if c >= 12]
    S_set = set(S_list)
    row_dict = {S: [] for S in S_list}

    print(f">>> total possible 6-vecs: {sum(u[1] for u in ucs)}")
    for row in all_rows:
        if len(set(row)) == len(row) and sum(row) in S_set:
            row_dict[sum(row)] += [row]

    # print(S_list)
    print(f">>> number of S to consider: {len(S_list)}")
    return row_dict
