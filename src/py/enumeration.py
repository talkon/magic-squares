from collections import Counter
import argparse
from typing import Union

Vec = tuple[int, int, int, int, int, int]

tup_bank = {}
primes = (2, 3, 5, 7, 11, 13, 17, 19, 23)


def gen_tup(s: int, n: int) -> list[tuple[int, ...]]:
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


def gen_rows(pows: tuple[int, ...], n: int) -> set[Vec]:
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


def split_by_sum(all_rows: set[Vec]) -> list[tuple[int, int]]:
    counts = {}
    unique_counts = {}
    for row in all_rows:
        s = sum(row)
        counts[s] = counts.get(s, 0) + 1
        if len(row) == len(set(row)):
            unique_counts[s] = unique_counts.get(s, 0) + 1

    ucs = sorted(unique_counts.items(), key=lambda x: x[1], reverse=True)
    return ucs


def rows_to_rowdict(all_rows: set[Vec]) -> dict[int, list[Vec]]:
    ucs = split_by_sum(all_rows)
    S_set = {S for S, c in ucs if c >= 12}
    row_dict = {S: [] for S in S_set}

    print(f">>> (enum) total possible 6-vecs: {sum(u[1] for u in ucs)}")
    for row in all_rows:
        if len(set(row)) == len(row) and sum(row) in S_set:
            row_dict[sum(row)].append(row)

    # reduction: if there's less than 12 vecs, we can't fill the grid
    for S, vecs in row_dict.items():
        while len(vecs) >= 12:
            count = Counter([elt for vec in vecs for elt in vec])
            if count.most_common()[-1][1] >= 2:
                break
            # reduction: we can't use vecs that have unique elts
            vecs = [vec for vec in vecs if all(count[elt] > 1 for elt in vec)]
        else:
            vecs = []
        row_dict[S] = vecs

    row_dict = {S: vecs for S, vecs in row_dict.items() if vecs}

    return row_dict


def write_enums(row_dict: dict[int, list[Vec]], file: str, SS: Union[list[int], None] = None):
    # rows are guaranteed to written be in ascending order of sum
    n_rows = 0
    with open(file, "w") as f:
        for s in sorted(row_dict.keys()):
            if (not SS) or (s in SS):
                for row in row_dict[s]:
                    print(" ".join(str(i) for i in row), file=f)
                    n_rows += 1
    print(f">>> wrote {n_rows} rows to {file} in ascending order of sum")

def write_counts(row_dict: dict[int, list[Vec]], file: str, SS: Union[list[int], None] = None):
    counts = {}
    for s in sorted(row_dict.keys()):
        if (not SS) or (s in SS):
            counts[s] = len(row_dict[s])
    with open(file, "w") as f:
        for s, c in sorted(counts.items(), key=lambda x: x[1], reverse=True):
            print(f"{s} {c}", file=f)
    print(f">>> wrote counts to {file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", metavar="file", type=str)
    parser.add_argument("--sum", metavar="S", type=int, nargs="*")
    parser.add_argument("-c", action="store_true")
    parser.add_argument("exponents", metavar="P", type=int, nargs="+")
    args = parser.parse_args()

    P = tuple(args.exponents)
    N = 6
    SS = args.sum
    all_rows = gen_rows(P, N)
    row_dict = rows_to_rowdict(all_rows)
    if args.c:
        write_counts(row_dict, args.file, SS)
    else:
        write_enums(row_dict, args.file, SS)
