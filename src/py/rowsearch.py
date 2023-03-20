import argparse
import multiprocessing
import random
import sys
import time

from arrangement import process
from enumeration import gen_rows, rows_to_rowdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sum", metavar="S", type=int, nargs="?")
    parser.add_argument("exponents", metavar="P", type=int, nargs="+")
    args = parser.parse_args()

    P = tuple(args.exponents)
    N = 6

    print(">>> Generating and simplifying rows")
    start = time.time()

    all_rows = gen_rows(P, N)
    row_dict = rows_to_rowdict(all_rows)
    if args.sum is not None:
        row_dict = {S: vecs for S, vecs in row_dict.items() if S == args.sum}

    print(f">>> number of S to consider: {len(row_dict)}")
    print(f">>> finished in {time.time()-start} seconds", flush=True)

    print(">>> Starting parallel search", flush=True)
    l = list(row_dict.items())
    random.shuffle(l)
    search_start = time.time()
    with multiprocessing.Pool(50) as p:
        solss = p.map(process, l)
        # solss = [process((327, row_dict[327]))]
        sols = []
        nms = []
        count = 0
        for c, ss, ns in solss:
            count += c
            if ss:
                sols += ss
            if ns:
                nms += ns
        print(
            f">>> finished search for P = {P} in {time.time()-search_start} seconds",
            flush=True,
        )
        print(f">>> total time: {time.time()-start} seconds", flush=True)
        print(f">>> total configurations: {count}", flush=True)
        print(f">>> solutions found: {len(sols)}", flush=True)
        print(f">>> solutions: {sols}", flush=True)
        print(f">>> near misses found: {len(nms)}", flush=True)
        print(f">>> near misses: {nms}", flush=True)
