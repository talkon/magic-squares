import argparse
import multiprocessing
import random
import time
from typing import Union

from arrangement import process
from enumeration import gen_rows, rows_to_rowdict


def log(type: str, message: str) -> None:
    print(f">>> ({type}) {message}", flush=True)


def run(
    N: int,
    P: tuple[int, ...],
    SS: Union[list[int], None] = None,
    *,
    deterministic: bool = False,
    serial: bool = False,
) -> None:
    log("enum", "Generating and simplifying rows")
    enum_start = time.time()
    all_rows = gen_rows(P, N)
    row_dict = rows_to_rowdict(all_rows)
    sp_sets = []
    for S, vecs in row_dict.items():
        if (not SS) or (S in SS):
            sp_sets.append((S, vecs))
    if deterministic:
        sp_sets.sort()
    else:
        random.shuffle(sp_sets)
    log("enum", f"number of S to consider: {len(sp_sets)}")
    log("enum", f"finished in {time.time()-enum_start} seconds")

    arrange_start = time.time()
    try:
        if serial:
            log("arrange", "Starting serial search")
            results = [process(S, vecs) for S, vecs in sp_sets]
        else:
            log("arrange", "Starting parallel search")
            with multiprocessing.Pool(50) as p:
                results = p.starmap(process, sp_sets)
    finally:
        log("arrange", f"finished in {time.time()-arrange_start} seconds")
        log("stats", f"total time: {time.time()-enum_start} seconds")

    configurations = sum(result.count for result in results)
    solutions = [sol for result in results for sol in result.solutions]
    near_misses = [nm for result in results for nm in result.near_misses]
    log("stats", f"total configurations: {configurations}")
    log("stats", f"solutions found: {len(solutions)}")
    log("stats", f"solutions: {solutions}")
    log("stats", f"near misses found: {len(near_misses)}")
    log("stats", f"near misses: {near_misses}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--deterministic", action="store_true")
    parser.add_argument("--serial", action="store_true")
    parser.add_argument("--sum", metavar="S", type=int, nargs="*")
    parser.add_argument("exponents", metavar="P", type=int, nargs="+")
    args = parser.parse_args()

    N = 6
    P = tuple(args.exponents)
    SS = args.sum
    run(N, P, SS, deterministic=args.deterministic, serial=args.serial)
