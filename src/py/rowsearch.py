import argparse
import multiprocessing
import random
import time
from typing import Union

from arrangement import process
from enumeration import gen_rows, rows_to_rowdict


def log(type: str, message: str) -> None:
    print(f">>> ({type}) {message}", flush=True)


def run(N: int, P: tuple[int], S: Union[int, None], deterministic: bool = False) -> None:
    log("enum", "Generating and simplifying rows")
    enum_start = time.time()
    all_rows = gen_rows(P, N)
    row_dict = rows_to_rowdict(all_rows)
    l = sorted(row_dict.items())
    if not deterministic:
        random.shuffle(l)
    log("enum", f"number of S to consider: {len(l)}")
    log("enum", f"finished in {time.time()-enum_start} seconds")

    arrange_start = time.time()
    if S is not None:
        log("arrange", "Starting serial search")
        results = [process(S, row_dict[S])]
    else:
        log("arrange", "Starting parallel search")
        with multiprocessing.Pool(50) as p:
            results = p.starmap(process, l)

    configurations = sum(result.count for result in results)
    solutions = [sol for result in results for sol in result.solutions]
    near_misses = [nm for result in results for nm in result.near_misses]

    log("arrange", f"finished search for P={P} in {time.time()-arrange_start} seconds")

    log("stats", f"total time: {time.time()-enum_start} seconds")
    log("stats", f"total configurations: {configurations}")
    log("stats", f"solutions found: {len(solutions)}")
    log("stats", f"solutions: {solutions}")
    log("stats", f"near misses found: {len(near_misses)}")
    log("stats", f"near misses: {near_misses}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--deterministic", action="store_true")
    parser.add_argument("--sum", metavar="S", type=int, nargs="?")
    parser.add_argument("exponents", metavar="P", type=int, nargs="+")
    args = parser.parse_args()

    N = 6
    P = tuple(args.exponents)
    S = args.sum
    deterministic = args.deterministic
    run(N, P, S, deterministic)
