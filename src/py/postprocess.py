import dataclasses
import argparse
from operator import mul
from functools import reduce
from itertools import permutations
from enum import IntEnum
from typing import Union, Callable, Any
from glob import glob, escape
from sys import stderr

from arrangement import Table, Vec
from enumeration import primes

class DiagonalType(IntEnum):
    NONE = 0
    S = 2
    P = 3
    SP = 7

groups: tuple[tuple[tuple[int]]] = (
    ((0, 1), (2, 3), (4, 5)),
    ((0, 1), (2, 4), (3, 5)),
    ((0, 1), (2, 5), (3, 4)),
    ((0, 2), (1, 3), (4, 5)),
    ((0, 2), (1, 4), (3, 5)),
    ((0, 2), (1, 5), (3, 4)),
    ((0, 3), (1, 2), (4, 5)),
    ((0, 3), (1, 4), (2, 5)),
    ((0, 3), (1, 5), (2, 4)),
    ((0, 4), (1, 2), (3, 5)),
    ((0, 4), (1, 3), (2, 5)),
    ((0, 4), (1, 5), (2, 3)),
    ((0, 5), (1, 2), (3, 4)),
    ((0, 5), (1, 3), (2, 4)),
    ((0, 5), (1, 4), (2, 3))
)

@dataclasses.dataclass
class DiagonalStats:
    S: int
    P: Union[tuple[int], None]
    P_val: int
    count_so_far: int
    max_nb: int
    S_count: int
    P_count: int
    SP_count: int
    score_counts: dict[int, int]
    best_score: int
    best_square: tuple[tuple[int]]

    @classmethod
    def pretty_print_header(cls) -> None:
        print(f"  {'P_tup':14} {'count_so_far':>13} {'P_val':>12} {'S':>4} {'max':>4}" +
              f" {'#S':>3} {'#P':>3} {'#SP':>3} {'best':>4}  {'score_counts'}")

    def pretty_print(self) -> None:
        P_str = ' '.join(str(x) for x in self.P) if self.P else "None"
        P_val = self.P_val if self.P_val else "None"
        print(f"  {P_str:14} {self.count_so_far:>13} {P_val:12} {self.S:4} {self.max_nb:4}" + 
              f" {self.S_count:3} {self.P_count:3} {self.SP_count:3} {self.best_score:4}  {sorted(self.score_counts.items())}")
    
    def pretty_print_best_square(self) -> None:
        P_str = ' '.join(str(x) for x in self.P) if self.P else "None"
        P_val = self.P_val if self.P_val else "None"
        print(f"  {'P_tup':<7}={P_str:>13}")
        print(f"  {'P_val':<7}={P_val:>13}")
        print(f"  {'S':<7}={self.S:>13}")
        print(f"  {'max_nb':<7}={self.max_nb:>13}")
        print(f"  {'score':<7}={self.best_score:>13}")
        print(f"  {'square':<7}=")
        for row in self.best_square:
            print('    ' + ' '.join([f"{i:4}" for i in row]))
@dataclasses.dataclass
class SStats:
    S: int
    num_vecs: int
    count: int
    num_sols: int
    best_score: int   

    @classmethod
    def pretty_print_header(cls) -> None:
        print(f"{'S':>9} {'nvecs':>5} {'count':>12} {'nsols':>5} {'best':>4}")

    def pretty_print(self) -> None:  
        print(f"{self.S:>9} {self.num_vecs:5} {self.count:12} {(self.num_sols if self.num_sols > 0 else ''):5} {(self.best_score if self.best_score >= 0 else ''):4}")

@dataclasses.dataclass
class PStats:
    P: Union[tuple[int], None]
    P_val: Union[int, None]
    S_stats: list[SStats]
    count: int = 0
    time: float = 0
    num_vecs: int = 0
    max_num_vecs: int = 0
    num_sols: int = 0
    best_score: int = -1
    num_S: int = 0
    max_S: int = 0
    total_S_count: int = 0
    total_P_count: int = 0
    total_SP_count: int = 0

    def update(self, diagonal_stats: DiagonalStats):
        self.num_sols += 1
        self.best_score = max(self.best_score, diagonal_stats.best_score)
        self.S_stats[-1].num_sols += 1
        self.S_stats[-1].best_score = max(self.S_stats[-1].best_score, diagonal_stats.best_score)
        self.total_S_count += diagonal_stats.S_count
        self.total_P_count += diagonal_stats.P_count
        self.total_SP_count += diagonal_stats.SP_count

    @classmethod
    def pretty_print_header(cls) -> None:
        print(f"  {'P':14} {'P_val':>12} {'num_vecs':>8} {'max':>5} {'max_S':>5} {'count':>13}" +
              f" {'time':>10} {'sols':>5} {'#S':>4} {'#P':>4} {'#SP':>4} {'best':>4}")

    def pretty_print(self) -> None:
        P_str = ' '.join(str(x) for x in self.P) if self.P else "None"
        P_val = self.P_val if self.P_val else "None"
        max_S = min(99999, self.max_S)
        print(f"  {P_str:14} {P_val:>12} {self.num_vecs:8} {self.max_num_vecs:5} {max_S:5} {self.count:13}" +
              f" {self.time:>10.2f} {self.num_sols:5} {self.total_S_count:4} {self.total_P_count:4} {self.total_SP_count:4}" + 
              f" {(self.best_score if self.best_score >= 0 else ''):4}")

    def pretty_print_SStats(self, verbose: int) -> None:
        if verbose >= 5:
            print("  Statistics for each S with a solution:")
            SStats.pretty_print_header()
            for S_stat in self.S_stats:
                if S_stat.num_sols > 0:
                    S_stat.pretty_print()
        if verbose >= 6:
            print("  Statistics for each S:")
            SStats.pretty_print_header()
            for S_stat in self.S_stats:
                S_stat.pretty_print()

@dataclasses.dataclass
class OverallStats:
    count: int
    time: float
    num_P: int
    num_S: int
    num_vecs: int
    max_num_vecs: int
    num_sols: int

    def pretty_print(self) -> None:
        print("Overall statistics:")
        print(f"  {'count':<13}={self.count:14}")
        print(f"  {'time':<13}={self.time:14.2f}")
        print(f"  {'num_P':<13}={self.num_P:14}")
        print(f"  {'num_S':<13}={self.num_S:14}")
        print(f"  {'num_vecs':<13}={self.num_vecs:14}")
        print(f"  {'max_num_vecs':<13}={self.max_num_vecs:14}")
        print(f"  {'num_sols':<13}={self.num_sols:14}")

@dataclasses.dataclass
class Solution:
    S: int
    P: int
    table: Table

    def __init__(self, table: Table):
        self.table = table
        self.S = sum(table.rows[0])
        self.P = reduce(mul, table.rows[0], 1)
        self.assert_consistent()

    def assert_consistent(self) -> None:
        assert len(self.table.rows) == 6, f"expected 6 rows in solution, found {len(self.table.rows)}"
        assert len(self.table.cols) == 6, f"expected 6 cols in solution, found {len(self.table.rows)}"

        for rc in self.table.rows + self.table.cols:
            assert sum(rc) == self.S, f"expected S={self.S}, but sum of {rc} is {sum(rc)}"
            assert reduce(mul, rc, 1) == self.P, f"expected P={self.S}, but product of {rc} is {sum(rc)}"
        
        elts = set()
        for r in self.table.rows:
            elts.update(r)
        assert len(elts) == 36, f"expected 36 distinct numbers in table, found {len(elts)}"
        
        for r in self.table.rows:
            for c in self.table.cols:
                num_intersect = len(set(r).intersection(c))
                assert num_intersect == 1, f"row {r} and col {c} has {num_intersect} elements in common, expected 1"

    def diagonal_stats(self, P_tup: tuple[int], count_so_far: int) -> DiagonalStats:
        rows = self.table.rows
        cols = self.table.cols

        best_score = -1
        best_perm_group = None

        scores = {}
        score_counts = {}

        S_count = 0
        P_count = 0
        SP_count = 0

        # score all diagonals
        for perm in permutations(range(6)):
            diag = tuple(list(set(rows[i]).intersection(cols[j]))[0] for i, j in enumerate(perm))
            correct_sum = (sum(diag) == self.S)
            correct_product = (reduce(mul, diag, 1) == self.P)
            if correct_sum and correct_product:
                scores[perm] = int(DiagonalType.SP)
                SP_count += 1
                S_count += 1
                P_count += 1
            elif correct_product:
                scores[perm] = int(DiagonalType.P)
                P_count += 1
            elif correct_sum:
                scores[perm] = int(DiagonalType.S)
                S_count += 1
            else:
                scores[perm] = int(DiagonalType.NONE)
      
        # count scores of diagonal pairs
        for p1 in permutations(range(6)):
            for g, p2 in zip(groups, other_diagonals(p1)):
                score = scores[p1] + scores[p2]
                score_counts[score] = score_counts.get(score, 0) + 1
                if score > best_score:
                    best_score = score
                    best_perm_group = (p1, g)
        
        # compute best square
        p, g = best_perm_group
        rc_to_perm_row = [g[0][0], g[1][0], g[2][0], g[2][1], g[1][1], g[0][1]]
        # reordered_rows = [rows[rc_to_perm_row[i]] for i in range(6)]
        # reordered_cols = [cols[p[rc_to_perm_row[i]]] for i in range(6)]
        best_square = tuple(
            tuple(
                list(set(rows[rc_to_perm_row[i]]).intersection(cols[p[rc_to_perm_row[j]]]))[0]
                for j in range(6)
            )
            for i in range(6)
        )
        
        max_nb = max(max(row) for row in rows)

        return DiagonalStats(self.S, P_tup, self.P, count_so_far, max_nb, S_count, P_count, SP_count, score_counts, best_score, best_square)

@dataclasses.dataclass
class CSearchStats:
    overall_stats: Union[OverallStats, None]
    solutions: list[Solution]
    diagonal_stats: list[DiagonalStats]
    P_stats: dict[Union[tuple[int], None], PStats]

    def assert_num_sols(self, expected_num_sols: int) -> None:
        self.compute_overall_stats()
        assert self.overall_stats.count > 0, f"Expected positive search count"
        assert self.overall_stats.num_sols == expected_num_sols, f"Expected {expected_num_sols} solutions, found {self.overall_stats.num_sols}"

    def compute_overall_stats(self) -> None:
        overall_stats = OverallStats(0, 0, 0, 0, 0, 0, 0)
        overall_stats.num_P = len(self.P_stats)
        for P_stat in self.P_stats.values():
            overall_stats.count += P_stat.count
            overall_stats.time += P_stat.time
            overall_stats.num_S += P_stat.num_S
            overall_stats.num_vecs += P_stat.num_vecs
            overall_stats.max_num_vecs = max(overall_stats.max_num_vecs, P_stat.max_num_vecs)
            overall_stats.num_sols += P_stat.num_sols
        self.overall_stats = overall_stats

    def pretty_print(self, sort: str, verbose: int) -> None:
        # print overall stats
        if self.overall_stats:
            self.overall_stats.pretty_print()
        else:
            print(">> overall stats have not been computed")
        if self.solutions:
            print("Best overall square:")
            best_dstat = max(self.diagonal_stats, key=lambda x: (sorted(x.score_counts.items(), reverse=True), -x.P_val, -x.S, -x.max_nb))
            best_dstat.pretty_print_best_square()

        # print P stats
        print("Statistics for each P, sorted by factorization of P:")
        sorted_P_stats_by_P = sorted(self.P_stats.values(), key=lambda x: (len(x.P), x.P[::-1]))
        PStats.pretty_print_header()
        for P_stat in sorted_P_stats_by_P:
            P_stat.pretty_print()

        if verbose >= 1 and len(self.P_stats) > 1:
            print("Statistics for each P, sorted by count:")
            sorted_P_stats = sorted(self.P_stats.values(), key=lambda x: x.count)
            PStats.pretty_print_header()
            for P_stat in sorted_P_stats:
                P_stat.pretty_print()

            print("Statistics for each P, sorted by P:")
            sorted_P_stats_by_P = sorted(self.P_stats.values(), key=lambda x: x.P_val)
            PStats.pretty_print_header()
            for P_stat in sorted_P_stats_by_P:
                P_stat.pretty_print()

        # print solutions and diagonal stats
        if self.solutions:
            sorted_dstats_by_score = sorted(
                self.diagonal_stats, 
                key=lambda x: (sorted(x.score_counts.items(), reverse=True), -x.P_val, -x.S, -x.max_nb),
                reverse=True
            )
            sorted_dstats_by_P = sorted(
                self.diagonal_stats, 
                key=lambda x: (x.P, x.S, x.max_nb)
            )
            if verbose == 0:
                print("Diagonal statistics for top 100 solutions, sorted by score:")
                DiagonalStats.pretty_print_header()
                for diagonal_stat in sorted_dstats_by_score[:100]:
                    diagonal_stat.pretty_print()
            if (verbose >= 1 and sort == "score") or (verbose >= 2):
                print("Diagonal statistics for each solution, sorted by score:")
                DiagonalStats.pretty_print_header()
                for diagonal_stat in sorted_dstats_by_score:
                    diagonal_stat.pretty_print()
            if (verbose >= 1 and sort == "P") or (verbose >= 2):
                print("Diagonal statistics for each solution, sorted by P:")
                DiagonalStats.pretty_print_header()
                for diagonal_stat in sorted_dstats_by_P:
                    diagonal_stat.pretty_print()
            if (verbose >= 3 and sort == "score") or (verbose >= 4):
                print("Best diagonals for each solution, sorted by score:")
                for diagonal_stat in sorted_dstats_by_score:
                    diagonal_stat.pretty_print_best_square()
            if (verbose >= 3 and sort == "P") or (verbose >= 4):
                print("Best diagonals for each solution, sorted by P:")
                for diagonal_stat in sorted_dstats_by_P:
                    diagonal_stat.pretty_print_best_square()

        # print S stats (very long!!)
        if verbose >= 5:
            print("Detailed statistics for each P:")
            for P_stat in sorted_P_stats_by_P:
                PStats.pretty_print_header()
                P_stat.pretty_print()
                P_stat.pretty_print_SStats(verbose)

# given the permutation for the first diagonal, returns possible permutations for the second diagonal
def other_diagonals(permutation: Vec) -> list[Vec]:
    out = []
    for grouping in groups:
        gmap = {g[0]: g[1] for g in grouping}
        gmap.update({g[1]: g[0] for g in grouping})
        out += [tuple(permutation[gmap[i]] for i in range(6))]
    return out

def parse_arrangement_output(
    file: str, 
    stats: CSearchStats = CSearchStats(None, [], [], {})
) -> CSearchStats:
    cur_sum = 0
    P: Union[tuple[int], None] = None    
    P_val: Union[int, None] = None   
    try:
        P = tuple()
        for x in file.split('.')[-2].split('/')[-1].split('_')[1:]:
            if x.isdigit():
                P += (int(x),)
            else:
                break
        #  print("inferred value of P:", P)
        P_val = reduce(mul, (p ** r for p, r in zip(primes, P)), 1)
        #  print("P_val:", P_val)
    except:
        pass
    P_stat: PStats = stats.P_stats.get(P, default=PStats(P, P_val, []))
    with open(file, "r") as f:
        while line := f.readline():
            split = line.split()
            if split == ['solution', 'found']:
                table = Table([], [], [], [])
                # parse rows
                for _ in range(6):
                    ri, *row = (int(n) for n in f.readline().split())
                    table.rows += [tuple(row)]
                    table.ris += [ri]
                # blank line
                f.readline()
                # parse cols
                for _ in range(6):
                    ci, *col = (int(n) for n in f.readline().split())
                    table.cols += [tuple(col)]
                    table.cis += [ci]
                # compute diagonal stats
                solution = Solution(table)
                diagonal_stats = solution.diagonal_stats(P, P_stat.count)
                # record solution and diagonal stats
                stats.solutions += [solution]
                stats.diagonal_stats += [diagonal_stats]
                # update P stats
                P_stat.update(diagonal_stats)
            elif split[0] == "sum":
                cur_sum = int(split[1])
                num_vecs = int(split[3])
                # update P stats
                P_stat.num_vecs += num_vecs
                P_stat.num_S += 1
                P_stat.max_S = max(cur_sum, P_stat.max_S)
                P_stat.max_num_vecs = max(num_vecs, P_stat.max_num_vecs)
                P_stat.S_stats += [SStats(cur_sum, num_vecs, 0, 0, -1)]
            elif split[:2] == ["num", "searched:"]:
                count = int(split[2])
                P_stat.count += count
                P_stat.S_stats[-1].count = count
            elif split[:2] == ["completed", "in"]:
                P_stat.time += float(split[2])
    stats.P_stats[P] = P_stat
    return stats

def parse_directory(glob_str: str) -> CSearchStats:
    stats = CSearchStats(None, [], [], {})
    print("Looking for files matching pattern", glob_str, file=stderr)
    for file in glob(glob_str):
        print("Processing", file, file=stderr)
        stats = parse_arrangement_output(file, stats)
    return stats

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--file', type=str)
    group.add_argument('--glob', type=str)
    parser.add_argument("--verbose", metavar="V", type=int, default=1)
    parser.add_argument("--expect-nsols", type=int, nargs=1)
    parser.add_argument("--sort", choices=['P', 'score'], nargs=1)
    args = parser.parse_args()

    sort = args.sort if args.sort else ("score" if args.glob else "P") 
    stats = parse_directory(args.glob) if args.glob else parse_arrangement_output(args.file) 
    if args.expect_nsols:
        stats.assert_num_sols(args.expect_nsols[0])
    stats.compute_overall_stats()
    stats.pretty_print(sort, args.verbose)




