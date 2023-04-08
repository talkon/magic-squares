import dataclasses
import argparse
from arrangement import Table, Vec
from operator import mul
from functools import reduce
from itertools import permutations
from enum import IntEnum

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
    P: int
    S_count: int
    P_count: int
    SP_count: int
    score_counts: dict[int, int]
    best_score: int
    best_square: tuple[tuple[int]]

    def pretty_print(self, verbose: int) -> None:
        print(f"{self.P:13} {self.S:4} {self.S_count:3} {self.P_count:3} {self.SP_count:3} {self.best_score:4}  {sorted(self.score_counts.items())}")
        if verbose >= 2:
            for row in self.best_square:
                print('    ' + ' '.join([f"{i:4}" for i in row]))
            

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

    def diagonal_stats(self) -> DiagonalStats:
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

        return DiagonalStats(self.S, self.P, S_count, P_count, SP_count, score_counts, best_score, best_square)

@dataclasses.dataclass
class CSearchStats:
    count: int
    time: float
    nsums: int
    nvecs: int
    max_nvecs: int
    nsols: int
    solutions: list[Solution]
    diagonal_stats: list[DiagonalStats]

    def pretty_print(self, verbose: int) -> None:
        print("Overall statistics:")
        print(f"  {'count':<10}={self.count:12}")
        print(f"  {'time':<10}={self.time:12}")
        print(f"  {'nsums':<10}={self.nsums:12}")
        print(f"  {'nvecs':<10}={self.nvecs:12}")
        print(f"  {'max_nvecs':<10}={self.max_nvecs:12}")
        print(f"  {'nsols':<10}={self.nsols:12}")
        if self.nsols:
            if verbose >= 1:
                print("Diagonal statistics for each solution:")
                print(f"{'P':>13} {'S':>4} {'#S':>3} {'#P':>3} {'#SP':>3} {'best':>4}  {'score_counts'}")
                for diagonal_stat in self.diagonal_stats:
                    diagonal_stat.pretty_print(verbose)
            best_dstat = max(self.diagonal_stats, key=lambda x: x.best_score)
            print("Best found square:")
            print(f"  {'P':<7}={best_dstat.P:13}")
            print(f"  {'S':<7}={best_dstat.S:13}")
            print(f"  {'score':<7}={best_dstat.best_score:13}")
            print(f"  {'square':<7}=")
            for row in best_dstat.best_square:
                print('    ' + ' '.join([f"{i:4}" for i in row]))



# given the permutation for the first diagonal, returns possible permutations for the second diagonal
def other_diagonals(permutation: Vec) -> list[Vec]:
    out = []
    for grouping in groups:
        gmap = {g[0]: g[1] for g in grouping}
        gmap.update({g[1]: g[0] for g in grouping})
        out += [tuple(permutation[gmap[i]] for i in range(6))]
    return out

def parse_arrangement_output(file: str) -> CSearchStats:
    stats = CSearchStats(0, 0, 0, 0, 0, 0, [], [])
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
                solution = Solution(table)
                stats.solutions += [solution]
                stats.diagonal_stats += [solution.diagonal_stats()]
                stats.nsols += 1
            elif split[0] == "sum":
                nvecs = int(split[3])
                stats.nvecs += nvecs
                stats.nsums += 1
                if nvecs > stats.max_nvecs:
                    stats.max_nvecs = nvecs
            elif split[:2] == ["num", "searched:"]:
                stats.count += int(split[2])
            elif split[:2] == ["completed", "in"]:
                stats.time += float(split[2])
    return stats

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str)
    parser.add_argument("--verbose", "-v", action="count", default=0)
    args = parser.parse_args()

    parsed = parse_arrangement_output(args.file)
    parsed.pretty_print(args.verbose)



