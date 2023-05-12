import dataclasses
import math
import time
from collections import Counter
from enum import Enum
from typing import Callable
import bisect


Vec = tuple[int]


class Judgment(Enum):
    NOTHING = 0
    NEAR_MISS = 1
    SOLUTION = 2


@dataclasses.dataclass
class Table:
    rows: list[Vec]
    cols: list[Vec]
    ris: list[int] = dataclasses.field(default_factory=list)  # row indices
    cis: list[int] = dataclasses.field(default_factory=list)  # column indices


@dataclasses.dataclass
class SearchStats:
    count: int
    start_time: float
    near_misses: list[Table]
    solutions: list[Table]


def elapsed(stats: SearchStats) -> float:
    """
    Return elapsed time since search start.
    """
    return time.time() - stats.start_time


PROCESS_ID = None


def log(type: str, message: str) -> None:
    print(f"[{PROCESS_ID}] ({type}) {message}", flush=True)


def make_vecs(vecs: list[Vec]) -> tuple[dict[int, int], list[Vec]]:
    """
    Relabel vecs by frequency. Returns label_to_elt dictionary and the new vecs.
    """
    count = Counter([elt for vec in vecs for elt in vec])

    # log(f"(Setup) reduced from {l1} vecs to {len(sp_sets)} vecs")

    elt_to_label = {}
    for i, (elt, _) in enumerate(count.most_common()):
        elt_to_label[elt] = i  # more common elt = smaller label

    vecs = sorted(
        [tuple(sorted([elt_to_label[i] for i in vec], reverse=True)) for vec in vecs],
        reverse=True,
    )

    label_to_elt = {label: elt for elt, label in elt_to_label.items()}
    return label_to_elt, vecs


def make_intersections(vecs: list[Vec]) -> dict[int, dict[int, int]]:
    """
    Precompute intersections[i][j] = len(set(vecs[i]) & set(vecs[j])).
    """
    intersections = {}
    for i, vec1 in enumerate(vecs):
        intersections[i] = {
            j: len(set(vec1) & set(vec2)) for j, vec2 in enumerate(vecs)
        }
    return intersections


def to_elts(label_to_elt: dict[int, int], table: Table) -> Table:
    """
    Translate a table of labels to a table of elements.
    """
    return Table(
        [tuple(label_to_elt[label] for label in row) for row in table.rows],
        [tuple(label_to_elt[label] for label in col) for col in table.cols],
    )


def judge(S: int, elt_table: Table) -> Judgment:
    """
    Judge a table of elements, and returns if it's a solution, a near miss, or
    neither.
    """
    score = len(elt_table.rows) * len(elt_table.cols)

    if score == 25:
        row_elts = {elt for row in elt_table.rows for elt in row}
        col_elts = {elt for col in elt_table.cols for elt in col}
        # missing elt, additively
        missing_elt = 6 * S - sum(row_elts | col_elts)
        last_row = [*(row_elts - col_elts), missing_elt]

        # this is slightly silly code, but we need to ensure solutions aren't counted twice
        if math.prod(elt_table.rows[0]) != math.prod(last_row):
            return Judgment.NEAR_MISS
    elif score == 36:
        return Judgment.SOLUTION
    # elif score > 25:
    #     raise Exception(f"{score} should never happen?")

    return Judgment.NOTHING


def recorder(
    stats: SearchStats, S: int, label_to_elt: dict[int, int], vecs: list[Vec]
) -> Callable[[Table], Judgment]:
    """
    Returns a function that records a found table, then updates stats.
    """

    def record(table: Table) -> Judgment:
        """
        Record a found table, updating stats if needed.
        """
        stats.count += 1
        if stats.count % 100000 == 0:
            nvecs = len(vecs)
            indices = sorted(table.ris + table.cis)
            log("Log", f"{stats.count} {elapsed(stats)} {nvecs} {indices}")

        elt_table = to_elts(label_to_elt, table)
        judgment = judge(S, elt_table)
        if judgment == Judgment.NEAR_MISS:
            log("near miss 5x5", f"{elt_table}")
            stats.near_misses.append(elt_table)
        elif judgment == Judgment.SOLUTION:
            log("SOLUTION  6x6", f"{elt_table}")
            stats.solutions.append(elt_table)
        return judgment

    return record


def searcher(
    vecs: list[Vec],
    intersections: dict[int, dict[int, int]],
    record: Callable[[Table], Judgment],
) -> Callable[[int, Table, set[int], tuple[int], tuple[int]], None]:
    """
    Returns a backtracking search function.
    """

    def search_aux(i: int, table: Table, unmatched: set[int], 
                   row_indices: tuple[int], col_indices: tuple[int]):
        """
        Backtracking: start searching from partial table, 
        with possible vec indices in row_indices and col_indices
        """
        if i < len(vecs):
            max_elt = vecs[i][0]
            if unmatched and max(unmatched) > max_elt:
                return

        

        if record(table) != Judgment.NOTHING:
            return

        new_row_indices = []
        for j in row_indices:
            vec = vecs[j]
            is_valid_row = (
                # can't overlap existing rows:
                (len(table.ris) == 0 or intersections[j][table.ris[-1]] == 0)
                # must overlap existing cols:
                and (len(table.cis) == 0 or intersections[j][table.cis[-1]] == 1) 
            )
            if(is_valid_row):
                new_row_indices.append(j)

        new_col_indices = []
        for j in col_indices:
            vec = vecs[j]
            is_valid_col = (
                len(table.rows) > 0
                # must overlap existing rows:
                and (len(table.ris) == 0 or intersections[j][table.ris[-1]] == 1)
                # can't overlap existing cols:
                and (len(table.cis) == 0 or intersections[j][table.cis[-1]] == 0)
            )
            if(is_valid_col):
                new_col_indices.append(j)

        if i==0:
            new_col_indices = new_row_indices
        new_row_indices = tuple(new_row_indices)
        new_col_indices = tuple(new_col_indices)


        temp_row_indices = new_row_indices
        temp_col_indices = new_col_indices

        max_unmatched = max(unmatched, default = 0)

        row_buckets = {}
        col_buckets = {}

        # group rows and cols by max value
        for j in new_row_indices:
            vec = vecs[j]
            if vec[0] < max_unmatched:
                continue
            if vec[0] not in row_buckets:
                row_buckets[vec[0]] = set()
            if vec[0] not in col_buckets:
                col_buckets[vec[0]] = set()
            row_buckets[vec[0]].add(j)

        for j in new_col_indices:
            vec = vecs[j]
            if vec[0] < max_unmatched:
                continue
            if vec[0] not in row_buckets:
                row_buckets[vec[0]] = set()
            if vec[0] not in col_buckets:
                col_buckets[vec[0]] = set()
            col_buckets[vec[0]].add(j)

        for max_el in row_buckets:
            if max_el == max_unmatched:
                for j in row_buckets[max_el]:
                    vec = vecs[j]
                    temp_row_indices = new_row_indices[bisect.bisect_right(new_row_indices, j):]
                    temp_col_indices = new_col_indices[bisect.bisect_right(new_col_indices, j):]
                    new_table = dataclasses.replace(
                            table, rows=table.rows + [vec], ris=table.ris + [j]
                        )
                    search_aux(j + 1, new_table, unmatched ^ set(vec), temp_row_indices, temp_col_indices)
            else:
                # first set added must be followed by another set added with same max value, 
                # if its max is greater than max_unmatched
                for j1 in row_buckets[max_el]:
                    for j2 in col_buckets[max_el]:
                        if i == 0 and j1 > j2:
                            continue
                        if intersections[j1][j2] == 1:
                            vec1 = vecs[j1]
                            vec2 = vecs[j2]
                            j = max(j1, j2)
                            temp_row_indices = new_row_indices[bisect.bisect_right(new_row_indices, j):]
                            temp_col_indices = new_col_indices[bisect.bisect_right(new_col_indices, j):]
                            new_table = dataclasses.replace(
                                    table, rows=table.rows + [vec1], ris=table.ris + [j1],
                                        cols=table.cols + [vec2], cis=table.cis + [j2]
                                )
                            search_aux(j + 1, new_table, unmatched ^ set(vec1) ^ set(vec2), temp_row_indices, temp_col_indices)

        if max_unmatched in col_buckets:
            for j in col_buckets[max_unmatched]:
                vec = vecs[j]
                temp_row_indices = new_row_indices[bisect.bisect_right(new_row_indices, j):]
                temp_col_indices = new_col_indices[bisect.bisect_right(new_col_indices, j):]
                new_table = dataclasses.replace(
                        table, cols=table.cols + [vec], cis=table.cis + [j]
                    )
                search_aux(j + 1, new_table, unmatched ^ set(vec), temp_row_indices, temp_col_indices)

    return search_aux


def process(S: int, vecs: list[Vec]) -> SearchStats:
    global PROCESS_ID
    PROCESS_ID = S
    stats = SearchStats(0, time.time(), [], [])

    label_to_elt, vecs = make_vecs(vecs)
    intersections = make_intersections(vecs)

    record = recorder(stats, S, label_to_elt, vecs)
    search_aux = searcher(vecs, intersections, record)

    # log("Searching", "")
    search_aux(0, Table([], []), set(), tuple(i for i in range(len(vecs))), tuple(i for i in range(len(vecs))))
    log(
        f"Completed in {elapsed(stats)} seconds",
        f"explored {stats.count} configs, {len(stats.near_misses)} near misses, {len(stats.solutions)} solutions",
    )

    return stats
