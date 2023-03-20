import dataclasses
import math
import time
from collections import Counter
from enum import Enum
from typing import Callable


Vec = tuple[int, int, int, int, int, int]


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
    Remove vecs that can't appear in the table, and relabel by frequency.
    Returns label_to_elt dictionary and the new vecs.
    """
    count = Counter()

    while len(vecs) >= 12:
        count = Counter([elt for vec in vecs for elt in vec])
        if count.most_common()[-1][1] >= 2:
            break
        vecs = [vec for vec in vecs if all(count[elt] > 1 for elt in vec)]
    else:
        return {}, []  # not enough vecs to fill grid

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

        if math.prod(elt_table.rows[0]) == math.prod(last_row):
            return Judgment.SOLUTION
        else:
            return Judgment.NEAR_MISS

    if score > 25:
        raise Exception(f"{score} should never happen?")

    return Judgment.NOTHING


def recorder(
    stats: SearchStats, S: int, label_to_elt: dict[int, int], vecs: list[Vec]
) -> Callable[[Table], None]:
    """
    Returns a function that records a found table, then updates stats.
    """

    def record(table: Table) -> None:
        """
        Record a found table, updating stats if needed.
        """
        stats.count += 1
        if stats.count % 100000 == 0:
            nvecs = len(vecs)
            indices = sorted(table.ris + table.cis)
            log("Log", f"{stats.count} {elapsed(stats)} {nvecs} {indices}")

        elt_table = to_elts(label_to_elt, table)
        match judge(S, elt_table):
            case Judgment.NOTHING:
                pass
            case Judgment.NEAR_MISS:
                log("near miss 5x5", f"{elt_table}")
                stats.near_misses.append(elt_table)
            case Judgment.SOLUTION:
                log("SOLUTION  6x6", f"{elt_table}")
                stats.solutions.append(elt_table)

    return record


def searcher(
    vecs: list[Vec],
    intersections: dict[int, dict[int, int]],
    record: Callable[[Table], None],
) -> Callable[[int, Table, set[int]], None]:
    """
    Returns a backtracking search function.
    """

    def search_aux(i: int, table: Table, unmatched: set[int]):
        """
        Backtracking: start searching from partial table, with vecs beginning
        from index i.
        """
        if i < len(vecs):
            max_elt = vecs[i][0]
            if unmatched and max(unmatched) > max_elt:
                return

        record(table)

        for j in range(i, len(vecs)):
            vec = vecs[j]

            is_valid_row = (
                # can't overlap existing rows:
                all(intersections[j][ri] == 0 for ri in table.ris)
                # must overlap existing cols:
                and all(intersections[j][ci] == 1 for ci in table.cis)
            )
            if is_valid_row:
                new_table = dataclasses.replace(
                    table, rows=table.rows + [vec], ris=table.ris + [j]
                )
                search_aux(j + 1, new_table, unmatched ^ set(vec))

            is_valid_col = (
                len(table.rows) > 0
                # must overlap existing rows:
                and all(intersections[j][ri] == 1 for ri in table.ris)
                # can't overlap existing cols:
                and all(intersections[j][ci] == 0 for ci in table.cis)
            )
            if is_valid_col:
                new_table = dataclasses.replace(
                    table, cols=table.cols + [vec], cis=table.cis + [j]
                )
                search_aux(j + 1, new_table, unmatched ^ set(vec))

    return search_aux


def process(input: tuple[int, list[Vec]]) -> tuple[int, list[Table], list[Table]]:
    global PROCESS_ID
    S, vecs = input
    PROCESS_ID = S
    stats = SearchStats(0, time.time(), [], [])

    label_to_elt, vecs = make_vecs(vecs)
    if not vecs:
        log(f"Completed in {elapsed(stats)} seconds", "by reduction")
        return 0, [], []

    intersections = make_intersections(vecs)

    record = recorder(stats, S, label_to_elt, vecs)
    search_aux = searcher(vecs, intersections, record)

    # log("Searching", "")
    search_aux(0, Table([], []), set())
    log(
        f"Completed in {elapsed(stats)} seconds",
        f"explored {stats.count} configs, {len(stats.near_misses)} near misses, {len(stats.solutions)} solutions",
    )

    return stats.count, stats.solutions, stats.near_misses
