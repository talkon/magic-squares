import dataclasses
import argparse
from arrangement import Table

@dataclasses.dataclass
class CSearchStats:
    count: int
    time: float
    solutions: list[Table]

def parse_arrangement_output(file: str) -> CSearchStats:
    stats = CSearchStats(0, 0, [])
    with open(file, "r") as f:
        while line := f.readline():
            split = line.split()
            if split == ['solution', 'found']:
                solution = Table([], [], [], [])
                # parse rows
                for _ in range(6):
                    ri, *row = (int(n) for n in f.readline().split())
                    solution.rows += [tuple(row)]
                    solution.ris += [ri]
                # blank line
                f.readline()
                # parse cols
                for _ in range(6):
                    ci, *col = (int(n) for n in f.readline().split())
                    solution.cols += [tuple(col)]
                    solution.cis += [ci]
                stats.solutions += [solution]
            elif split[:2] == ["num", "searched:"]:
                stats.count = int(split[2])
            elif split[:2] == ["completed", "in"]:
                stats.time = float(split[2])
    return stats

#TODO: verify correctness of found solutions
#TODO: check if solutions have diagonals

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=str)
    args = parser.parse_args()

    parsed = parse_arrangement_output(args.file)
    print(parsed)




