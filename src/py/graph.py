import sys

from enumeration import gen_rows, rows_to_rowdict
from arrangement import make_vecs

'''
Code for generating LAD-formatted graphs
Example uses:
    `python graph.py pattern`
    `python graph.py target 327 10 4 3 2`
Output graphs can be used with Glasgow subgraph solver like this:
    `(path-to-glasgow-subgraph-solver) --format lad 6x6.lad_pattern 6x6_P_10_4_3_2_S_327.lad_target
'''

def make_pattern(N: int) -> None:
    with open(f'{N}x{N}.lad_pattern', 'w') as f:
        print(N*N + 2*N + 3, file=f)
        for i in range(N*N):
            print(f'2 {N*N + i // N} {N*N + N + i % N}', file=f)
        for i in range(N):
            print(' '.join(map(str, [N+1] + list(range(N*i, N*(i+1))) + [N*N + 2*N + 2])), file=f)
        for i in range(N):
            print(' '.join(map(str, [N+1] + list(range(i, N*N+i, N)) + [N*N + 2*N + 2])), file=f)
        # make dummy triangle to force correct direction of bipartite graph
        print(f'2 {N*N + 2*N + 1} {N*N + 2*N + 2}', file=f)
        print(f'2 {N*N + 2*N} {N*N + 2*N + 2}', file=f)
        print(' '.join(map(str, [2*N+2] + list(range(N*N, N*N+2*N+2)))), file=f)
        

def make_target(S: int, P: tuple[int, ...], N: int) -> None:
    all_rows = gen_rows(P, N)
    row_dict = rows_to_rowdict(all_rows)

    filt_rows = row_dict[S]
    inv_map, vecs, n_vecs = make_vecs(S, filt_rows, start=0)

    # print(vecs)

    n_nums = max(max(n for n in vec) for vec in vecs) + 1
    # vertices 0,...,n-1 correspond to numbers
    # vertices n,... correspond to vecs
    # print(n_nums)

    num_to_vec = {n: [] for n in range(n_nums)}
    for i, vec in enumerate(vecs):
        for n in vec:
            num_to_vec[n] += [i + n_nums]
    n_vtxs = n_nums + n_vecs

    fname = f"{N}x{N}_P_" + ('_'.join(map(str, P))) + f"_S_{S}.lad_target"
    with open(fname, 'w') as f:
        print(n_vtxs + 3, file=f)
        for n in range(n_nums):
            print(' '.join(map(str, [len(num_to_vec[n])] + num_to_vec[n])), file=f)
        for vec in vecs:
            print(' '.join(map(str, [N+1] + list(vec) + [n_vtxs + 2])), file=f)
        # make dummy triangle to force correct direction of bipartite graph
        print(f'2 {n_vtxs + 1} {n_vtxs + 2}', file=f)
        print(f'2 {n_vtxs} {n_vtxs + 2}', file=f)
        print(' '.join(map(str, [n_vecs + 2] + list(range(n_nums, n_vtxs + 2)))), file=f)

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "pattern":
        make_pattern(6)
    if mode == "target":
        S = int(sys.argv[2])
        P = tuple([int(i) for i in sys.argv[3:]])
        N = 6
        make_target(S, P, N)

        


    


