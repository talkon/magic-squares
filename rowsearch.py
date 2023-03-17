import time
import random
import sys
from multiprocessing import Pool

def process(input):

  S, filt_rows = input
  start = time.time()

  # print(f"[{S}] (Starting)", flush=True)

  # filt_rows = row_dict[S] 
  # sorted([row for row in all_rows if sum(row) == S and len(set(row)) == len(row)], reverse=True)
  l1 = len(filt_rows)
  # print(filt_rows)

  eltc = {}
  eltcl = []

  while True:
    if len(filt_rows) <= 11:
      # print(f"[{S}] (Eliminated by reduction)")
      print(f"[{S}] (Completed in {time.time()-start} seconds by reduction)", flush=True)
      return 0, [], []
    eltc = {}
    for row in filt_rows:
      for elt in row:
        eltc[elt] = eltc.get(elt, 0) + 1
    eltcl = sorted(eltc.items(), key=lambda x: x[1], reverse=True)
    if eltcl[-1][1] >= 2:
      break
    else:
      filt_rows = [row for row in filt_rows if all(eltc[elt] > 1 for elt in row)]
  
  # print(f"[{S}] (Setup) reduced from {l1} vecs to {len(filt_rows)} vecs", flush=True)

  rename_map = {}
  inv_map = {}
  for i, elt in enumerate(eltcl):
    rename_map[elt[0]] = i
    inv_map[i] = elt[0]

  vecs = sorted([tuple(sorted([rename_map[i] for i in row], reverse=True)) for row in filt_rows], reverse=True)
  # for vec in vecs:
  #   print(vec)
  nvecs = len(vecs)

  # print(">>> finished in", time.time()-start, "seconds")

  def make_intersections(vecs):
    d = {}
    for i in range(len(vecs)):
      dd = {}
      for j in range(len(vecs)):
          dd[j] = len(set(vecs[i]).intersection(set(vecs[j])))
      d[i] = dd
    return d

  # print(">>> Precomputing intersections")
  # start = time.time()
  intersections = make_intersections(vecs)
  # print(">>> finished in", time.time()-start, "seconds")

  count = [0]
  search_start = time.time()
  sols = []
  nms = []

  def record(rows, cols, ris, cis, count, sols, nms):

    count[0] += 1 
    if count[0] % 100000 == 0:
      rs = [tuple(inv_map[e] for e in row) for row in rows]
      cs = [tuple(inv_map[e] for e in col) for col in cols]
      print(f"[{S}] (Log) {count[0]} {time.time()-start} {len(vecs)} {sorted(ris+cis)}", flush=True)

    score = len(rows) * len(cols)

    if score == 25:
      rs = [tuple(inv_map[e] for e in row) for row in rows]
      cs = [tuple(inv_map[e] for e in col) for col in cols]
      missing_elt = 6 * S - sum(set().union(*rs).union(*cs))
      last_row = list(set().union(*rs).difference(set().union(*cs))) + [missing_elt]
      pl = 1
      p1 = 1
      for e1, el in zip(rs[0], last_row):
        pl = pl * el
        p1 = p1 * e1
      # print(e, p1, pl, flush=True)
      # if p1 == pl:
      #   print(f"[{S}] (SOLUTION  6x6) {(rs, cs)}", flush=True)
      #   sols += [(rs, cs)]
      if p1 != pl:
        print(f"[{S}] (Near miss 5x5) {(rs, cs)}", flush=True)
        nms += [(rs, cs)]
    
    if score == 36:
      rs = [tuple(inv_map[e] for e in row) for row in rows]
      cs = [tuple(inv_map[e] for e in col) for col in cols]
      print(f"[{S}] (SOLUTION  6x6) {(rs, cs)}", flush=True)
      sols += [(rs, cs)]

  def search_aux(i, rows, cols, ris, cis, unmatched):
    # print(i, rows, cols)
    if (i < nvecs):
      mr = vecs[i][0]
      if unmatched and max(unmatched) > mr:
        return
      # union_rows = set().union(*rows)
      # union_cols = set().union(*cols)
      # if not union_cols.issuperset(sorted(e for e in union_rows if e > mr)):
      #   return
      # if not union_rows.issuperset(sorted(e for e in union_cols if e > mr)):
      #   return
    
    record(rows, cols, ris, cis, count, sols, nms)

    for j in range(i, nvecs):
      vec = vecs[j]

      is_valid_row = all(intersections[j][ri] == 0 for ri in ris) and all(intersections[j][ci] == 1 for ci in cis) 
      if is_valid_row:
        um = unmatched.symmetric_difference(vec)
        search_aux(j+1, rows + [vec], cols, ris + [j], cis, um)
      
      is_valid_col = len(rows) > 0 and all(intersections[j][ri] == 1 for ri in ris) and all(intersections[j][ci] == 0 for ci in cis) 
      if is_valid_col:
        um = unmatched.symmetric_difference(vec)
        search_aux(j+1, rows, cols + [vec], ris, cis + [j], um)

  # print(f"[{S}] (Searching)", flush=True)
  search_aux(0, [], [], [], [], set())
  # print(f"[{S}] (Log) {count[0]} {time.time()-search_start} {len(vecs)}")
  print(f"[{S}] (Completed in {time.time()-start} seconds: explored {count[0]} configs, {len(nms)} near misses, {len(sols)} solutions)", flush=True)

  return count[0], sols, nms

if __name__ == '__main__':
  tup_bank = {}
  primes = (2, 3, 5, 7, 11, 13, 17, 19, 23)

  P = tuple([int(i) for i in sys.argv[1:]])
  N = 6

  def gen_tup(s, n):
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

  def gen_rows(pows, n):
    rows = {(1,) * n}
    for p, pow in zip(primes, pows):
      ptups = gen_tup(pow, n)
      new_rows = set()
      for row in rows:
        for tup in ptups:
          new_rows.add(tuple(
            sorted([elt * (p ** r) for elt, r in zip(row, tup)], reverse=True)
          ))
      rows = new_rows
    return rows

  print(">>> Generating and simplifying rows")
  start = time.time()

  all_rows = gen_rows(P, N)

  counts = {}
  unique_counts = {}
  for row in all_rows:
    s = sum(row)
    counts[s] = counts.get(s, 0) + 1
    if len(row) == len(set(row)):
      unique_counts[s] = unique_counts.get(s, 0) + 1

  ucs = sorted(unique_counts.items(), key=lambda x: x[1], reverse=True)
  S_list = [S for S, c in ucs if c >= 12]
  S_set = set(S_list)
  row_dict = {S: [] for S in S_list}

  print(f">>> total possible 6-vecs: {sum(u[1] for u in ucs)}")

  for row in all_rows:
    if len(set(row)) == len(row) and sum(row) in S_set:
      row_dict[sum(row)] += [row]

  # print(S_list)
  print(f">>> number of S to consider: {len(S_list)}")

  print(f">>> finished in {time.time()-start} seconds", flush=True)

  print(">>> Starting parallel search", flush=True)
  l = list(row_dict.items())
  random.shuffle(l)
  search_start = time.time()
  with Pool(50) as p:
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
    print(f">>> finished search for P = {P} in {time.time()-search_start} seconds", flush=True)
    print(f">>> total time: {time.time()-start} seconds", flush=True)
    print(f">>> total configurations: {count}", flush=True)
    print(f">>> solutions found: {len(sols)}", flush=True)
    print(f">>> solutions: {sols}", flush=True)
    print(f">>> near misses found: {len(nms)}", flush=True)
    print(f">>> near misses: {nms}", flush=True)
    
