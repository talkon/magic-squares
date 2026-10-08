# Additive-multiplicative magic squares: literature and calibration data (October 2026)

Lens: literature and calibration. This note collects what is known about
additive-multiplicative ("add-mult") magic and semi-magic squares of each order,
published counting results that bear on a heuristic, and constructions. It then
runs three cheap calibration checks against exact counts and our own squares.

**How the facts were obtained, and how far to trust them.** Direct page fetches
were blocked by the egress proxy for every relevant host (multimagie.com, arxiv.org,
oeis.org, mathworld, wikipedia, combinatorics.org, umn.edu). The facts therefore
come from four places, each tagged:

* **[V]** verified here: the square was re-checked arithmetically (all line sums
  and products, distinct entries), or the number was recomputed (factorizations,
  vector counts with `bin/enumerate --vec-size n`).
* **[S]** web-search summaries of the named page. Usually the numbers are quoted
  literally. Two numbers were garbled and are flagged where they appear.
* **[R]** the project's report (`project-report.pdf`, 2023) or the README.
* **[M]** my background knowledge, mostly bibliographic details. These are
  plausible but unchecked here.

Wherever a quoted square was available, its numbers were confirmed [V].

---

## 1. Status by order

S = magic sum, P = magic product, "Max" = largest entry. "Semi-magic" means
rows and columns only. S_min(P) is the smallest sum of n distinct divisors of P
with product P. N is the number of such n-sets with sum S, computed here with
`bin/enumerate --vec-size n` [V].

| n | add-mult magic square | smallest known semi-magic | notes |
|---|---|---|---|
| 2 | impossible (trivial) | impossible | |
| 3 | **impossible** (elementary, see §3) | not checked here. The §3 argument needs the diagonals, so it does not cover semi-magic squares | |
| 4 | **impossible**: L. Morgenstern 2007 proof [S]; also in Cilleruelo–Luca 2010 [S] | **exists**. Morgenstern 2007, 54 examples with Max < 256 [S]. Smallest S = 247, P = 3,369,600 = 2^7 3^4 5^2 13, Max 156 [S][V]. Smallest Max = 110, S = 325, P = 39,916,800 = 11! [S][V] | S=247: S_min 172 (S/S_min 1.44), **N = 8**: the 4 rows + 4 cols are the only 4-sets [V] |
| 5 | **open**. Prize: Boyer's Main Enigma 6, €1000 + champagne [S] | **exists**. Morgenstern 2007, "20 semi-magic examples with Max nb < 276" [S]. Smallest Max 182, S = 476, P = 3,424,861,440 = 2^8 3^5 5 7 11^2 13 [S][V]; it even has one diagonal with sum 476 (S-type) [V] | S_min = 404 (1.18), N = 78 [V] |
| 6 | **open**. Prize: Small Enigma 6a, €500 + champagne [S] | **exists**. Morgenstern 2007 (first ones) [S]. Smallest S = **289** (Max 135, P = 1,596,672,000 = 2^11 3^4 5^3 7 11; S+P type). Smallest Max = **105** (S = 316, P = 7,264,857,600 = 2^9 3^4 5^2 7^2 11 13). Smallest P = **508,032,000** = 2^10 3^4 5^3 7^2 (S = 327, Max 200) [S][R][V] | N(S) = 483, 516, 451 for these three; S/S_min = 1.40, 1.19, 1.91 [V] |
| 7 | **exists**: S. Miquel, Aug 2016, S = 465, P = 150,885,504,000 = 2^10 3^7 5^3 7^2 11, Max 252, own Rust program on an i7-920 PC [S], "around 600 hours" [R] | T. Shirakawa, Apr–May 2010 (first 7x7 semi-magic): Max 168, S = 310, P = 14,529,715,200 = 2^10 3^4 5^2 7^2 11 13 (S+S type) [S][R]. Oct 2010: "nearly magic" SP+P, Max 154, S = 380, P = 14! = 87,178,291,200 [S][R] | magic: S_min 279 (1.67), N = 8,963. semi: S_min 200 (1.55), N = 3,834 [V] |
| 8 | **exists**: W. W. Horner 1955, S = 840, P = 2,058,068,231,856,000, Max 261 [S]. Boyer Nov 2005: S = 600, P = 67,463,283,888,000, Max 225; and S = 760, P = 51,407,948,592,000, Max 333 [S]. Shirakawa Feb 2013: Max 117, P = 15! = 1,307,674,368,000 (S not retrieved) [S] | – | the 2013 squares appear on Boyer's add-mult 8th/9th-order page, but one summary called them "multiplicative" |
| 9 | **exists**: Horner 1952, S = 848, P = 5,804,807,833,440,000, Max 290 [S]. Boyer 2005 improved it [S]. Shirakawa Feb 2013: Max 153, P = 266,765,571,072,000 [S] | – | |
| 10, 11, 14, 22, 26 | first ones by Shirakawa, Apr–Jun 2010. 10x10: S = 37,800, P ≈ 1.60e33, Max 9963 [S] | – | |
| 12, 13, 14, 15, 17, 19, 21 | best known by Shirakawa (Aug 2013 news) [S] | – | |
| 16 | Horner 1955 (orders 8 and 16) [S] | – | |
| 18 | Zhu Jiacheng 1996, from four 3x3 magic rectangles, 324 distinct entries [S] | – | J. Statist. Plann. Inference 51(3) 331–337 |

Smallest known add-mult magic square overall: Miquel's 7x7 [S][R][V].

### 1.1 Partial results ("types") for 6x6

The report's notation: S / P / SP = one traversal (a diagonal after permuting rows
and columns) with the magic sum / product / both; "+" joins the two diagonals.

| type | first / smallest | source |
|---|---|---|
| 0 (semi-magic) | Morgenstern 2007; smallest P S=327 | [S][R] |
| S+P | Morgenstern 2007, S = 289, P = 2^11 3^4 5^3 7 11 | [R][S] |
| S+S | Morgenstern 2008, S = 360, P = 2^10 3^6 5^3 7 13 = 8,491,392,000 | [R][V factorization] |
| S, P, P+P, SP | 427,869 / 178,153 / 721 / 519 squares in the 2023 search (758,949 total) | [R] |
| SP+S | 9 squares in 2023 (smallest P: S = 577); 1 in Oct 2026 (P = 2^12 3^6 5^3 7^2 11 17, S = 836) | [R] |
| SP+P | 1 square in 2023 (P = 2^13 3^5 5^3 7^2 13, S = 632); 1 in Oct 2026 (P = 2^16 3^5 5^4 7^2, S = 849, re-checked here) | [R][V] |
| SP+SP (magic) | none known | |

For 7x7 the ladder went S+S, then SP+P (Shirakawa 2010), then SP+SP (Miquel 2016).

### 1.2 Counts of known squares (calibration of "how many exist")

* 4x4 semi-magic: 54 with Max < 256 (Morgenstern) [S].
* 5x5 semi-magic: 20 with Max < 276 (Morgenstern) [S]. A summary also says "there
  is no semi-magic example with Max nb < 91". It is ambiguous whether that refers
  to 5x5 or 7x7.
* 6x6 semi-magic: 758,949 (2023 search, ~10,000 CPU-h) [R]. 5,821 more near
  S_min in the Oct 2026 runs (6,498 distinct in the data used below) [R]. Ours
  include an independent rediscovery of Morgenstern's Max-105 square, S = 316 [V].
* 7x7 magic: at least 1 (Miquel 2016). I could not find how many there are, or
  any later or smaller one.

### 1.3 Multiplicative-only magic squares (lower bounds on P)

Smallest possible magic product, OEIS A114060 [S]. The values are 216 (3x3),
5040 (4x4; Borkovitz & Hwang 1983 proved it minimal [S]), 302,400 (5x5),
25,945,920 (6x6), 3,632,428,800 (7x7; entries from 1..91, Boyer [S]),
670,442,572,800 (8x8) and 140,792,940,288,000 (9x9).

A113026 lists the magic products of 6x6 multiplicative magic squares: 25945920,
26611200, 28828800, ... [S].

So the smallest 6x6 add-mult semi-magic P (5.08e8) is 20x the smallest 6x6
multiplicative magic P. The 7x7 add-mult magic P (1.51e11) is 42x the 7x7
multiplicative minimum.

---

## 2. Prizes and the open problem as posed

Boyer's "Enigmas on magic squares" were announced in 2010 by Plus magazine:
12 prizes totalling €8,000 plus 12 bottles of champagne [S]. The add-mult ones are:

* **Main Enigma 6** (€1000 + 1 bottle): a 5x5 add-mult magic square of distinct
  positive integers, or a proof of impossibility [S].
* **Small Enigma 6a** (€500 + 1 bottle): the same for 6x6 [S].

One summary rendered these as "61,000 and 6,500", a PDF extraction garbling of
"€1,000 and €500".

As of the searches (Oct 2026) there is no announced 5x5 or 6x6 solution and no
impossibility proof [S].

---

## 3. Impossibility proofs

* **3x3** (elementary; my derivation, consistent with [S]). In a 3x3 additive magic
  square every line through the centre c is {c−x, c, c+x}, and the magic sum is 3c.
  In a 3x3 multiplicative magic square the same argument in log space gives
  P = c^3. A line through the centre then has product c(c^2 − x^2) = c^3, so x = 0
  and the entries repeat.
* **4x4**. Lee Morgenstern 2007, "4x4 magic/semi-magic impossibility proof"
  (multimagie.com/English/Morgenstern05.htm) [S]. Linear identities of 4x4
  additive magic squares and their multiplicative (log) analogues give a pair of
  cells with a+b = c+d and ab = cd. Then {a,b} = {c,d} ("duplication lemma"), so
  the entries are not distinct [S].
  Cilleruelo & Luca (2010) also show 4x4 add-mult magic squares cannot have
  distinct entries, and ask "Are there additive-multiplicative magic squares of
  order r = 5 with distinct entries?" [S].
* **5x5, 6x6**: no impossibility result, and no argument ruling out specific
  P, except what Weisenberg (below) proposes.

---

## 4. Published work on the 5x5 / 6x6 problem and on heuristics

* **D. Weisenberg**, "Some Thoughts on the Search for 5×5 and 6×6
  Additive-Multiplicative Magic Squares", arXiv:2311.06326 (10 Nov 2023), and
  Minnesota J. Undergrad. Math. 8(1) (2024) [S]. Abstract: a square "can be
  described by a form determined by the prime factorizations of its entries";
  identifying these forms "might be helpful in finding such a square or ruling
  out specific magic products" [S]. In our terms the forms are the per-prime
  exponent matrices, i.e. lattice points of e_q·(Birkhoff or magic polytope).
  The search summaries show no counting heuristic or probability estimate in it.
* **Cilleruelo & Luca**, "On multiplicative magic squares", Electron. J. Combin.
  17 (2010) #N8 [S]. It gives a lower bound on max − min for an r×r multiplicative
  magic square with distinct entries: for r = 3, x_M − x_m ≥ x_m^{3/4}, which is
  sharp up to 1+o(1) for an infinite family. It gets the true order of magnitude
  for r = 4, plus the 4x4 add-mult impossibility [S]. This is the only rigorous
  "size" result found. It bounds the spread of entries, not the existence.
* **Borkovitz & Hwang**, "Multiplicative magic squares", Discrete Math. 47 (1983)
  1–11 [S]: constructions and the minimal products for 3x3 and 4x4.
* **The project's own report** (Boonsiriseth, Quines, Raphael, May 2023) [R]:
  - traversal rates 1 in 830 (S), 1 in 2600 (P), 1 in 1.05M (SP);
  - 5400 diagonal pairs, giving P(magic | semi-magic) ≈ 5400/(830^2·2600^2) ≈ 1 in
    860 million, so ~1,200 CPU-years at 45 s per square (the abstract says ~400);
  - SP traversals are 2.05x more common than independence predicts;
  - empirically p_S ∝ 1/S and p_P ∝ 1/τ(P).
  The README forecast refines this to E ≈ 0.03–0.045 magic squares in 10 CPU-years.
* **I found no published heuristic** for the number of add-mult magic squares
  or for the expected size of the smallest 5x5 or 6x6 one. The searches covered
  multimagie, arXiv, OEIS, MathWorld and Plus. The only probability-type
  statements are the report's and Boyer's informal "nobody knows".

---

## 5. Counting results that a heuristic can lean on

### 5.1 Exponent matrices ↔ Birkhoff polytope (multiplicative structure)

For P = ∏ q^{e_q}, an n×n multiplicative semi-magic square of positive integers
(not necessarily distinct) with line product P is the same thing as one
nonnegative integer matrix E_q with all line sums e_q for each prime q. So

  #{multiplicative semi-magic squares with product P} = ∏_q H_n(e_q),
  #{multiplicative magic squares with product P}      = ∏_q M_n(e_q),

where H_n(t) counts the lattice points of t·B_n (B_n is the Birkhoff polytope)
and M_n(t) those of the magic-square polytope [M; standard]. Distinctness then
removes a fraction that vanishes as the exponents grow.

* H_n(t) is a polynomial of degree (n−1)^2 (Stanley, Ehrhart theory). Examples:
  H_3(t) = C(t+2,2) + 3·C(t+3,4) (MacMahon), H_n(1) = n!. **Beck & Pixton**,
  "The Ehrhart polynomial of the Birkhoff polytope", Discrete Comput. Geom. 30
  (2003) 623–637, give H_n for n ≤ 9 and vol(B_10) in arXiv:math/0305332 [S].
* **Canfield & McKay**, "The asymptotic volume of the Birkhoff polytope",
  arXiv:0705.2422 (Online J. Anal. Comb. 4, 2009) [S]. The formula [M] is
  vol B_n = exp(−(n−1)^2 ln n + n^2 − (n−½) ln 2π + 1/3 + o(1)).
* **Canfield & McKay**, "Asymptotic enumeration of integer matrices with large
  equal row and column sums", Combinatorica 30 (2010) 655–680,
  arXiv:math/0703600 [S]. M(m,s;n,t) ≈ C(n+s−1,s)^m C(m+t−1,t)^n / C(mn+λmn−1, λmn)
  with λ = s/n = t/m (Good's estimate) [S], times e^{1/2+o(1)} [M]. Also
  Barvinok & Hartigan, arXiv:0709.3810, on contingency tables and transportation
  polytopes [S].
* **Magic-square polytope** (diagonals included). M_n(t) is a quasi-polynomial
  of degree (n−1)^2 − 2 = n^2 − 2n − 1, two less than H_n: the two diagonal
  constraints cost t^{−2}.
  - Beck, Cohen, Cuomo, Gribelyuk, "The number of 'magic' squares, cubes, and
    hypercubes", Amer. Math. Monthly 110 (2003) 707–717, arXiv:math/0201013 [S].
  - Ahmed, De Loera, Hemmecke, "Polyhedral cones of magic cubes and squares",
    arXiv:math/0201108 [S].
  - Beck & van Herick, "Enumeration of 4×4 magic squares", arXiv:0907.3188 [S].
  - Beck & Zaslavsky, "Six little squares and how their numbers grow",
    J. Integer Seq. 13 (2010), arXiv:1004.0282 [S].
  - The Ehrhart series was known for n ≤ 6 (CTEuclid, then Todd-polynomial
    methods, arXiv:2304.13323) and has now been extended to n = 7, 8
    (arXiv:2607.21338, Jul 2026) [S].
  - Exact n = 3 case [V, by hand]: M_3(t) = 2t^2/9 + 2t/3 + 1 when 3 | t and 0
    otherwise. So M_3/H_3 → (2t^2/27)/(t^4/8) = 0.593/t^2, averaged over t.
    The Gaussian traversal heuristic of §7.1, with uniform-Dirichlet entry
    variance, gives (n+1)/(2π t^2) = 0.637/t^2: 7% high.

### 5.2 Counting additive magic squares: growth with the entry bound

* With entries in [1, B] and S free, magic squares form a lattice cone of
  dimension n^2 − 2n. Semi-magic squares form one of dimension n^2 − 2n + 2. So
  #magic ~ c B^{n^2−2n} and #magic/#semi-magic ~ B^{−2}: each diagonal costs ~1/B,
  like p_S ~ 1/S [M; linear algebra].
* **Circle method.** Rome & Yamagishi, "On the existence of magic squares of
  powers", arXiv:2406.09364, Res. Number Theory 11 (2025) 91 [S]: n×n magic
  squares of d-th powers exist for all n ≥ n_0(d), and of squares for all n ≥ 4.
  They describe this as the first use of the circle method on magic squares.
  Flores, "A circle method approach to K-multimagic squares", J. London Math.
  Soc. (2025), doi:10.1112/jlms.70290, arXiv:2406.08161 [S], proves
  M_{K,N}(P) ~ c·P^{N(N−K(K+1))} for N > 2K(K+1). For K = 1 this is the
  N^2 − 2N above. arXiv:2411.01091 extends it to distinct entries [S].
  There is **no analogue for multiplicative constraints**, which are not
  polynomial equations of bounded degree in the entries.
* **One multiplicative equation** [M, heuristic; n = 2 is classical]. The number
  of solutions of x_1⋯x_n = y_1⋯y_n with all variables ≤ B should be
  ≍ B^n (log B)^{(n−1)^2}: parametrize by an n×n "factor matrix" z_ij with
  x_i = ∏_j z_ij and y_j = ∏_i z_ij. For n = 2 this is c·B^2 log B. So one
  product equation costs ~B^{−n} up to the large (log B)^{(n−1)^2} factor.
  Naive independence over the 2n+1 product equations of a square then gives a
  negative exponent of B. Solutions exist anyway, which shows the count is
  dominated by smooth, highly composite P. That is why searches fix P with many
  small prime factors. **Any analytic growth estimate must go through the
  distribution of smooth numbers / τ(P), not a box count.**

### 5.3 Exact counts of normal squares (entries 1..n^2), used in §7.1

| n | magic squares (up to symmetry) | semi-magic incl. magic (up to symmetry, A271103) |
|---|---|---|
| 3 | 1 | 9 |
| 4 | 880 | 68,688 |
| 5 | 275,305,224 (Schroeppel 1973) [M] | 579,043,051,200 [S] |
| 6 | 17,753,889,197,660,635,632 (H. Mino, 2024, ~80,000 GPU-h) [S]; earlier estimate (1.7745 ± 0.0016)e19 (Pinn & Wieczerkowski 1998, cond-mat/9804109) [S] | 94,590,660,245,399,996,601,600 (Ripatti, arXiv:1807.02983) [S] |

Sources: OEIS A271103 / A271104 [S] and W. Trump, "How many magic squares are
there?", https://www.trump.de/magic-squares/howmany.html [S].

---

## 6. Constructions, Latin squares, and Euler's 36 officers

* **Horner** (Scripta Mathematica, 1952 for 9x9 and 1955 for 8x8 [S]; volume and
  page numbers 18:300–303 and 21:23–27 are [M]). Per Keedwell–Dénes-style
  summaries: "Horner (1952) showed how to construct addition-multiplication magic
  squares of any odd order and in Horner (1955) he obtained ... orders 8 and 16.
  Both papers make use of Latin squares" [S]. "Any odd order" cannot literally
  include 3 (impossible) or 5 and 7 (unknown before 2016). It presumably means
  composite or large odd orders such as 9 = 3·3.
* The classical multiplicative construction [M] puts M_ij = a^{A_ij} b^{B_ij},
  where A and B are orthogonal diagonal Latin squares. Every line contains every
  symbol of A and of B once, so the line products are equal. Orthogonality makes
  the entries distinct. Add-mult constructions (Horner, Zhu, Boyer, Shirakawa)
  combine such Latin structure with extra conditions that make the sums equal
  as well.
* **Order 6 has no Graeco-Latin square.** Euler (1782) posed the 36 officers
  problem. Tarry (1900) proved impossibility by exhaustion (9,408 cases), and
  Stinson gave a 3-page proof (J. Combin. Theory A 36 (1984) 373–376) [S].
  Bose, Shrikhande and Parker (1959–60) disproved Euler's conjecture for every
  other n ≡ 2 mod 4 > 6 [M]. Orthogonal diagonal Latin squares exist for all n
  except 2, 3, 6 [M] (Brown, Cherry, Most, Parker, Wallis 1992). The "quantum 36
  officers" exist (Rather et al., PRL 128 (2022) 080507; arXiv:2204.06800) [S],
  but that does not help with integers.
  **Consequence:** at n = 6 the Latin-square routes to multiplicative structure
  (exponent matrices E_q = permutation-matrix sums organized by two orthogonal
  Latin squares) are unavailable. The 4k+2 orders 10, 14, 18, 22, 26 do have
  add-mult examples (Zhu 1996; Shirakawa 2010).
  That fits the data in §7: 6x6 semi-magic squares behave as generic random
  objects, with Poisson traversal counts and Gaussian traversal sums. No
  structured subfamily is known that would make diagonals cheap.
* Every known 5x5, 6x6 and 7x7 add-mult object (Morgenstern 2007–08,
  Shirakawa 2010, Miquel 2016, our 2023/2026 searches) was found by **computer
  search**, not by algebraic construction [S][R].
* Off-topic but nearby: Mehat, "A proper Euler magic matrix of order 6"
  (arXiv:2608.15318, 2026) is about orthogonal integer matrices whose squared
  entries are magic. It does not construct add-mult squares [S].

---

## 7. Calibration checks done here (cheap; < 3 CPU-min total)

### 7.1 Does P(magic | semi-magic) ≈ p_S^2 hold for normal squares?

For any semi-magic square all row and column means are S/n. The sum over a
uniformly random traversal therefore has exactly Var = n^2/(n−1) · var(entries).
The local-CLT density at S is p_S ≈ (2π Var)^{−1/2}. If the two main diagonals of
a random semi-magic square behaved like independent random traversals, then
P(magic | semi-magic) = p_S^2. Exact counts (§5.3) give [V]:

| n | S | p_S (Gaussian) | p_S^2 | observed #magic/#semi | obs / p_S^2 |
|---|---|---|---|---|---|
| 3 | 15 | 0.0728 | 5.3e-3 | 0.111 | 21 |
| 4 | 34 | 0.0375 | 1.40e-3 | 0.0128 | 9.1 |
| 5 | 65 | 0.0221 | 4.90e-4 | 4.75e-4 | **0.97** |
| 6 | 111 | 0.0143 | 2.05e-4 | 1.88e-4 | **0.92** |

From n = 5 on, the naive "two independent Gaussian diagonals" model is right to
3–8%, with no fitted parameter. At n = 3 and 4, structured squares (Latin/group
based, with lattice-concentrated traversal sums) dominate and inflate the ratio
9–21x. For 6x6 add-mult squares the analogue of that structure is ruled out
twice over: by §6, and by the Poisson s/p-counts in research/ideas.md. So
pair-correlation factors of ~1 (the forecast uses S+S 1.1–1.2) are the right
order.

### 7.2 Gaussian (local CLT) prediction of traversal rates for our add-mult squares

Script: `gauss_check.py` and `gauss_full.py` (in this directory). Data: 6,498
distinct 6x6 semi-magic squares from `forecast/state/units` and
`verify-forecast/*` (P near S_min, N ≈ 1000–8000).

For each square the script computes the exact covariance, over the 720
traversals, of the vector (traversal sum, exponent sum of each prime of P):
Cov = Xc·Xcᵀ/(n−1). It then predicts the density at the target lattice point
(S, e_2, e_3, ...) as (2π)^{−d/2} det(Cov)^{−1/2}. The exact counts were
recomputed and match the stored s_count / p_count for all squares [V].

| quantity (summed over 6,498 squares) | observed | prediction | obs/pred |
|---|---|---|---|
| S-traversals | 8,089 | 1-D Gaussian: 8,256 | **0.98** |
| P-traversals | 2,008 | product of exact per-prime marginals (primes independent): 946 | 2.12 |
| | | product of 1-D Gaussians per prime: 1,220 | 1.65 |
| | | joint Gaussian with cross-prime covariance: 2,893 | **0.69** |
| SP-traversals | 14 | independence on each square's observed s, p (Σ s·p/720): 3.29 | **4.3** (90% CI 2.6–6.7) |
| | | joint Gaussian (sum + all exponents): 17.1 | **0.82** (90% CI 0.49–1.28) |

By number of primes in P (4 / 5 / 6), P-traversals obs/joint-Gaussian are 0.84 /
0.69 / 0.56, and SP-traversals obs vs joint-Gaussian are 3/2.8, 8/7.9, 3/6.4.

Takeaways:

* **The sum is understood.** p_S = (2π·(36/5)·var(entries))^{−1/2} is right to 2%.
  For a fixed P, var(entries) ∝ S^2 roughly, which reproduces the report's
  empirical p_S ∝ 1/S.
* **The product is not independent across primes.** Big entries carry many
  prime factors, so exponent sums of different primes are positively correlated
  (factor ~2.1 over independence). The full Gaussian overshoots by ~1.4x on
  average, and by more with more primes. Primes with e_q = 1 are exactly
  "permutation-matrix" variables, and P(exactly one hit) = 264/720 = 0.367 is far
  from Gaussian in the joint tail. The truth lies between the two bounds,
  ~0.7x the full Gaussian.
* **ρ has a mechanism and a predictor.** Traversals with the right exponents
  also have a more typical sum, because sum and log-size are positively
  correlated. The joint Gaussian puts the S–P association at 17.1/4.98 = **3.4**.
  The observed association relative to each square's own marginals is 4.3
  (2.6–6.7). The joint Gaussian predicts the absolute SP count to 0.82
  (0.49–1.28) with no fitted parameter. A per-square analytic estimator of
  P(magic) is then ≈ 5400 · pair · p_SP(Gaussian)^2 × ~0.8^2. It can be
  evaluated on every square, including the 758,949 legacy squares if their grids
  are available. That removes the need to fit ρ from 14 events.

### 7.3 Ladder across orders: where the first objects appeared

Computed with `bin/enumerate --vec-size n` [V]:

| n | object | S | S_min(P) | S/S_min | N vectors at S |
|---|---|---|---|---|---|
| 4 | smallest-S semi-magic | 247 | 172 | 1.44 | 8 (exactly its rows + cols) |
| 4 | smallest-Max semi-magic | 325 | 319 | 1.02 | 8 |
| 5 | Morgenstern semi-magic | 476 | 404 | 1.18 | 78 |
| 6 | smallest-S semi-magic (S+P) | 289 | 206 | 1.40 | 483 |
| 6 | smallest-Max semi-magic | 316 | 265 | 1.19 | 516 |
| 6 | smallest-P semi-magic | 327 | 171 | 1.91 | 451 |
| 7 | smallest semi-magic (S+S) | 310 | 200 | 1.55 | 3,834 |
| 7 | **magic** (Miquel) | 465 | 279 | 1.67 | 8,963 |

Use as anchors:

* The smallest semi-magic squares sit at N ≈ 8, 78, ~500 and ~3,800 for
  n = 4, 5, 6, 7.
* At n = 7 the first magic square appeared at S_magic / S_semi-min = 465/310 = 1.5,
  N_magic / N_semi ≈ 2.3, after ~600 single-PC hours.
* Count n^2 cells against 2(2n+1) constraints (2n+1 for the sums, 2n+1 for
  the products). For n = 5, 6, 7 that is 25 vs 22, 36 vs 26 and 49 vs 30, so
  the free fraction is 0.12, 0.28 and 0.39. This is the usual reason smaller
  orders are harder [R].
* The 6x6 records cluster at S = 289–327 with N ≈ 450–520, and nobody has
  improved on Morgenstern's 2007 records (S = 289, Max = 105, P = 5.08e8). Our
  search rediscovered the Max-105 square.

---

## 8. Corrections to the task brief

* "Known smallest 6x6 semi-magic: Morgenstern's P = 2^10 3^4 5^3 7^2, S = 327" is
  the **smallest-P** record. The smallest **S** is 289 (P = 2^11 3^4 5^3 7 11,
  S+P type, Max 135). The smallest **Max** is 105 (S = 316). All three are
  Morgenstern 2007 [S][R][V].
* "None is known for order 6 (nor 5)": correct as of Oct 2026 [S]. Add: 3x3 and
  4x4 are proved impossible; 7x7 exists (2016); 8x8 and 9x9 have existed since the
  1950s; many larger orders are known.

---

## 9. References (URLs as found; most could not be opened directly)

1. C. Boyer, Smallest additive-multiplicative magic square. http://www.multimagie.com/English/SmallestAddMult.htm (fr: /Francais/SmallestAddMult.htm)
2. C. Boyer, Additive-multiplicative magic squares, 8th and 9th order. http://www.multimagie.com/English/Multiplicative8_9.htm ; 10th order and up: http://www.multimagie.com/English/AddMult10_.htm
3. C. Boyer, Smallest multiplicative magic squares. http://www.multimagie.com/English/Multiplicative.htm ; 6th and 7th order: http://www.multimagie.com/English/Multiplicative6_7.htm
4. C. Boyer, Unsolved multimagic problems. http://www.multimagie.com/English/Problems.htm ; Enigmas: http://www.multimagie.com/English/Enigmas.htm
5. L. Morgenstern, 4x4 magic/semi-magic impossibility proof. http://www.multimagie.com/English/Morgenstern05.htm
6. multimagie news: Nov 2007 http://multimagie.com/English/News0711.htm ; Jul 2008 http://multimagie.com/English/News0807.htm ; Feb 2009 http://multimagie.com/English/News0902.htm ; Aug 2013 http://www.multimagie.com/English/News1308.htm ; Horner biography http://www.multimagie.com/English/Horner.htm
7. Plus magazine (2010), Win money with magic squares. https://plus.maths.org/os/latestnews/may-aug10/magic/index
8. D. Weisenberg, Some Thoughts on the Search for 5×5 and 6×6 Additive-Multiplicative Magic Squares, arXiv:2311.06326; Minnesota J. Undergrad. Math. 8(1) (2024). https://arxiv.org/abs/2311.06326 , https://pubs.lib.umn.edu/index.php/mjum/article/view/6015
9. J. Cilleruelo, F. Luca, On multiplicative magic squares, Electron. J. Combin. 17 (2010) #N8. https://www.combinatorics.org/ojs/index.php/eljc/article/view/v17i1n8
10. D. Borkovitz, F. K. Hwang, Multiplicative magic squares, Discrete Math. 47 (1983) 1–11. https://www.sciencedirect.com/science/article/pii/0012365X83900675
11. W. W. Horner, Addition-multiplication magic squares, Scripta Math. 18 (1952); Addition-multiplication magic squares of order 8, Scripta Math. 21 (1955). [M for volume/pages]
12. Zhu Jiacheng, A construction of addition-multiplication magic square of order 18, J. Statist. Plann. Inference 51(3) (1996) 331–337. https://www.sciencedirect.com/science/article/abs/pii/0378375895000968
13. K. Boonsiriseth, CJ Quines, S. Raphael, Searching for Additive-Multiplicative Magic Squares (2023), /home/user/magic-squares/project-report.pdf ; https://github.com/talkon/magic-squares
14. OEIS A114060, A113026, A271103, A271104. https://oeis.org/A114060 etc.
15. A. Ripatti, On the number of semi-magic squares of order 6, arXiv:1807.02983.
16. K. Pinn, C. Wieczerkowski, Number of magic squares from parallel tempering Monte Carlo, Int. J. Mod. Phys. C 9 (1998) 541, arXiv:cond-mat/9804109; W. Trump, How many magic squares are there? https://www.trump.de/magic-squares/howmany.html (Mino 2024 exact 6x6 count).
17. M. Beck, D. Pixton, The Ehrhart polynomial of the Birkhoff polytope, Discrete Comput. Geom. 30 (2003) 623–637, arXiv:math/0202267; The volume of the 10th Birkhoff polytope, arXiv:math/0305332.
18. E. R. Canfield, B. D. McKay, The asymptotic volume of the Birkhoff polytope, arXiv:0705.2422; Asymptotic enumeration of integer matrices with large equal row and column sums, Combinatorica 30 (2010), arXiv:math/0703600. A. Barvinok, J. A. Hartigan, arXiv:0709.3810.
19. M. Beck, M. Cohen, J. Cuomo, P. Gribelyuk, The number of "magic" squares, cubes, and hypercubes, Amer. Math. Monthly 110 (2003) 707–717, arXiv:math/0201013. M. Ahmed, J. De Loera, R. Hemmecke, arXiv:math/0201108. M. Beck, A. van Herick, arXiv:0907.3188. M. Beck, T. Zaslavsky, arXiv:1004.0282. Ehrhart series of magic squares of orders 7 and 8, arXiv:2607.21338.
20. N. Rome, S. Yamagishi, On the existence of magic squares of powers, arXiv:2406.09364, Res. Number Theory 11 (2025) 91. D. Flores, A circle method approach to K-multimagic squares, J. London Math. Soc. (2025), arXiv:2406.08161; Existence of K-multimagic squares and magic squares of kth powers with distinct entries, arXiv:2411.01091.
21. G. Tarry (1900); D. R. Stinson, A short proof of the nonexistence of a pair of orthogonal Latin squares of order six, J. Combin. Theory A 36 (1984) 373–376; Bose–Shrikhande–Parker (1960) [M]; S. A. Rather et al., Thirty-six entangled officers of Euler, PRL 128 (2022) 080507 [M]; arXiv:2204.06800.
