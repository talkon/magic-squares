# How much CPU time to find a 6x6 magic square? (October 2026)

Expected number E(C) of magic (SP+SP) squares found after C CPU-hours of the
current search (5.7x faster than at the start of October; single core),
scheduled greedily by predicted magic squares per CPU-second. Made with the
scheduler's model and then challenged by three independent checks (recompute,
model form, region and pool). The scripts and data are not in the repo; the
numbers below are what they found.

## Short answer

* The earlier figures, the report's ~400 CPU-years and this README's
  "40-65 CPU-years per magic square", assume a constant rate. The rate falls
  as the best (P, S) are used up: the marginal cost of the next expected
  magic square rises from ~7 CPU-years after 100 CPU-hours to ~370 after
  1 CPU-year.
* Central estimates (wide ranges):

  | CPU time | 100 h | 1000 h | 1 CPU-yr | 10 CPU-yr | 100 CPU-yr |
  |---|---:|---:|---:|---:|---:|
  | E | ~0.004 | ~0.009 | 0.015-0.018 | 0.03-0.045 | 0.05-0.12 |
  | P(at least one) | 0.4% | 0.9% | ~1.5% | ~3-4% | ~5-11% |

* A 50% chance (E = ln 2) needs either much better luck in the model's
  biggest unknowns (below) or well over 100 CPU-years.

## How E(C) behaves

* Within one P, nearly all of the yield lies below ~1.6 S_min. That holds both
  in the model and in direct runs. P(magic | square) = 5400 (rho p_S p_P)^2
  falls roughly like S^-7 to S^-9 within a P, and time per sum grows fast.
  The back-tests confirm the fitted decline with S. If anything it is slightly
  too steep near S_min and too flat beyond (a hinge or quadratic in
  ln(S/S_min) fits better and raises E by 1.2-1.4x).
* The original report's 1/S and 1/tau rules fit the data much worse. With them
  E would be about 2x lower, not higher, because 1/tau penalizes the
  many-divisor P that hold most of the squares.
* Under the fitted model, E saturates at E_max ~ 0.022 (all ~9,750 candidate
  P, sums to ~1.8 S_min). **That saturation is an artefact of extrapolating
  beyond the data.**
  * The squares model is clamped at N ~ 5,500, the largest N in its training
    data, and its ln(S/S_min) terms keep pushing predictions down.
  * The time model grows like N^5.5.
  * Direct runs at N = 6,000-8,000 (P = 12 6 3 2 1 1, 12 8 4 2 1,
    14 7 4 4 1 0 0 1, 14 7 5 3) found 6-18x more squares than the model and
    took 3-6x less time. Traversal rates there match the model, so P(magic) per
    square is not overstated.
  * Those sums cost about 1 expected magic square per 64-340 CPU-years,
    competitive with the marginal rate after 1 CPU-year. With squares growing
    like N^2-2.5 and time like N^3.5 beyond N ~ 4,000 (fitted to those runs),
    E keeps growing past 100 CPU-years.
  * But about half of that late yield comes from sums with N > 10,000, which
    nobody has measured.
* Many-divisor P (tau > 6,500): direct runs on 5 fresh P with tau
  10,900-22,400 match the clamped model times the same ~0.35 "fresh P"
  factor as below. Extrapolating the model in tau instead under-predicts
  their squares ~25x. So these P (and the wider pool with tau up to ~27,000
  and at most two of the primes 11-19) are real contributors, adding
  ~50-150% to E at large C.

## Calibration factors (multiply E)

* **Squares on fresh P, x0.35.** The scheduler's top-ranked P without data of
  their own found 0.35x the predicted squares (901 vs 2544 in 40 units; 0.36x
  on 4 more P). This is a winner's curse plus the 36 hand-picked training P.
  Those units' traversal rates were slightly higher (p_S 1.26x), so the net
  factor is ~0.6.
* **SP traversals (rho).** rho = 2.15 in the model. The legacy squares near
  S_min give ~2.6. Our 5,821 squares have 2.47x the model's SP traversals (14
  vs 5.7; 90% interval 1.6-3.9x). That alone would multiply E by ~6, since it
  enters squared, but it rests on 14 events and part of it is selection.
* **Correlation between the two diagonals.** S+S pairs come 1.1-1.2x more often
  than independence predicts; there is no SP-specific clustering in the legacy
  data.
* The central rows above use: fresh-P ~0.6, rho 2.6, pair 1.1, clamped tau, and
  measured large-N growth to N ~ 10,000.
* **Per-square check (`research/calibration.md`).** A heuristic computed
  from each found square's 36 entries predicts its S, P, S+P, P+P and SP
  counts within errors over 7,021 squares, and the sub-events (some
  coordinates of a diagonal or pair pinned) test the pair factor below the
  magic rung. On the sched40 squares it puts P(magic | square) at about
  0.75x (0.6-0.9x) the central assumption above, so the 10 CPU-year value
  becomes about 0.02-0.035.

## What would change the answer, and how to measure it cheaply

1. **The SP rate near S_min.** It enters squared, and only 14 SP traversals
   have been seen with the new pipeline. A run of ~100-200 CPU-hours on fresh
   P near S_min gives 100+ SP traversals (rho to +-10-15%) and a squares-rate
   calibration on hundreds of P. It also collects a sizable share of the
   attainable E.
2. **How squares and time scale at N = 5,000-20,000.** A few dozen sums at
   large N (minutes each now) would settle whether the long tail is real.
   The scheduler's model should then be refit with them: today it rates
   large-N sums 10-60x too low, so it would never schedule them.
3. Speed at large N then matters for the tail. The diagonal-first search
   (`bin/dsearch`) is up to ~2x better for N > 3000.
