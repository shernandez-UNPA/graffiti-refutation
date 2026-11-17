# graffiti-refutation

Code and data for the computational study in

> S. Hernández, *Computers and Conjectures in Graph Theory: Proof, Discovery and Refutation*, submitted to *Theory and Applications of Graphs*.

The repository has two parts that share the **same** definitions of the 35 Graffiti
conjectures of Table 1 in Roucairol and Cazenave (2025) and of the invariants they use
(`de/src/conjectures.jl`, `de/src/invariants.jl`), so that heuristic and exhaustive results
are directly comparable.

| Folder | Content |
|---|---|
| `de/` | Differential evolution (DE) for counterexample search: random search, DE/rand/1/bin, memetic DE and variants with a smoothed objective, with a random-key encoding and an order curriculum. Produces Table 2 of the paper. |
| `exhaustive/` | Exhaustive verification: every graph of a class (connected, trees, triangle-free, girth ≥ 5) up to a given order, generated with nauty, is evaluated in parallel. Produces the minimum orders reported in the paper. |

Both parts use only the Julia standard library (no packages to install). Tested with Julia 1.13.
Code comments and the detailed READMEs (`README_es.md`) are in Spanish.

## Main results

| Graffiti | Class | Smallest order found by DE | Minimum order (exhaustive) |
|---|---|---|---|
| 29 | connected | 7 | 7 (unique) |
| 30 | connected | 9 | 9 (unique) |
| 137 | connected | 31 | > 10 |
| 139 | connected | 24 | > 10 |
| 166 | connected | 5 | 5 (unique: C₅) |
| 189 | connected | 4 | 3 (K₃; DE starts at n = 4) |
| 197 | triangle-free | 14 | 14 (unique: the path P₁₄) |
| 289 | girth ≥ 5 | 19 | **18** (six graphs) |
| 301 | tree | 14 | 14 (unique) |
| 302 | tree | 16 | 16 |
| 321 | triangle-free | 4 | 4 |
| 711 | connected | 5 | 5 |
| 715 | tree | 12 | 12 |

Exhaustive ranges: connected 3–10 vertices (11 989 762 graphs), trees 4–22 (9 114 282),
triangle-free 4–14 (466 452 571), girth ≥ 5 5–18 (2 298 165 270). No open conjecture of
the benchmark has a counterexample within these ranges. The full table is in
`exhaustive/resultados/graffiti/tabla.md`; the counterexamples are in the `.g6` files
(graph6/sparse6) under `exhaustive/resultados/graffiti/<class>/nNN/`.

## Reproducing

### DE study (`de/`)

```bash
cd de
julia test/runtests.jl
julia -t 8 scripts/experiment.jl --only refuted --runs 30 --time 300
julia scripts/summary.jl
```

See `de/README_es.md` for all options (open conjectures, smoothed variants, Wilcoxon tests).
The raw runs are in `de/results/runs.csv` and the counterexamples in `de/results/counterexamples/`.

### Exhaustive verification (`exhaustive/`)

Requires nauty (`geng`, `gentreeg`; on Ubuntu: `sudo apt install nauty`). On Windows, use WSL.

```bash
cd exhaustive
julia test.jl
bash correr_graffiti.sh any 3 10          # class, first order, last order
bash correr_graffiti.sh tree 4 22
bash correr_graffiti.sh k3free 4 13
bash correr_por_partes.sh k3free 14 10    # large orders, in resumable portions
bash correr_graffiti.sh girth5 5 17
bash correr_por_partes.sh girth5 18 50
julia tabla_graffiti.jl --de ../de/results/summary.csv
```

On an 8-core laptop (AMD Ryzen 7 5700U) the triangle-free graphs on 14 vertices take
about 1 h 45 min; the girth ≥ 5 graphs on 18 vertices take of the order of 15 h.

## License

MIT (see `LICENSE`).

## Citation

See `CITATION.cff`. If you use this code, please cite the paper above and the archived
version of this repository on Zenodo.
