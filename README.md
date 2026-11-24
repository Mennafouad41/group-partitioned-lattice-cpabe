# Group-partitioned lattice CP-ABE: code and benchmarks

Code and raw results for the paper **"A Novel Zero-Trust Communication Protocol for Data Access Control in Medical Metaverse Applications"** (Fouad, Taha, Kotb).

This is a research prototype. The scheme is a group-partitioned lattice CP-ABE built on top of the Zhang-Zhang scheme as implemented in PALISADE-ABE. The security argument is in the paper and its supplementary material.

## Files

| File | What it is |
|---|---|
| `lattice_daemon3.cpp` | The scheme code (group-wise additive sharing, G groups, default relation). |
| `tree_bench.cpp` | Benchmark on the policies of a tree file. Runs the baseline (Zhang-Zhang = the same code with G = 1) and our scheme (G = 4). |
| `tree_med32.txt` | The medical-metaverse access tree: 16 policies, 32 attribute positions, 4 groups. |
| `tree_results_1794476410.csv` | Results used in Table 2 (mean over the 16 policies). The other `tree_results_*.csv` files are earlier runs. |
| `hybrid_tree_bench.cpp`, `hybrid_tree_bench.py` | Benchmark of the pairing-based TKN20 and AC17 (FAME) on the same tree. |
| `*-hybrid_results_*.csv` | Results of that benchmark (TKN20: `863372b4-...`, AC17: `d719af81-...`). |
| `run.log` | Console log of a run. |

## Requirements

* [PALISADE-ABE](https://gitlab.com/palisade) installed in `/usr/local` (headers in `/usr/local/include/palisade-abe`, `/usr/local/include/palisade`, `/usr/local/include/palisade/core`).
* For the pairing benchmark: the Hybrid-PQ-CP-ABE-Library, built in `~/Hybrid-PQ-CP-ABE-Library`.
* g++ with C++17. Tested on Ubuntu (WSL2).

## Lattice benchmark (baseline vs. our scheme)

```
g++ -O2 -std=c++17 -o tree_bench \
  -I/usr/local/include/palisade-abe -I/usr/local/include/palisade -I/usr/local/include/palisade/core \
  -L/usr/local/lib tree_bench.cpp -lPALISADEabe -lPALISADEbinfhe -lPALISADEpke -lPALISADEcore
./tree_bench tree_med32.txt 15
```

The arguments are the tree file and the number of repetitions. The program prints a table and writes `tree_results_<time>.csv`. If the PALISADE libraries are not found, run `export LD_LIBRARY_PATH=/usr/local/lib`.

## Pairing benchmark (TKN20, AC17)

```
LIB=$HOME/Hybrid-PQ-CP-ABE-Library
g++ -O2 -std=c++17 -I$LIB/include hybrid_tree_bench.cpp -L$LIB/lib -lhybrid-pq-cp-abe -Wl,-rpath,$LIB/lib -o hybrid_tree_bench
./hybrid_tree_bench tree_med32.txt 15
```

It writes `hybrid_results_<time>.csv`.

## Results reported in the paper (mean over the 16 policies, ms; ciphertext in KB)

| Scheme | KeyGen | Enc | Dec | Ciphertext |
|---|---|---|---|---|
| Baseline (Zhang-Zhang, PALISADE-ABE) | 629 | 545 | 6.5 | 40,089 |
| Ours | 1,034 | 380 | 5.0 | 27,804 |
| TKN20 (pairing) | 151 | 92 | 25.9 | 5.0 |
| AC17 (pairing) | 70 | 86 | 49.0 | 11.6 |

Measured on a Dell Latitude 7490 (Intel Core i7-8650U @ 1.90 GHz, 16 GB RAM, Windows 10 Pro with WSL Ubuntu).

## License

BSD 2-Clause (see `LICENSE`). The code is derived from the PALISADE-ABE library (Copyright (c) 2019, NJIT, BSD 2-Clause, provided "as is" without warranty); the NJIT notice applies to the parts derived from it.
