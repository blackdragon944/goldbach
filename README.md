# goldbach

Code accompanying the note "Goldbach partitions and the Sierpiński triangle" (arXiv: coming soon).

For the definitions, statements and results see the note; this repository contains only the
programs, tests, analysis scripts and small data files.

```sh
git clone https://github.com/blackdragon944/goldbach.git
cd goldbach
```

## Requirements

* a C compiler with POSIX threads: gcc or clang, or `zig cc` (`pip install ziglang`), which also
  works on Windows;
* Python 3.10+ with numpy ≥ 2.0 (the scripts use `np.bitwise_count`), matplotlib and pytest.

Tested on Windows 10 with Python 3.13, numpy 2.2, matplotlib 3.10, pytest 8.3 and zig 0.16.

## Building `src/astar.c`

```sh
# Linux / macOS
mkdir -p bin
cc -O2 -march=native -pthread -o bin/astar src/astar.c -lm

# Windows
pip install ziglang
mkdir bin
python -m ziglang cc -O2 -march=native -o bin/astar.exe src/astar.c
```

Modes (details at the top of `src/astar.c`):

| Command | What it computes |
|---|---|
| `astar scan N [threads] [dumpN dumpfile]` | a\*(n) for all even 8 ≤ n ≤ N: exceptions, records, histograms; optional uint16 dump of a\*(n) for n ≤ dumpN |
| `astar count A D outfile` | G_A(n) for all even n ≤ D, written as uint32 |
| `astar carries N [threads]` | for a = 3, the number of even 10 ≤ n ≤ N without a pair with at most c carries, c = 0..3 |
| `astar ext N [threads]` | a\*(n), the Chen-type a\*_C(n), the binary vertex rule a_bin(n) and the explicit vertex a_k(n) for all even 10 ≤ n ≤ N |
| `astar tries X [outfile]` | number of vertices tried until the first covering one, for the even n ≤ X not covered by a = 3 |
| `astar mult full N0 N1 [early]`, `astar mult sample N0 N1 count seed` | covering multiplicity M(p, n) of Goldbach pairs |

## Tests

```sh
python -m pytest tests
```

Tests that need `bin/astar` are skipped if it has not been built.

## Scripts and the note

Run everything from the repository root. The scripts write CSV files to `data/`, LaTeX macros and
table rows to `output/tex/` and figures to `output/figures/`.

| Script | Output | In the note |
|---|---|---|
| `src/analysis.py` | `output/figures/T3_n14.pdf` | Figure 1 |
| `src/analysis.py` | `output/tex/table_fixed_a.tex` | Table 1 |
| `src/analysis.py` | `output/figures/astar_distribution.pdf` | Figure 2 |
| `src/analysis.py` | `output/tex/table_records.tex` | Table 2 |
| `src/analysis.py` | `output/figures/records_growth.pdf` | Figure 3 |
| `src/analysis2.py` | `output/figures/multiplicity.pdf` | Figure 4 |
| `src/analysis.py` | `output/figures/astar_popcount_sparse.pdf` | Figure 5 |
| `src/analysis2.py` | `output/tex/table_binrule.tex` | Table 3 |
| `src/analysis2.py` | `output/tex/table_G3_popcount.tex` | Table 4 |
| `src/analysis.py` | `output/figures/G3_vs_n_and_heuristic.pdf` | Figure 6 |
| `src/analysis.py` | `output/tex/table_carries.tex` | Table 5 |
| `src/analysis2.py` | `output/tex/table_chen_records.tex` | Table 6 |
| `src/analysis.py`, `src/analysis2.py` | `output/tex/numbers*.tex` | all other numbers in the text |

`*_comma.tex` variants (decimal comma) are written next to each `.tex` file.
`src/analysis2.py` needs the outputs of the C runs listed at its top (for example
`astar ext 60000000000 4 > data/ext_6e10.txt`); the outputs used for the note are included in `data/`.

`src/goldbach_ref.py` is a simple Python reference implementation used by the tests.
`sierpinski.py` is an interactive, zoomable Sierpiński triangle (matplotlib).

## Running time and memory up to 6·10^10

Measured on a 4-core Intel i5-2500K (3.3 GHz) with 4 threads:

| Run | Time | Memory |
|---|---|---|
| `astar scan 60000000000 4` | sieve 44 s + scan 1202 s (about 21 min) | about 3.75 GB (bit sieve, N/16 bytes) |
| `astar ext 60000000000 4` | about 2 h 54 min | about 3.75 GB |

## Data files not included (too large)

| File | Size | How to get it |
|---|---|---|
| `data/astar_upto_1e8.bin` | 100 MB | `bin/astar scan 100000000 4 100000000 data/astar_upto_1e8.bin` (`src/analysis.py` does this itself if the file is missing) |
| `data/G3_upto_1e6.csv` | 6.7 MB | written by `python src/analysis.py` |
| `data/mult_early_1e6.txt` | 72 MB | `bin/astar mult full 10 1000000 1 > data/mult_early_1e6.txt` |

## License

MIT, see [LICENSE](LICENSE).
