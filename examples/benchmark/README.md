# MPInt Benchmark Suite

Measures the wall-clock performance of core `MPInt` operations across a range of
operand sizes. The primary purpose is to track the impact of optimisations over
time: run the benchmark before a change, commit the result to a log file, apply
the change, run again, and compare.

## Quick start — shell script

The easiest way to record results is the provided shell script, which handles
header creation, file appending, and comment passing automatically:

```sh
# First baseline
./examples/benchmark/run_bench.sh "no optimisations yet"

# After a change
./examples/benchmark/run_bench.sh "automatic mul dispatch added"

# Only multiplication benchmarks, release build, 10 iterations
./examples/benchmark/run_bench.sh \
  --comment="NTT mul, PR #42" \
  --bench=mul,mul_karatsuba,mul_fft \
  --iter=10

# Debug build (Pony asserts enabled, no compiler optimisations)
./examples/benchmark/run_bench.sh --debug "debug baseline"

# Re-initialise the log file (keeps only the header)
./examples/benchmark/run_bench.sh --header-only

# Built-in help
./examples/benchmark/run_bench.sh --help
```

Results accumulate in `examples/benchmark/bench.csv`. The script creates the
file with a CSV header on the first run, then appends on every subsequent run.

## Manual usage

Each run prints one CSV row per `(benchmark, size)` pair to **stdout**.
Redirect or append the output to a persistent log file:

```sh
# Initialise a new log (write the header once)
./build/release/benchmark --header-only > bench.csv

# Append a release baseline
./build/release/benchmark --no-header --comment="baseline" --iter=10 >> bench.csv

# Implement an optimisation, then record the new result
./build/release/benchmark --no-header --comment="after dispatch opt" --iter=10 >> bench.csv

# Record a debug build for comparison (asserts enabled, no optimiser)
./build/debug/benchmark --no-header --comment="debug baseline" --iter=10 >> bench.csv
```

### Columns

| Column        | Type   | Description                                                   |
|---------------|--------|---------------------------------------------------------------|
| `ts_s`        | int    | Unix timestamp (seconds) when the run started                 |
| `build`       | string | `debug` or `release` — Pony `--debug` compile flag           |
| `bench`       | string | Benchmark name (see table below)                              |
| `n_digits`    | int    | Operand size in base-65536 digits (1 digit = 16 bits)         |
| `n_iter`      | int    | Number of timed iterations                                    |
| `elapsed_ns`  | int    | Total wall-clock nanoseconds for all `n_iter` iterations      |
| `ns_per_iter` | int    | `elapsed_ns / n_iter` — per-operation cost                    |
| `checksum`    | int    | XOR of result bytes; detects silent correctness regressions   |
| `comment`     | string | Free-text note from `--comment=TEXT`; empty string if omitted |

> **`build` vs directory name.** The `build` field reflects the Pony compiler
> `--debug` flag, not the output directory. `make build-examples` (default
> `config=debug`) produces a binary with `build=debug`; `make config=release
> build-examples` produces `build=release`.

## Benchmarks

| Name            | Method called           | Current complexity | Tracks optimisation |
|-----------------|-------------------------|--------------------|---------------------|
| `mul`           | `a * b`                 | O(n²)              | #1 dispatch         |
| `mul_karatsuba` | `a.karatsuba_mul(b)`    | O(n^1.585)         | reference           |
| `mul_fft`       | `a.fast_mul(b)`         | O(n log n)         | reference / #2 NTT  |
| `string`        | `a.string()`            | O(n²)              | #3 string           |
| `from_string`   | `MPInt.from_string(s)?` | O(n²)              | #3 string           |
| `pow2`          | `a * a`                 | O(n²)              | squaring proxy      |
| `divrem`        | `a.divrem(b)`           | O(n²)              | division            |
| `isqrt`         | `(a²+1).isqrt()`        | O(n² log n)        | sqrt / Newton       |

`mul` and `mul_fft` checksums are equal to `mul_karatsuba` checksums at the same
size: this cross-check catches algorithm bugs introduced during refactoring.

## CLI flags

| Flag                  | Default         | Description                                            |
|-----------------------|-----------------|--------------------------------------------------------|
| `--bench=n1,n2,...`   | all benchmarks  | Comma-separated list of benchmark names to run         |
| `--size=n1,n2,...`    | 10,100,500,2000 | Operand sizes in base-65536 digits                     |
| `--iter=N`            | 5               | Timed iterations per `(bench, size)` pair              |
| `--max-mem=N`         | 8               | Memory limit in GB; skip pairs that would exceed it    |
| `--comment=TEXT`      | (empty)         | Free-text note written to the `comment` CSV column     |
| `--no-header`         | off             | Suppress the CSV header line (for appending to a log)  |
| `--header-only`       | off             | Print CSV header and exit (to initialise a new log)    |

## Automatic size guards

Some benchmarks are automatically skipped for large sizes to prevent
multi-minute runtimes or excessive memory use. A warning is printed to stderr
and the CSV row is omitted.

| Benchmark       | Auto-skipped when | Reason                                      |
|-----------------|-------------------|---------------------------------------------|
| `mul`           | n > 500           | O(n²) schoolbook; > 500 takes minutes       |
| `isqrt`         | n > 200           | Squares operand first, then O(n² log n)     |
| any             | est. RAM > limit  | `--max-mem` guard (default 8 GB)            |

### Memory guidance (64-bit, 8 GB budget)

Each base-65536 digit occupies 2 bytes. Peak allocation per run:

- **`mul_fft`**: two FFT arrays of `next_pow2(2n)` F64 values ≈ 32 × n bytes.
  Safe up to n ≈ 250 million digits at 8 GB.
- **All other benchmarks**: O(n) bytes; the time budget is the binding
  constraint long before memory becomes an issue.

## Operand reproducibility

Both operands are always built from fixed seeds (`0xDEADBEEF` / `0xCAFEBABE`)
using the `Rand` PRNG, so repeated runs at the same size produce the same
checksums. Log entries with identical `(bench, n_digits, checksum)` triples
can be directly compared.

Record the `ponyc` version alongside the log; it is not embedded in the binary.
The `build` column (debug / release) is included in every CSV row.

## Analysis with Python / pandas

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('bench.csv')
df['datetime'] = pd.to_datetime(df['ts_s'], unit='s')

# Compare ns_per_iter across runs for all multiplication variants
mul = df[df['bench'].isin(['mul', 'mul_karatsuba', 'mul_fft'])]
pivot = mul.pivot_table(
    values='ns_per_iter',
    index='n_digits',
    columns=['build', 'bench', 'ts_s'],
)
pivot.plot(logy=True, logx=True, marker='o', figsize=(10, 6))
plt.title('MPInt multiplication: ns per iteration vs operand size')
plt.xlabel('Operand size (base-65536 digits)')
plt.ylabel('ns / iteration (log scale)')
plt.tight_layout()
plt.show()
```

## Build

```sh
make build-examples          # debug build → build/debug/benchmark
make config=release build-examples  # release build → build/release/benchmark
```
