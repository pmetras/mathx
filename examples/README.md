# Examples

Each subdirectory is a standalone Pony program that demonstrates one or more
features of the `mathx` library. All examples are built with:

```sh
make build-examples          # debug build  → build/debug/<name>
make config=release build-examples  # release → build/release/<name>
```

---

## [`benchmark`](benchmark/)

**Purpose:** Measure the wall-clock performance of `MPInt` operations to track
the impact of optimisations over time.

Emits one CSV row per `(benchmark, size)` pair to stdout. Rows are designed to
be appended to a persistent log file and analysed with a spreadsheet or Python /
pandas. Records Unix timestamp, `debug`/`release` build type, operand size,
iteration count, elapsed nanoseconds, a checksum, and a free-text **comment**
that describes what changed between runs (e.g. `"no optimisations yet"` →
`"added Karatsuba dispatch"`).

The `run_bench.sh` script automates the workflow: it creates `bench.csv` with a
header on the first run and appends on every subsequent run.

**Benchmarks:** `mul` (schoolbook dispatch), `mul_karatsuba`, `mul_fft`,
`string`, `from_string`, `pow2`, `divrem`, `isqrt`.

**Key flags:** `--bench=`, `--size=`, `--iter=`, `--comment=`,
`--no-header`, `--header-only`, `--max-mem=`.

**Classes:** `MPInt`

---

## [`benford`](benford/)

**Purpose:** Verify [Benford's Law](https://en.wikipedia.org/wiki/Benford%27s_law)
on mathematical sequences.

Counts leading-digit frequencies in Fibonacci numbers, factorials, and a linear
recurrence series, then compares them against the theoretical Benford
distribution. Uses `MPInt` for exact arithmetic and `MPFloat` for the
high-precision reference values.

**Classes:** `MPInt`, `MPFloat`

---

## [`collatz_modular`](collatz_modular/)

**Purpose:** Investigate the [Collatz conjecture](https://en.wikipedia.org/wiki/Collatz_conjecture)
using native integers and modular arithmetic.

Applies the Collatz function to increasing values of N using `ULong` and the
`Modular` helpers for odd/even tests and arithmetic. Runs indefinitely;
stop with Ctrl+C.

**Classes:** `Modular`

---

## [`collatz_mpint`](collatz_mpint/)

**Purpose:** Investigate the Collatz conjecture using arbitrary-precision integers.

Same algorithm as `collatz_modular` but uses `MPInt` so N can grow beyond the
native `I64` range. Demonstrates tail-recursive actor message passing to keep
the Pony runtime responsive during a long-running loop.

**Classes:** `MPInt`

---

## [`collatz_parallel`](collatz_parallel/)

**Purpose:** Parallelise the Collatz investigation across multiple Pony actors.

Splits the search space among a pool of actors, each computing Collatz sequences
with `MPInt` and reporting back to a coordinator. Demonstrates Pony's
actor model for CPU-bound concurrent workloads.

**Classes:** `MPInt`

---

## [`constants`](constants/)

**Purpose:** Compute well-known mathematical constants at arbitrary precision.

Calculates π via Borwein's quadratic algorithm and Euler's number e via its
Taylor series using `MPFloat`. Also demonstrates the low-level `MPF` bindings
to the MPFR library for maximum performance. Requires `libgmp-dev` and
`libmpfr-dev` libraries installed on the computer.

**Classes:** `MPFloat`, `MPF` (GMP/MPFR FFI)

---

## [`factstats`](factstats/)

**Purpose:** Analyse the digit distribution of factorials up to 1000.

Computes each factorial with `MPInt` and counts the frequency of each decimal
digit (0–9). Prints a table with per-row counts and running totals, checking
empirically that all digits occur with roughly equal frequency — consistent
with Benford's Law for the trailing-digit distribution.

**Classes:** `MPInt`

---

## [`fractal`](fractal/)

**Purpose:** Draw an ASCII Mandelbrot fractal using `Complex` arithmetic.

Iterates the recurrence z → z² + c over a grid of complex-plane coordinates
using `Complex[F64]`, and maps the divergence speed to ASCII characters.
A compact demonstration of the `Complex` API.

**Classes:** `Complex`

---

## [`platform_limits`](platform_limits/)

**Purpose:** Display the floating-point limits of the current platform.

Runs the MACHAR algorithm (ported from Fortran/C) to determine machine epsilon,
minimum/maximum normalised values, radix, and other IEEE 754 parameters for
`F32` and `F64`. Compares the computed values against the `FLimits` class.

**Classes:** `FLimits`
