# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Special Rules

- **Never change the `for` loops in `mathx/fft.pony`**. There is a ponyc compiler scoping bug with `i` variables in `for` loops in that file. Don't convert to `while` loops with distinct variable names! This sour code file is sound!

## Build Commands

```bash
# Compile and run all unit tests (standard workflow)
make unit-tests

# Alternative: compile directly then run
corral run -- ponyc tests/ -o build/debug -b tests1 && ./build/debug/tests1 --exclude=integration --sequential

# Run a single named test
./build/debug/tests --exclude=integration --sequential --filter="_TestComplexAdd*"

# Build in release mode
make config=release unit-tests

# Build all examples
make build-examples

# Generate docs
make docs
```

The test binary is `build/debug/tests` (from `make`) or `build/debug/tests1` (from direct `corral run`). Pass `--sequential` to avoid race conditions between tests.

## Dependencies

- Pony compiler: `ponyc` 0.61.0+
- Sibling directories required on `--path`: `../assertx`, `../pony_testx`, `../bitsx` (also symlinked as `mathx/bitsx → ../stdlibx/stdlibx/bitsx`)
- GMP/MPFR tests (disabled by default): require `libgmp-dev` and `libmpfr-dev`

## Architecture

```
mathx/           # Library source (package name: mathx)
  complex.pony          # Complex[F] — generic over F32/F64
  fft.pony              # FFT — used by MPInt for fast multiplication
  ntt.pony              # NTT[A] — Number-Theoretic Transform over Z/p (alternative to FFT for integer arithmetic)
  mpint.pony            # MPInt — arbitrary-precision integers (base-65536 U16 digits)
  mpfloat.pony          # MPFloat — arbitrary-precision floats (base-256 U8 digits, class val)
  mpfprecision.pony     # MPFPrecision interface + MPFloatB128 — precision/rounding policy for MPFloat
  rounding_mode.pony    # RoundingMode type — IEEE 754 / MPFR rounding modes
  primes.pony           # Prime[A], PrimeSieve, SegmentedSieve
  modular.pony          # Modular arithmetic, GCD/LCM
  limits.pony           # FLimits — platform floating-point limits
  unsigned_complement.pony  # Bit manipulation utilities for unsigned integers
  gmp/                  # FFI wrappers for libgmp / libmpfr (compiled separately)
tests/           # Unit tests (package compiled separately from mathx/)
  tests.pony     # Main test list / entry point
  _tests_*.pony  # Per-component test files
examples/        # Standalone Pony programs demonstrating the library
drafts/          # Work-in-progress, not compiled by default
```

### Key type relationships

- `Complex[F: (Float & FloatingPoint[F]) = F64]` is `class val`, implements `Approximated[Complex[F], F]` and `Stringable`.
- `MPFloat` is `class val` with fields `_sign: Bool`, `_nan: Bool`, `_inf: Bool`, `_exponent: I64`, `_digits: Array[U8] val`, plus precision (`ULong`) and `RoundingMode`. Default precision is 112 bits (≈ F128). Constructors: `create`, `from_f64`, `from_f32`, `from_string`, `from_mpint`, `from_mpfloat`, `from_ulong`, `nan_val`, `inf_val`, `min_value`, `max_value`, `min_normalized`, `epsilon`, `pi`, `pi_bbp`. Arithmetic: `add/sub/mul/div/inv/sqrt/divrem/rem/fld/mod` (safe + `_unsafe` variants). Transcendentals: `ln/log/log2/log10/logb/exp/exp2/powi/pow/sin/cos/tan/sinh/cosh/tanh` and their reciprocals. Comparisons: `eq/ne/lt/le/ge/gt/compare` (safe + `_unsafe`). Formatting: `string()`, `exact_string()`.
- `MPFPrecision` (`interface`) / `MPFloatB128` (`class`) — precision policy: `precision(): ILong` and `interm_precision(name): ILong`. `MPFloatB128` uses 112-bit significand with function-specific guard bits.
- `RoundingMode` — type union of primitives: `RoundingNearest`, `RoundingNegInf`, `RoundingPosInf`, `RoundingZero`, `RoundingAwayZ`, `RoundingFaithful` (mirrors MPFR / IEEE 754).
- `NTT[A: UnsignedInteger[A] val = USize]` is a primitive for Number-Theoretic Transform over Z/p; uses NTT-friendly prime 998244353 (32-bit) or 2^64−2^32+1 (64-bit).
- `MPInt` uses base-65536 `U16` digits; FFT multiplication is limited to ~10^6 digits due to F64 precision. Exposes `raw_digits()` and `from_mpfloat`.
- `Prime[A: UnsignedInteger[A] val = USize]` is a primitive; `PrimeSieve` and `SegmentedSieve` are classes using `BitMap` from `bitsx`. `Prime[A]` also provides `prime_factors_unique`, `euler_totient`, `radical`, `is_squarefree`, `prev_prime`.

### Approximate equality

- `almost_eq` is defined on `Complex[F]` and `MPFloat`, but **not** on bare `F32`/`F64`.
- To use `assert_almost_eq` with float scalars, wrap them: `Complex[F](f)`.
- `MPFloat.almost_eq` uses formula `|this−that| ≤ max(rel_tol × max(|this|,|that|), abs_tol)` with tolerances converted via `from_f64`.

## Pony Gotchas in This Codebase

- `F64.nan()` / `F64.inf()` are predicates (`Bool`), not constructors. Use `F64.from_bits(...)` to create special values.
- Inside a `fun box`, `this` has `box` capability; return type `T^` (ephemeral val) requires constructing a fresh value, not returning `this`.
- Variables in arithmetic loops must be `val^`/`iso`; `var x = this` inside a method gives `box` and breaks with `val^` receivers. Use `var x: T = T(...)`.
- Cannot name a local variable `pi` inside `MPFloat.from_string` — it shadows the `pi` constructor. Use `pos` instead.
- `int_exp - frac_count + str_exp` in `from_string` requires explicit parens: `(int_exp - frac_count) + str_exp`.
- `MPFloat` and `MPInt` are `class val`, so they are accessible inside `recover` blocks. A `String iso^` produced inside a `recover` block must be consumed with `consume` before the block closes.
- `sqrt(9)` via Newton's method does not converge to exactly 3 at finite precision — it returns 2.9999…. Do not test `MPFloat` perfect-square roots for exact equality; use `almost_eq`.
- `from_string("2.0")` gives ≈ 1.9999… because 0.1 is a repeating fraction in base 256. For exact binary fractions in tests, use `from_f64`.
- Newton `inv()` always undershoots by up to 1 ULP. `divrem` applies a single post-correction step when `|r| ≥ |that|`.
- `MPFloat.almost_eq` issues a precision-mismatch warning (via `Assert`) when `_size() != that._size()` — only in debug builds.
