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
  complex.pony   # Complex[F] — generic over F32/F64
  fft.pony       # FFT — used by MPInt for fast multiplication
  mpint.pony     # MPInt — arbitrary-precision integers (base-65536 U16 digits)
  mpfloat.pony   # MPFloat — arbitrary-precision floats (base-256 U8 digits, class val)
  primes.pony    # Prime[A], PrimeSieve, SegmentedSieve
  modular.pony   # Modular arithmetic, GCD/LCM
  limits.pony    # FLimits — platform floating-point limits
  unsigned_complement.pony  # Bit manipulation utilities for unsigned integers
  gmp/           # FFI wrappers for libgmp / libmpfr (compiled separately)
tests/           # Unit tests (package compiled separately from mathx/)
  tests.pony     # Main test list / entry point
  _tests_*.pony  # Per-component test files
examples/        # Standalone Pony programs demonstrating the library
drafts/          # Work-in-progress, not compiled by default
```

### Key type relationships

- `Complex[F: (Float & FloatingPoint[F]) = F64]` is `class val`, implements `Approximated[Complex[F], F]` and `Stringable`.
- `MPFloat` is `class val` with fields `_sign`, `_nan`, `_inf`, `_exponent: I64`, `_digits: Array[U8] val`. Uses Newton's method for `inv()`, `sqrt()`, and `pi()`.
- `MPInt` uses base-65536 `U16` digits; FFT multiplication is limited to ~10^6 digits due to F64 precision.
- `Prime[A: UnsignedInteger[A] val = USize]` is a primitive; `PrimeSieve` and `SegmentedSieve` are classes using `BitMap` from `bitsx`.

### Approximate equality

`almost_eq` is defined on `Complex[F]` but **not** on bare `F32`/`F64`. To use `assert_almost_eq` with float scalars, wrap them: `Complex[F](f)`.

## Pony Gotchas in This Codebase

- `F64.nan()` / `F64.inf()` are predicates (`Bool`), not constructors. Use `F64.from_bits(...)` to create special values.
- Inside a `fun box`, `this` has `box` capability; return type `T^` (ephemeral val) requires constructing a fresh value, not returning `this`.
- Variables in arithmetic loops must be `val^`/`iso`; `var x = this` inside a method gives `box` and breaks with `val^` receivers. Use `var x: T = T(...)`.
- Cannot name a local variable `pi` inside `MPFloat.from_string` — it shadows the `pi` constructor. Use `pos` instead.
- `int_exp - frac_count + str_exp` in `from_string` requires explicit parens: `(int_exp - frac_count) + str_exp`.
