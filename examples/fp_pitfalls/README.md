# Floating-Point Pitfalls: F64 vs MPFloat

This example demonstrates nine classic IEEE 754 floating-point failure modes,
showing each F64 failure side-by-side with an `MPFloat` mitigation.

## Building and running

From the repository root:

```bash
make build-examples
./build/debug/fp_pitfalls
```

Or compile directly:

```bash
corral run -- ponyc examples/fp_pitfalls/ -o build/debug -b fp_pitfalls
./build/debug/fp_pitfalls
```

## Sections

### 1. Decimal representation error

`0.1` has no exact representation in binary (or base-256). Adding `0.1` ten
times gives a sum that is not exactly `1.0` in F64. The rounding error
(~5.5e-17) compounds over 10 additions.

Three `MPFloat` input paths are compared:

| Input path | Error in `sum – 1` |
|---|---|
| `from[F64](0.1)` | ~5.5e-17 (inherits F64's approximation) |
| `from_string("0.1")` | ~10^-77 (decimal parsed directly to 256 bits) |

**Takeaway**: use `from_string` when the source value is a decimal literal;
`from[F64]` only when you already have an F64 value.

### 2. Catastrophic cancellation

`sqrt(x+1) - sqrt(x)` returns 0 for large `x` in F64 because `x` and `x+1`
round to the same bit pattern, making the two square roots identical. The
mathematically equivalent stable form `1 / (sqrt(x) + sqrt(x+1))` avoids the
near-cancellation. With `MPFloat` at 256 bits, `x` and `x+1` remain distinct
and the direct subtraction agrees with the stable form.

### 3. Non-associativity

F64 addition is not associative. With `a = 2e16`, `b = 1.0 - 2e16`, `c = 0.5`
(exact sum = 1.5):

- `(a + b) + c` → `0.5` (the 1.0 in `b` is lost when `b` is stored as F64)
- `a + (b + c)` → `0.0` (both 1.0 and 0.5 are lost)

Two different wrong answers. At 256 bits `b = 1.0 - 2e16` is represented exactly
and both groupings give `1.5`.

**Note**: non-associativity is intrinsic to any finite-precision arithmetic;
MPFloat greatly reduces but does not eliminate it for values that exceed its own
precision.

### 4. Quadratic formula instability

The standard formula `x = (-b ± sqrt(b²-4ac)) / 2a` suffers catastrophic
cancellation when `b² >> 4ac` (one root near zero). For the equation
`x² - 10000.0001x + 0.0001 = 0`, the small root `x₂ ≈ 1e-8` computed by
the standard formula loses most significant digits. The numerically stable
alternative `x₂ = c / (a·x₁)` avoids the near-cancellation. At 256 bits,
`MPFloat` carries enough guard digits that even the direct formula agrees with
the stable form.

### 5. Muller's recurrence

The recurrence `x(n) = 108 - (815 - 1500/x(n-2)) / x(n-1)` with `x(0) = 4`,
`x(1) = 4.25` mathematically converges to 5, but any tiny rounding error is
amplified each step toward the repelling fixed point at 100. In F64 the sequence
visibly diverges to 100 within 20 steps. At 256 bits `MPFloat` stays near 5
throughout.

**Reference**: J.-M. Muller, *Elementary Functions*, Birkhäuser, 1997.

### 6. Harmonic series partial sum

The harmonic series `H(n) = 1 + 1/2 + 1/3 + … + 1/n` diverges in mathematics
but converges in F64: once `1/k` falls below the rounding unit of the running
total (~10^-16), adding it has no effect. The F64 sum stalls permanently around
`n ≈ 10^15`. This section shows that even at `N = 100_000`, the forward and
backward F64 sums already differ, exposing how summation order matters.
`MPFloat` at 128 bits (~38 decimal digits) keeps every term accurate down to
`1/N` with no silent drops.

### 7. Rounding mode pitfalls

IEEE 754 defines five rounding modes, but F64 provides no per-operation API for
directed rounding. Two pitfalls:

1. **Interval arithmetic**: a guaranteed `[lo, hi]` enclosure requires
   `RoundingNegInf` for the lower bound and `RoundingPosInf` for the upper
   bound. `MPFloat` accepts a rounding mode in its constructor and propagates it
   through all operations. The example computes a verified enclosure for `sqrt(2)`
   and confirms `lo² ≤ 2 ≤ hi²`.

2. **Accumulated directed error**: summing `1/3` six times with `RoundingZero`
   always truncates, giving a sum provably less than 2; with `RoundingPosInf`
   every step rounds up, giving a sum provably greater than 2. These are
   verified lower and upper bounds.

### 8. Subnormal number pitfalls

IEEE 754 subnormal (denormalized) numbers fill the gap between zero and the
smallest normal value (`2^-1022 ≈ 2.225e-308` for F64). They trade mantissa
bits for extended range:

- `min_normal / 2` has only 52 significant bits instead of 53.
- `(x + 1.0) - 1.0 ≠ x` when `x` is subnormal: adding a subnormal to `1.0`
  loses all its mantissa bits (the ULP of 1.0 is ~2e-16, far larger than the
  subnormal value ~5e-324).
- Many CPUs flush subnormals to zero in hardware (FTZ/DAZ modes), causing
  results to differ silently between platforms.

**MPFloat has no subnormals**: every value is either exactly zero or a fully
normalized number with all `prec` bits significant, regardless of magnitude.
The exponent is `I64`, so the range is approximately `10^(±10^19)`. The example
stores `1e-320` (in the F64 subnormal range) with full 128-bit precision and
verifies that `(x/2)·2/2 = x/2` exactly.

### 9. Overflow and underflow

F64 overflows silently to `±∞` and underflows silently to `±0`, corrupting
downstream calculations. Three traps:

1. **Intermediate overflow in `hypot`**: `sqrt(a² + b²)` with `a = 1e200`
   overflows at `a²` even though the final result `~1.41e200` fits in F64.
   The stable form `a·sqrt(2)` avoids the squaring. `MPFloat` computes
   `sqrt(a² + a²)` directly without overflow.

2. **Factorial overflow**: `171!` overflows F64 to `+∞`; `170!` is the last
   representable value. `MPFloat` at 256 bits computes `200!` correctly.

3. **Underflow breaks reciprocals**: dividing the smallest positive F64
   (`≈5e-324`) by 2 underflows to 0; its reciprocal then becomes `+∞`.
   `MPFloat` stores `5e-350` (far below F64's minimum) as a finite value;
   dividing by 2 and taking the reciprocal gives the correct finite result.

## Key `MPFloat` features used

| Feature | Description |
|---|---|
| `MPFloat.from[F64](v, prec)` | Convert from F64, preserving its value |
| `MPFloat.from_string(s, prec)?` | Parse a decimal string directly, avoiding F64 approximation |
| `MPFloat.from[USize](n, prec)` | Convert from an integer type |
| `MPFloat(prec, RoundingMode)` | Construct with a specific rounding mode |
| `v.almost_eq(w, rel_tol, abs_tol)` | Approximate equality with relative and absolute tolerances |
| `v.is_infinite()` | Test for `±∞` |
| `v.string()` | Convert to a decimal string (base-256 scientific notation) |
| `v.exact_string()` | Exact decimal mantissa, decimal exponent, and inexact flag |

## Caveats

- `MPFloat.from[F64](0.1)` inherits F64's approximation error — always use
  `from_string` when the input is a decimal literal.
- `MPFloat` reduces but does not eliminate non-associativity; values that exceed
  its own precision still round.
- `string()` shows very small values (exponent below about -300) using base-256
  notation that may display as `0e-N`; the value is stored correctly. Use
  `exact_string()` for a decimal representation at any magnitude.
- The Muller recurrence shows that amplification instability eventually defeats
  any finite precision — MPFloat delays but does not eliminate divergence if
  precision is set too low.
