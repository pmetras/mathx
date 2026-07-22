# Mathematical Constants Example

This example demonstrates the calculation of well-known mathematical constants
using arbitrary-precision floating-point numbers (`MPFloat` and `gmp.MPFloat`)
and the MPFR library bindings.

## Features

- **Borwein's Quadratic Algorithm:** Implementation of an algorithm that
  converges quadratically to π.
- **MPFR Integration:** Demonstrates using both the high-level `MPFloat`
  wrapper and the low-level `MPF` bindings to the MPFR library.
- **Series Calculation:** Calculates the base of the natural logarithm ($e$)
  using its Taylor series expansion.
- **Precision Control:** Shows how to specify and manage bit-precision for
  floating-point calculations.
- **Step-by-step comparison:** `pi_compare` runs both implementations in
  lockstep and prints aligned intermediate values so precision loss is
  visible at a glance.

## Usage

Build the example using the project's `Makefile`:

```bash
make build-examples
```

Run the binary:

```bash
./build/debug/constants
```

## Constants Calculated

1. **Pi ($\pi$):**
   - Calculated using `gmp.MPFloat.pi()` (direct MPFR call) — reference value.
   - Calculated using the Borwein quadratic formula with `MPFloat`.
   - Calculated using the Borwein quadratic formula with `gmp.MPFloat`.
   - Calculated using the Borwein quadratic formula with raw `MPF` bindings
     (lower-level API, minimum memory allocation).
2. **Euler's Number ($e$):**
   - Calculated using the series $\sum_{n=0}^{\infty} \frac{1}{n!}$ up to
     1 000 terms or until convergence.

## Accuracy Observations (`pi_compare`)

`pi_compare` runs the same Borwein quadratic iteration in lockstep for both
`MPFloat` and `gmp.MPFloat` at 10 000-bit precision (≈ 3 010 decimal digits).
It prints the last 35 significant digits of every intermediate per iteration,
right-aligned so differences are immediately visible:

```
── Iteration 10 ──────────────────────────────────
  pi_n+1  M: ...748940907186494231961567945240311    agree=3009 vs ref  |  3009 vs GMP
  pi_n+1  G: ...7489409071864942319615679452086      agree=3011 vs ref
  pi ref:    ...748940907186494231961567945208
```

### What the output shows

- Every intermediate (`sqrt(x)`, `x_n+1`, `y_n+1`) agrees between `MPFloat`
  and `gmp.MPFloat` to **full precision minus at most 1 digit** in every
  iteration. Both implementations correctly round each individual operation.
- The divergence from the true π appears only in `pi_n+1`, and only at the
  very last digits. After iteration 10 (the final one):
  - `gmp.MPFloat` agrees with the reference to **3011 digits**.
  - `MPFloat` agrees with the reference to **3009 digits**.
  - The gap is **2 decimal digits** (about 7 bits) out of 3 010.

### Why the gap exists

Each operation is correctly rounded to ≤ 1 ULP at output precision. Over the
~60 arithmetic operations in 12 Borwein iterations the rounding errors
accumulate. GMP internally uses 64-bit guard limbs (~64 guard bits); MPFloat
uses `⌊log₂(p)⌋ + 1` guard bytes (11 bytes = 88 bits at 10 000-bit precision).
Both are sufficient for a single operation to be faithfully rounded, but GMP's
rounding decisions happen to align with the true value 2 digits further than
MPFloat's in this particular computation.

This is an inherent property of round-to-nearest arithmetic over a chain of
operations: two correctly-rounded implementations can disagree at the last
1–2 digits. Eliminating it entirely requires Ziv's strategy (retry at higher
precision when rounding is ambiguous), which is not implemented here.

### Guard bytes

`_MPFAlgo._guard_bytes` controls the working precision added beyond the output
precision for each operation. The guard grows logarithmically with precision:

| Operation          | Guard                                     |
|--------------------|-------------------------------------------|
| `add`, `sub`       | `max(2, ⌊log₂ p⌋ + 1)` bytes             |
| `mul`, `inv`,`sqrt`| `max(4, ⌊log₂ p⌋ + 1)` bytes             |
| `ln`               | `max(6, ⌊log₂ p⌋ + 1)` bytes             |
| `exp`, `trig`, `pi`| `max(8, ⌊log₂ p⌋ + 1)` bytes             |

At 10 000 bits: 11 guard bytes for all basic operations.
At 1 000 000 bits: 17 guard bytes — the error budget scales correctly even
for extremely high-precision calculations.
