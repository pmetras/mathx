# Mathematical Constants Example

This example demonstrates the calculation of well-known mathematical constants using arbitrary-precision floating-point numbers (`MPFloat`) and the MPFR library bindings.

## Features

- **Borwein's Quadratic Algorithm:** Implementation of an algorithm that converges quadratically to Pi.
- **MPFR Integration:** Demonstrates using both the high-level `MPFloat` wrapper and the low-level `MPF` bindings to the MPFR library.
- **Series Calculation:** Calculates the base of the natural logarithm ($e$) using its Taylor series expansion.
- **Precision Control:** Shows how to specify and manage bit-precision for floating-point calculations.

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

1.  **Pi ($\pi$):** 
    - Calculated using `MPFloat.pi()` (direct MPFR call).
    - Calculated using the Borwein quadratic formula with `MPFloat`.
    - Calculated using the Borwein quadratic formula with raw `MPF` bindings for maximum performance and minimum memory allocation.
2.  **Euler's Number ($e$):**
    - Calculated using the series $\sum_{n=0}^{\infty} \frac{1}{n!}$ up to 1000 terms or until convergence.
