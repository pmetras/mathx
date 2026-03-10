# Collatz Conjecture with Modular Arithmetic

This example investigates the [Collatz conjecture](https://en.wikipedia.org/wiki/Collatz_conjecture) (also known as the Syracuse conjecture) using standard `ULong` integers and the `Modular` arithmetic utilities from `mathx/modular.pony`.

## The Conjecture

The Collatz function is defined for any positive integer $n$ as:
-   If $n$ is even: $n = n / 2$
-   If $n$ is odd: $n = 3n + 1$

The conjecture states that for any starting value $N \ge 1$, the sequence will always reach 1.

## Features

-   **Modular Arithmetic Utility:** Uses `Modular[ULong]` to perform calculations. 
    -   Even check: `Modular[ULong].sub_mod(n, 0, 2) == 0`.
    -   Next odd step: $3n+1$ is calculated modulo `ULong.max_value()` for safe wrapping.
-   **Step and Peak Tracking:** Calculates both the number of steps and the maximum value reached for each sequence.
-   **Periodic Recaps:** Prints progress every power of 10, highlighting the most complex $N$ and the highest peak reached.
-   **Efficiency:** Uses fixed-size `ULong` for faster calculations compared to arbitrary precision, suitable for initial ranges of $N$. Note that for very large $N$, the sequence might overflow $2^{64}-1$.

## Usage

Build the example using the project's `Makefile`:

```bash
make build-examples
```

Run the binary:

```bash
./build/debug/collatz_modular
```

To stop the calculation, use `Ctrl+C`.
