# Collatz Conjecture with MPInt

This example investigates the [Collatz conjecture](https://en.wikipedia.org/wiki/Collatz_conjecture) (also known as the Syracuse conjecture) using arbitrary-precision integers (`MPInt`).

## The Conjecture

The Collatz function is defined for any positive integer $n$ as:
-   If $n$ is even: $n = n / 2$
-   If $n$ is odd: $n = 3n + 1$

The conjecture states that for any starting value $N \ge 1$, the sequence will always reach 1.

## Features

-   **Arbitrary Precision:** Uses `MPInt` for all calculations, allowing $N$ to grow beyond standard integer limits.
-   **Step Counting:** Calculates the number of steps to reach 1 for each $N$.
-   **Periodic Recaps:** Every time $N$ reaches a power of 10, the program prints the number with the maximum steps found in the previous decade.
-   **Performance-Conscious:** Uses an asynchronous message-passing loop to prevent blocking the Pony runtime while maintaining responsiveness to `Ctrl+C`.

## Usage

Build the example using the project's `Makefile`:

```bash
make build-examples
```

Run the binary:

```bash
./build/debug/collatz_mpint
```

To stop the calculation, use `Ctrl+C`.
