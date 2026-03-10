# Parallel Collatz Conjecture with MPInt

This example parallelizes the investigation of the [Collatz conjecture](https://en.wikipedia.org/wiki/Collatz_conjecture) using Pony's actor model and `MPInt` (Arbitrary-precision integers).

## The Conjecture

The Collatz function is defined for any positive integer $n$ as:
-   If $n$ is even: $n = n / 2$
-   If $n$ is odd: $n = 3n + 1$

The conjecture states that for any starting value $N \ge 1$, the sequence will always reach 1.

## Features

-   **Parallel Execution:** Distributes the calculation of Collatz sequences across multiple worker actors, maximizing CPU utilization on multi-core systems.
-   **Actor-Based Load Balancing:** 
    -   A `Master` actor manages the range of $N$ and coordinates workers.
    -   `Worker` actors perform the intensive `MPInt` calculations independently.
-   **Arbitrary Precision:** Uses `MPInt` to handle sequences that grow beyond standard 64-bit limits.
-   **Progress Tracking:** Periodic recaps every power of 10 show the $N$ with the most steps and the highest value reached globally.
-   **Non-Blocking Design:** Prevents runtime saturation by processing numbers in batches and maintaining responsiveness to system events.

## Usage

Build the example using the project's `Makefile`:

```bash
make build-examples
```

Run the binary:

```bash
./build/debug/collatz_parallel
```

To stop the calculation, use `Ctrl+C`.
