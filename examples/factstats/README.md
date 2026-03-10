# Factorial Statistics Example

This example demonstrates the use of `MPInt` (Arbitrary-precision integers) to calculate factorials of numbers from 0 up to 1000 and analyze the distribution of digits in their decimal representation.

## Features

- **Arbitrary Precision:** Uses `MPInt` to calculate $N!$ for $N \in [0, 1000]$.
- **Digit Analysis:** Counts occurrences of each decimal digit (0-9) in each factorial result.
- **Statistical Summary:** 
  - Provides cumulative sums and relative percentages of each digit.
  - Displays intermediate summaries and re-prints the header every 50 lines for readability.
- **Formatted Output:** 
  - Truncates large factorial results to the 120 most significant digits to keep the table readable.
  - Displays the total digit count for each factorial.

## Usage

Build the example using the project's `Makefile`:

```bash
make build-examples
```

Run the binary:

```bash
./build/debug/factstats
```

## Implementation Details

The program uses `MPInt` for exact integer arithmetic. Factorials grow very quickly; for instance, $1000!$ has 2568 digits. The `MPInt.string()` method converts these large numbers into decimal strings for analysis. The output is formatted into a fixed-width table for easy inspection of how digits are distributed as $N$ increases.
