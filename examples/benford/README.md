# Benford's Law Example

This example demonstrates the verification of [Benford's Law](https://en.wikipedia.org/wiki/Benford%27s_law) on various mathematical sequences using the `mathx` library's `MPInt` and `MPFloat` classes.

Benford's Law, also known as the Newcomb-Benford Law, observes that in many naturally occurring sets of numerical data, the leading digit is likely to be small. For example, the number 1 appears as the leading digit about 30% of the time, while 9 appears as the leading digit less than 5% of the time.

See: [The mathematics of Benford's law: a primer ](https://digitalcommons.calpoly.edu/cgi/viewcontent.cgi?article=1090&context=rgp_rsr)

## Sequences Analyzed

1.  **Fibonacci Numbers:** $F_{n} = F_{n-1} + F_{n-2}$, starting with $F_0 = 0, F_1 = 1$.
2.  **Factorials:** $n! = 1 × 2 × ... × n$.
3.  **Linear Recurrence:** $x_{n+1} = 2x_n + 1$, tested with:
    *   5 random starting values.
    *   A specific starting value: $x_0 = 9.94962308959395941219332124109226$.

## How it Works

For each sequence, the program calculates 1,000 terms. It then extracts the first five significant digits of each term and calculates the frequency of digits 0-9 at each of these five positions. The results are presented in tables, comparing the observed frequencies against the theoretical values predicted by Benford's Law.

## Running the Example

To build the example, run:

```bash
make build-examples
```

To run the example:

```bash
./build/debug/benford
```
