# Mathematical library

The library includes various mathematical algorithms.

* `FLimits` class has functions to determines the platform/compiler limits of floating point numbers.
* `Primes` class deals with prime numbers and includes the Miller-Rabin probabilistic primality test.
* `Complex` is a primitive implementing complex numbers arithmetic.
* `FFT` implements the Fast Fourier Transform algorithm.
* `MPInt` and `MPFloat` are multiple precision integer and floating point implementations. A more complete implementation using [GMP](https://gmplib.org/) and [MPFR](https://mpfr.org/) libraries through FFI is also available.

## Setup

You need to install `libgmp-dev` and `libmpfr-dev` packages if you want to use GMP and MPFR

## References

The `docs/` directory contains copies of the research papers behind this implementation:

- [What Every Computer Scientist Should Know About Floating-Point Arithmetic](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html) (Goldberg)
- [MPFR: A Multiple-Precision Binary Floating-Point Library With Correct Rounding](https://www.mpfr.org/algorithms.pdf) (Fousse, Hanrot, Lefèvre, Pélissier, Zimmermann — INRIA RR-5753)
- [Multidigit Multiplication for Mathematicians](https://cr.yp.to/papers/m3-20010811-retypeset-20220327.pdf) (Bernstein)

## TODO

### `MPInt`
- [ ] List private methods that are never used. List public methods too. Do we have synonyms? Can we reduce the API scope?
- [ ] Replace `from_mpfloat` into `from[A: ((Number | MPInt | MPFloat) & Real[A] val)]` when `MPFloat` implements `Real[MPFloat]`. This will require reviewing the math types hierarchy...

### `MPFloat`

- [ ] Add support for base in `from_string`.
- [ ] Add implementation of number-theoritic transform. Check which one is more performant. Can it work on larger *digit*s?
- [ ] Probably need `is_positive`, to demonstrate that `is_nan == not is_positive and not is_negative`?
- [ ] Recognize `_` digit-separator in `from_string`.
- [ ] Add `mpint` to convert into a `MPInt`. The type of rounding must be thought about and if the function must be partial...
- [ ] Evaluate if we rename `raw_digits` and `???` to `bits` and `from_bits` to be compatible with `F64` or `F32`? The signature is different...
- [ ] Write extensive tests for `from_string` to check coverage.
- [ ] Check all constructors to see if precision is defined correctly and compatible with GMPFR: If there is an overflow interpreting the value and it can't be coded with the desired precision, the number becomes `-inf` or `inf`.
- [ ] Correct `almost_eq` function to use `MPFloat` instead of `F64`.
- [ ] Test that digit separators `_` are accepted in `from_string`.
- [ ] Option to add digit separators `_` in `exact_string`.
- [ ] Write tests for all `F64` or `F32` corner cases, like cancellation, loss of precision, etc.

### `MPFRep`

- [ ] Do we introduce `zero` and `one` constructors/constants?
- [ ] Remove the `\do_not_use\` methods.
- [ ] Check the transformations between `MPFRep` and `MPInt`.

### `_MPFAlgo`

- [ ] Instead of creating GMPF classes from MPF with the same name, what about having a `_GMPFAlgo` that is a mapping to MPF and that can be substituted to `_MPFAlgo`?

### Anomalies

- [ ] Review the type hierarchy - See Zim.
- [ ] Boolean function `nan`, `finite` and `infinite` should be renamed `is_nan`, `is_finite` an `is_infinite` to stay consistent with usage.
- [ ] Add boolean function `is_integer` to `FloatingPoint`.
- [ ] Why are `SignedInteger.bitwidth` and `SignedInteger.bytewidth` returning different types (respectively `A` and `USize`)?
- [ ] In trait `FloatingPoint`, the methods `precision2` and `precision10` return `U8`. That type is not adapted for multiprecision types where precision can be >> 256. Should be `ULong`
- [ ] Same comment for `min_exp2`, `min_exp10`, `max_exp2` and `max_exp10` that return `U16`. Should be `ULong`

# Miscellaneous

- [ ] `powi` signature must be changed.
- [ ] `tag` methods like `radix`, `precision2`, `precision10`, `min_exp2`, `min_exp10`, `max_exp2`, `max_exp10` have meaning only for fixed precision (width) floats.
