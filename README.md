# Mathematical library

The library includes various mathematical algorithms.

* `FLimits` class has functions to determines the platform/compiler limits of floating point numbers.
* `Primes` class deals with prime numbers and includes the Miller-Rabin probabilistic primality test.
* `Complex` is a primitive implementing complex numbers arithmetic.
* `FFT` implements the Fast Fourier Transform algorithm.
* `MPInt` and `MPFloat` are multiple precision integer and floating point implementations. A more complete implementation using [GMP](https://gmplib.org/) and [MPFR](https://mpfr.org/) libraries through FFI is also available.

## Setup

You need to install `libgmp-dev` and `libmpfr-dev` packages if you want to use GMP and MPFR

## TODO

### `MPFloat`

[ ] Add support for rounding to `MPFloat`.
[ ] Add support for base in `from_string`.
[ ] Add implementation of number-theoritic transform. Check which one is more performant. Can it work on larger *digit*s?
[x] Keep the base of the digits (256.0) and related information in fields, like what was done in `MPInt`
[x] Remove `clone`? Do we need it?
[ ] Probably need `is_positive`, to demonstrate that `is_nan == not is_positive and not is_negative`?
[x] `from_string` must recognize "   nan   " or "   0.00000   ".
[x] `from_string` must accept numbers larger than `F64`.
[ ] In `from_string`, why operation are done with `F64` instead of `U32`, for instance, and converted to `F64.from_bits` in the end?
[ ] Recognize `_` digit-separator in `from_string`.
[x] Add `from_mpint` to create a new `MPFloat`from an `MPInt`.
[ ] Add `mpint` to convert into a `MPInt`. The type of rounding must be thought about and if the function must be partial...
