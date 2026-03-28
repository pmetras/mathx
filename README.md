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
[ ] Evaluate id we rename `raw_digits` and `???` to `bits` and `from_bits` to be compatible with `F64` or `F32`? The signature is different...
[ ] Write extensive tests for `from_string` to check coverage.
[ ] Check all constructors to see if precision is defined correctly and compatible with GMPFR: If there is an overflow interpreting the value and it can't be coded with the desired precision, the number becomes `-inf` or `inf`.
[ ] Optimize working directly on the `_digits` array instead of creating a new `MPFloat` at each operation.
[ ] Evaluate if worth replacing some `while ... end` loops by `for i in Range...` loops, for easier readability.
[ ] Correct `almost_eq` function to use `MPFloat` instead of `F64`.
[ ] See if pi calculation by Kudnovsky calculation can be optimized by using `MPInt`.
[ ] Test that digit separators `_` are accepted in `from_string`.
[ ] Option to add digit separators `_` in `exact_string`.
[x] See if `i128`, `i64`, `i32`, `i16` and `i8` can be optimized to prevent creating a new `MPFloat`, using `I128.from_bits`.

### Anomalies
[ ] What the use of `FloatingPoint` functions `ldexp` and `frexp` that seems not to use the value of `this` but the paramter `x`?
[ ] Boolean function `nan`, `finite` and `infinite` should be renamed `is_nan`, `is_finite` an `is_infinite` to stay consistent with usage.
[ ] Add boolean function `is_integer` to `FloatingPoint`.
[ ] Why are `SignedIteger.bitwidth` and `SignedInteger.bytewidth` returning different types (respectively `A` and `USize`)?
[ ] In trait `FloatingPoint`, the methods `precision2` and `precision10` return `U8`. That type is not adapted for multiprecision types where precision can be >> 256. Should be `ULong`
[ ] Same comment for `min_exp2`, `min_exp10`, `max_exp2` and `max_exp10` that return `U16`. Should be `ULong`
[ ] `frexp` has a wrong signature in `F32`, `F64` and `FloatingPoint`. It is currently `frexp(): (A, U32)`. It should be `frexp(): (A, ILong)` or `frexp(): (A, I32)` using an integer for the exponent instead of an unsigned. See https://llvm.org/docs/LangRef.html#llvm-frexp-intrinsic