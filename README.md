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

### `MPInt`
[x] `karatsuba_mul` USize underflow crash: when both operands are > 128 digits but one is less than half the size of the other (e.g. 140-digit × 300-digit), the capacity hint `that._digits.size() - half` wraps around on USize subtraction → OOM crash. Fix: use clamped subtraction in lines 855 and 869.
[x] `from_string` exponent builds a literal zero-padded string of length `exp`, then processes it character-by-character with `_short_mul`. For large exponents (e.g. `1@1000000`) this allocates O(exp) memory and runs O(exp²) operations. Should use `digit_shl` + `pow(base, exp)` instead.
[x] `karatsuba_mul` threshold uses `or` instead of `and` (line 842): falls back to schoolbook if *either* operand is ≤ 128 digits, contradicting the docstring ("both numbers must be > 2048-bits"). Also masks the underflow bug above for some size combinations.
[x] `fast_mul`: no bounds check that inputs are within the ~1M-digit safe range for F64 precision; silently returns wrong results for larger inputs. Also discards any remaining carry after the digit-conversion loop instead of asserting it is zero.
[x] Two private constructors with inconsistent zero-normalization: `_create` normalizes `-0 → +0`; `_from_array` does not, allowing an invalid `_negative=true, _digits=[0]` state.
[x] `divrem` single-digit path (lines 900–921) copies `_digits` twice and calls `_short_div` twice to obtain quotient and remainder independently; a single call returns both.
[x] `is_zero()` contains an unreachable `_digits.size() == 0` guard (line 2004); the invariant guarantees `_digits.size() ≥ 1`.
[x] `isqrt()` silently incorrect for numbers > ~10^19 digits: `bitwidth().usize()` overflows via `ilong()`, corrupting the Newton initial guess.
[x] Make NTT really generic and not assume that `A` is `USize` or platform-dependant.
[ ] Complete implementation of `from[A: ((Number | MPInt | MPFloat) & Real[A] val)]` for better types support.
[x] Add `from_u64_array`, `from_u32_array` and `from_u8_array` to `BitMap` or evaluate if we can write `BitMap.from_array[A: (U8 | U16 | U32 | U64)]`. In `MPInt`, replace calls to `BitMap.from_u16_array` with `BitMap.from_u32_array`. Of course, do the same for `to_u16_array`.
[x] Correct the `raw_digits` docstring.
[ ] List private methods that are never used. List public methods too. Do we have synonyms? Can we reduce the API scope?

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
[x] What the use of `FloatingPoint` functions `ldexp` and `frexp` that seems not to use the value of `this` but the paramter `x`?
[ ] Boolean function `nan`, `finite` and `infinite` should be renamed `is_nan`, `is_finite` an `is_infinite` to stay consistent with usage.
[ ] Add boolean function `is_integer` to `FloatingPoint`.
[ ] Why are `SignedIteger.bitwidth` and `SignedInteger.bytewidth` returning different types (respectively `A` and `USize`)?
[ ] In trait `FloatingPoint`, the methods `precision2` and `precision10` return `U8`. That type is not adapted for multiprecision types where precision can be >> 256. Should be `ULong`
[ ] Same comment for `min_exp2`, `min_exp10`, `max_exp2` and `max_exp10` that return `U16`. Should be `ULong`
[x] `frexp` has a wrong signature in `F32`, `F64` and `FloatingPoint`. It is currently `frexp(): (A, U32)`. It should be `frexp(): (A, ILong)` or `frexp(): (A, I32)` using an integer for the exponent instead of an unsigned. See https://llvm.org/docs/LangRef.html#llvm-frexp-intrinsic
[x] `F32` and `F64` `min` implementation is not IEEE 754 compliant as it does not propagate NaN values:
    * NaN.min(5) → returns 5
    * 5.min(NaN) → returns NaN
    * NaN.min(NaN) → returns NaN
    We should have:
    * NaN.min(5) → returns NaN
    * 5.min(NaN) → returns NaN
    * NaN.min(NaN) → returns NaN