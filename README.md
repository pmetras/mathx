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
[x] Complete implementation of `from[A: ((Number | MPInt | MPFloat) & Real[A] val)]` for better types support.
[x] Add `from_u64_array`, `from_u32_array` and `from_u8_array` to `BitMap` or evaluate if we can write `BitMap.from_array[A: (U8 | U16 | U32 | U64)]`. In `MPInt`, replace calls to `BitMap.from_u16_array` with `BitMap.from_u32_array`. Of course, do the same for `to_u16_array`.
[x] Correct the `raw_digits` docstring.
[ ] List private methods that are never used. List public methods too. Do we have synonyms? Can we reduce the API scope?
[ ] Replace `from_mpfloat` into `from[A: ((Number | MPInt | MPFloat) & Real[A] val)]` when `MPFloat` implements `Real[MPFloat]`.

### `MPFloat`

[x] Add support for rounding to `MPFloat`.
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
[ ] Evaluate if we rename `raw_digits` and `???` to `bits` and `from_bits` to be compatible with `F64` or `F32`? The signature is different...
[ ] Write extensive tests for `from_string` to check coverage.
[ ] Check all constructors to see if precision is defined correctly and compatible with GMPFR: If there is an overflow interpreting the value and it can't be coded with the desired precision, the number becomes `-inf` or `inf`.
[ ] Optimize working directly on the `_digits` array instead of creating a new `MPFloat` at each operation.
[ ] Evaluate if worth replacing some `while ... end` loops by `for i in Range...` loops, for easier readability.
[ ] Correct `almost_eq` function to use `MPFloat` instead of `F64`.
[ ] See if pi calculation by Kudnovsky calculation can be optimized by using `MPInt`.
[ ] Test that digit separators `_` are accepted in `from_string`.
[ ] Option to add digit separators `_` in `exact_string`.
[x] See if `i128`, `i64`, `i32`, `i16` and `i8` can be optimized to prevent creating a new `MPFloat`, using `I128.from_bits`.
[ ] Write tests for all `F64` or `F32` corner cases, like cancellation, loss of precision, etc.

[ ] Make all `from_*` constructors private and replace with a generic `from[A]` constructor where `A` is `Number | MPInt | MPFloat`
[ ] Change the default precision to 128 to remain consistent with `from[A]`.


### `MPFRep`

[x] Do we really want to keep `_MPFBase`? No! `MPFRep` already has all methods (but 2) of `_MPFBase`.
[ ] Why is `MPFRef` class defined as `val`? Wasn't the goal to have a `ref` class that can be used to implement in-place operations and have `MPFloat` that provides the `val` interface for the developer?
This is the follow-up step of the plan: optimizations. A class ref MPFRep would mean:
* `_digits` becomes `Array[U8] ref` (mutable in place)
* Operations like `_trunc`, `_add_mag` could mutate instead of allocate
* `MPFloat` would still be `class val`, but it would own an `iso` or `val` snapshot of a `ref MPFRep` underneath
[ ] Do we introduce `zero` and `one` constructors/constants?
[ ] Remove the `\do_not_use\` methods.
[ ] Check the transformations between `MPFRep` and `MPInt`.
[ ] Path to exact rounding...


### `MPFContext`

[x] What are the rules to decides which operations are into `MPFContext` vs `_MPFAlgo`? What is the reasons why are `sqrt`, `ln` or `log`, `exp` in `MPFContext`? What are the advantages of having these operations into `MPFContext`? Wouldn't it be simpler to have a lighweight `MPFContext` and all code into `_MPFAlgo`, saving a delegating call to `MPFContext.op(params)`, replaced by a direct call to `_MPFAlgo.op(context, params)`? Also, the goal of `_MPFAlgo` is to be able to add new operations algorithms in a single place, without having to touch `MPFContext`.


### `_MPFAlgo`

[ ] Instead of creating GMPF classes from MPF with the same name, what about having a `_GMPFAlgo` that is a mapping to MPF and that can be substituted to `_MPFAlgo`?
[x] Treatment of special values must be consistent. `_mul` must do it like `_add`. What is the strategy? Helper methods `_*` are optimized and don't check for special values but only public-facing methods? Explain the rule that must consistently be applied. Check all methods.
[x] Check that helper methods are used consistently. They must not be inlined in the code. Also, there must not be duplicate helper methods.
[ ] Correct Chudnovsky algorithm for pi. 


### Anomalies

[ ] Review the type hierarchy - See Zim.
[x] What the use of `FloatingPoint` functions `ldexp` and `frexp` that seems not to use the value of `this` but the paramter `x`?
[ ] Boolean function `nan`, `finite` and `infinite` should be renamed `is_nan`, `is_finite` an `is_infinite` to stay consistent with usage.
[ ] Add boolean function `is_integer` to `FloatingPoint`.
[ ] Why are `SignedInteger.bitwidth` and `SignedInteger.bytewidth` returning different types (respectively `A` and `USize`)?
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


# Question

[ ] Note that there are presently 10 failing tests.
[ ] Pony does not support overload. `powi` signature must be changed.
[ ] Should the precision be a compile-time parameter or a run-time parameter?
[ ] Is there any advantage at having the precision a run-time parameter? In which situation would it be used?
[ ] Are the guard constants per operation dependant of magnitude or precision of the operands?
[ ] Tests will need to be adapted.


[ ] `FloatingPoint` defines `precision2`, `min_exp2`, etc. as tag methods (callable without an instance), implying fixed values for all instances. `MPFloat` is parametric in precision — different instances have different precisions. This is a fundamental mismatch.
[ ] What is `FloatingPoint.radix`? What is the supposed `U8`returned value?
[ ] How to implement `FloatingPoint` to unsigned values, according to rounding mode?
[ ] `tag` methods like `radix`, `precision2`, `precision10`, `min_exp2`, `min_exp10`, `max_exp2`, `max_exp10` have meaning only for fixed precision (width) floats.
1. add_unsafe, sub_unsafe, mul_unsafe, div_unsafe, neg_unsafe, eq_unsafe, etc. — default implementations use +~, -~, *~, /~
These all use Pony's built-in wrapping/unsafe arithmetic operators (+~, ==~, etc.), which are compiler intrinsics over machine words. The trait provides them as defaults. MPFloat must override every single one — but it already does, so this is not a blocker, just a note that the defaults are completely wrong and must never be inherited.

2. hash() and hash64() — default implementations in Real[A] call usize() and u64()

fun hash64(): U64 =>
  var x = u64()   // truncates to 64 bits — loses all precision beyond F64 range
  x = (not x) + (x << 21)
  ...
This is a correctness problem: two MPFloat values that differ only in digits beyond the 64th bit would hash to the same value, violating the contract a == b → a.hash() == b.hash() only accidentally. MPFloat must override both hash() and hash64() with an implementation that hashes over all _digits.

3. string() — default implementation in FloatingPoint converts to F64 first

fun string(): String iso^ =>
  _ToString._f64(f64())
This would silently truncate every MPFloat to double precision in string output. MPFloat already overrides string(), so again not a blocker — but if the override is ever removed, this default is a silent precision trap.

4. min_value() and max_value() constructors — inherited from Real[A]
These are required abstract constructors. For a fixed-precision type they mean "the most negative / most positive representable value." For MPFloat, the exponent is an I64, so the actual bounds are enormous but finite. The question is: what precision should min_value() and max_value() use? There is no this to query, so they can only use the default precision. This is semantically odd — MPFloat.min_value() would be a 112-bit float, but a user working at 10000-bit precision would never use it.

5. epsilon() and min_normalised() constructors
epsilon() means the smallest value e such that 1 + e != 1. For a fixed type this is a constant. For MPFloat it depends on the precision of the instance doing the computation — but as a no-arg constructor it has no instance. Same problem as tag methods: must return a nominal value for the default precision.

6. ldexp(x: A, exponent: I32) and frexp(): (A, I32) — exponent clamped to I32
I32 allows exponents up to ±2 billion. MPFloat uses I64 for its exponent — but more importantly, it is in base 256, so the actual binary exponent is _exponent * 8. A value like MPFloat with _exponent = I64.max_value() has a binary exponent of ~7.4 × 10^19, which does not fit in I32. ldexp/frexp would silently corrupt such values. For any realistic computation this is unlikely to matter, but it is a formal incompatibility.

7. from[B: (Number & Real[B] val)](a: B) — generic conversion constructor
This constructor must accept any Number & Real value — U8, I128, F32, etc. — and produce an MPFloat. The precision of the result is unspecified. This requires a decision: use default precision? Use precision derived from B? There is no clean answer from the trait's perspective.

The fundamental conclusion is: FloatingPoint[A] was designed for IEEE 754 machine types (F32, F64), where all these values are compile-time constants. Implementing it on a variable-precision type means accepting that several of its constructors and tag methods are nominal rather than exact. The plan's approach — implement the trait for interoperability, document the limitations — is the right pragmatic position. The more correct long-term answer would be a FloatingPoint redesign that separates the "numeric operations" interface from the "fixed-format metadata" interface, but that is a Pony standard library change, not something this project can drive.

