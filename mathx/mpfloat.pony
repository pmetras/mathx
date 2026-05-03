// Multi-precision floating point numbers — thin public wrapper.
//
// `MPFloat` is the stable public API.  Internally it bundles two orthogonal
// components:
//
//   - `_rep: MPFRep`      — the pure representation (sign, exponent, digits)
//   - `_ctx: MPFContext`  — the precision and rounding policy (no arithmetic)
//
// All arithmetic calls `_MPFAlgo.op(_ctx, ...)` directly.  `_MPFAlgo`
// handles special values, runs the algorithm at working precision, and rounds
// the result to output precision via `_ctx._round_to`.

use "../assertx"
use "../formatx"

use "collections"
use "debug"


class val MPFloat is (Formattable & Stringable)
  """
  MPFloat represents real numbers with arbitrary precision.

  https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic

  Contrarily to fixed-precision floating-point numbers, MPFloat supports
  unlimited precision limited only by available memory and the needs of the
  user.

  Numbers are stored in normalised form:

    value = (-1)^_sign × 0.d₀d₁d₂… × 256^_exponent

  where d₀…d_{n-1} are base-256 digits (U8) stored in big-endian order in
  `_digits` (most-significant digit at index 0), and `_exponent` is a signed
  base-256 exponent. For a non-zero finite number `d₀ ≠ 0` (i.e. the
  representation is normalised). Zero has an all-zero digit array.

  Special values: NaN (`_nan = true`) and ±infinity (`_inf = true`) are
  supported, mirroring the GMP/MPFR binding in `mathx/gmp/gmpfloat.pony`.

  `MPFloat` is a `class val` — globally immutable, just like `MPInt`. All
  arithmetic operations return new `MPFloat` instances.

  The multiplication of two `MPFloat` uses a fast-multiplication algorithm
  using real Fast Fourier Transform, allowing floats with around 10e6 decimal
  digits. For higher precision, Number-Theoretic Transforms would be needed.

  Both `from_string` and `string`/`exact_string` use full multi-precision
  arithmetic; precision is limited only by `prec` bits, not by F64.

  The precision `prec` is the number of bits for the mantissa (significand).
  The precision of a `MPFloat` is set when the instance is created. It is
  expected that all `MPFloat`s involved in operations use the same
  precision. If that's not the case, warning messages are printed when compiled
  in debug mode, and results can be not the one expected by the developer. As a
  consequence, the precision of the result of the operation `this op that` or
  `this.op(that, ...)` is set to the precision of `this`.

  The default precision is 128 bits for the mantissa, giving ~38 significant
  decimal digits. This covers the full value range of `U128` and `I128`
  without truncation, which is why 128 bits is preferred over the 112-bit
  mantissa of IEEE 754 binary128.

  Usage of `MPFloat` class assumes that floats have the same precision and it
  prints warnings, when compiled in debug mode, on operations involving
  instances with different precisions. This is because accuracy can degrade
  when mixing operands with different precisions. But nothing prevent you from
  doing it, and there's even a constructor to create a new `MPFloat` from
  another one and changing its precision.

  ## Rounding guarantees

  Every operation rounds its result to the output precision under
  `rounding_mode()`. Two levels of guarantee are provided:

  | Operation                        | Guarantee           | Notes                       |
  | -------------------------------- | ------------------- | --------------------------- |
  | `add`, `sub`, `mul`              | Correct (≤ 0.5 ULP) | IEEE 754 §5.4: exact then round once  |
  | `div`, `inv`, `sqrt`             | Correct (≤ 0.5 ULP) | Newton + one refinement step          |
  | `ln`, `log`, `log2`, `log10`     | Correct (≤ 0.5 ULP) | Ziv's iteration                      |
  | `exp`, `exp2`                    | Correct (≤ 0.5 ULP) | Ziv's iteration                      |
  | `powi`, `pow`                    | Correct (≤ 0.5 ULP) | Ziv's iteration                      |
  | `logb`                           | Faithful (≤ 1 ULP)  | `ln/ln` — two roundings              |
  | `sin`, `cos`, `tan`              | Faithful (≤ 1 ULP)  | argument in radians                  |
  | `csc`, `sec`, `cot`              | Faithful (≤ 1 ULP)  | derived from sin/cos                 |
  | `sinh`, `cosh`, `tanh`           | Faithful (≤ 1 ULP)  |                                      |
  | `sech`, `csch`, `coth`           | Faithful (≤ 1 ULP)  |                                      |

  **Correctly rounded** — the result is bit-for-bit identical to rounding the
  infinite-precision mathematical value under `rounding_mode()`. This matches
  the guarantee that MPFR provides. For `add`/`sub`/`mul` this follows from
  IEEE 754 §5.4: the intermediate result is exact and a single rounding step
  is applied. For transcendentals it is achieved via Ziv's iteration: working
  precision is doubled and the computation retried whenever the rounding
  direction cannot be determined unambiguously from the guard bits.

  **Faithfully rounded** — the result differs from the exact value by at most
  1 ULP in either direction. This is the level provided by all trigonometric
  and hyperbolic functions (`sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`, and
  their reciprocals), `logb`, and `pi`.
  """

  let _rep: MPFRep val
  """
  The pure representation: sign, exponent, and digit array.
  All structural queries are delegated here.
  """

  let _ctx: MPFContext
  """
  The precision and rounding policy.  Passed as first argument to `_MPFAlgo`
  operations; not an arithmetic dispatcher.
  """


  new val _from(r: MPFRep val, c: MPFContext) =>
    """
    Private canonical constructor used by all arithmetic results.
    Callers supply an already-rounded `MPFRep` from `_MPFAlgo` and the context.
    """
    _rep = r
    _ctx = c


  new val _create(
    sgn: Bool,
    nan: Bool,
    inf: Bool,
    expn: I64,
    digits: Array[U8] val,
    rnd: RoundingMode = RoundingNearest)
  =>
    """
    Private constructor retained for internal use that still needs a raw
    MPFRep + explicit fields (e.g. transcendental constructors that build the
    final value from intermediate results). The precision is inferred from
    the digit array length.
    """
    _rep = MPFRep._create(sgn, nan, inf, expn, digits)
    _ctx = MPFContext(digits.size() * 8, rnd)


  new val create(prec: USize = 128, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a positive zero with `prec` bits of mantissa (significand).

    This is the default constructor: `MPFloat()` gives a 128-bit positive zero,
    and `MPFloat(n)` gives a size-n bits positive zero (all digits zero).
    """
    let p = (prec + 7) / 8
    _rep = MPFRep.create(p)
    _ctx = MPFContext(prec, rnd)


  new val nan_val() =>
    """
    Create a Not-a-Number value.
    """
    _rep = MPFRep.nan_val()
    _ctx = MPFContext(128, RoundingNearest)


  new val inf_val(positive: Bool = true) =>
    """
    Create an infinite value. Pass `positive = false` for −∞.
    """
    _rep = MPFRep.inf_val(positive)
    _ctx = MPFContext(128, RoundingNearest)


  new val from_string(s: String = "",
                      prec: USize = 128,
                      base: U8 = 10,
                      rnd: RoundingMode = RoundingNearest) ? =>
    """
    Create a new `MPFloat` by parsing the string `s` with `prec` bits
    of precision (default 128, ≈38 significant decimal digits). Raises
    an error if `s` is not a recognised floating-point representation.

    The `base` parameter (default 10) selects the numeral base; currently
    only base 10 is implemented. The `rnd` parameter (default nearest) is
    accepted for API compatibility with GMP; rounding support is a future TODO.

    Accepted formats:
    - `""`, `"0"`, `"0.0"`, `"+0"`, `"+0.0"` → +0
    - `"-0"`, `"-0.0"` → −0
    - `"nan"`, `"NaN"`, `"@NaN@"` → NaN
    - `"+inf"`, `"@Inf@"`, `"inf"` → +∞
    - `"-inf"`, `"-@Inf@"` → −∞
    - `[+-]d+[.d*][e|E|@[+-]d+]` → finite decimal

    The exponent part is always in decimal. For base larger than 10, the exponent
    character is `'@'`. For instance, you must write `1.234@56` instead of
    `1.234e56` if you are using base 16. For bases lower or equal to 10, one can
    use `e`, `E` or `@` to indicate the exponent part of the number.

    Leading and trailing whitespace is stripped. The character `_` can be used
    to separate groups of digits. Unlike the `F64` or `F32` paths, the
    decimal mantissa is parsed directly into multi-precision arithmetic, so
    the full `prec` bits of precision are exploited regardless of how many
    significant digits the string contains. For special values NaN and infinites,
    case is significative.
    """
    if base != 10 then error end  // TODO: implement parsing for bases other than 10
    let c = MPFContext(prec, rnd)
    _rep = _MPFAlgo.from_string(c, s)?
    _ctx = c


  new val from[A: ((Number | MPInt | MPFloat) & Real[A] val)](n: A,
                                                             prec: USize = 128,
                                                             rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` from any numeric value `n` with `prec` bits of
    precision (default 128, ~38 significant decimal digits) and rounding mode
    `rnd` (default nearest).

    Accepts:
    - `MPFloat` — precision-change copy
    - `MPInt` — exact integer conversion
    - `F64` — preserves fractional part and special values (NaN, ±∞)
    - `F32` — preserves fractional part and special values (NaN, ±∞)
    - `U8`, `U16`, `U32`, `U64`, `U128`, `ULong`, `USize` — exact integer
    - `I8`, `I16`, `I32`, `I64`, `I128`, `ILong`, `ISize` — exact integer

    The 128-bit default covers the full range of `U128`/`I128` without
    truncation. Integer paths are exact up to `prec` bits; float paths are
    rounded to the nearest representable value.
    """
    let p = (prec + 7) / 8
    _ctx = MPFContext(prec, rnd)
    iftype A <: MPFloat then
      let f: MPFloat = n
      _rep = MPFRep.from_mpfloat_rep(f._rep, p)
    elseif A <: F64 then
      _rep = MPFRep.from_f64(n.f64(), p)
    elseif A <: F32 then
      _rep = MPFRep.from_f32(n.f32(), p)
    elseif A <: MPInt then
      let m: MPInt = n
      _rep = MPFRep.from_mpint(m, p)
    elseif A <: (I8 | I16 | I32 | I64 | I128 | ILong | ISize) then
      // Signed integer types: convert via I128 (exact for all).
      _rep = MPFRep.from_mpint(MPInt.from[I128](n.i128()), p)
    elseif A <: (U8 | U16 | U32 | U64 | U128 | ULong | USize) then
      // Unsigned integer types (U8..U64, ULong, USize, U128): convert via U128.
      _rep = MPFRep.from_mpint(MPInt.from[U128](n.u128()), p)
    else
      Debug("[MPFloat.from] Unknown new type. Constructor failed! Adapt MPFloat code.")
      _rep = MPFRep
    end


  new \do_not_use\ min_normalized(prec: USize = 128, rnd: RoundingMode = RoundingNearest) =>
    """
    The smallest normalized floating point number.

    As don't have the notion of normalized number in the sense of the binary
    floating point representation, the smallest normalized `MPFloat` has been
    set to `0.0`.
    """
    let p = (prec + 7) / 8
    _rep = MPFRep.min_normalized(p)
    _ctx = MPFContext(prec, rnd)


  new epsilon(prec: USize = 128, rnd: RoundingMode = RoundingNearest) =>
    """
    Create the *machine epsilon* for the given precision: the smallest positive
    `MPFloat` ε such that `1 + ε ≠ 1` in base 256, i.e. `ε = 256^(1 − p_digits)`.

    For the default precision of 128 bits (16 bytes, ≈ 38 significant decimal digits):
    `ε ≈ 1.47 × 10^{−38}`.

    This is the base-256 analogue of `F64.epsilon() ≈ 2.22 × 10^{−16}`.
    In general `MPFloat.epsilon(p)` decreases as `p` grows — the user controls
    precision and thus controls the machine epsilon.

    Example use in tests: compare with absolute tolerance `2 × epsilon`:
    ```
    let eps = MPFloat.epsilon(p)
    h.assert_true(got.almost_eq(expected, MPFloat.from[F64](1e-31), MPFloat.from[F64](1e-31)))
    ```
    """
    let p = (prec + 7) / 8
    _rep = MPFRep.epsilon(p)
    _ctx = MPFContext(prec, rnd)


  new val min_value() =>
    """
    The minimum value is −∞ (`-inf`).
    """
    _rep = MPFRep.min_value()
    _ctx = MPFContext(128, RoundingNearest)


  new val max_value() =>
    """
    The maximum value is +∞ (`+inf`).
    """
    _rep = MPFRep.max_value()
    _ctx = MPFContext(128, RoundingNearest)


  new val pi(prec: USize = 128, rnd: RoundingMode = RoundingNearest) =>
    """
    The pi constant, calculated with the specified `prec` accuracy (number of
    bits, default 128 giving ~38 decimal digits). The rounding
    mode `rnd` is not used yet (TODO).

    Uses Machin's formula:
      π = 16·arctan(1/5) − 4·arctan(1/239)
    where each arctan is computed by the Taylor series
      arctan(x) = x − x³/3 + x⁵/5 − x⁷/7 + ⋯
    Both arguments are small (< 0.25), giving fast geometric convergence.
    """
    let c = MPFContext(prec, rnd)
    _rep = _MPFAlgo.pi(c)
    _ctx = c


  new val pi_chudnovsky(prec: USize = 128, rnd: RoundingMode = RoundingNearest) =>
    """
    Compute π using the Chudnovsky algorithm at the given precision.

    TODO: Fix bug
    ⚠ Warning: this implementation currently gives only ~4 correct digits —
    see `_MPFAlgo._pi_chudnovsky` for the known bug description.
    """
    let c = MPFContext(prec, rnd)
    _rep = _MPFAlgo.pi_chudnovsky(c)
    _ctx = c


  new val pi_bbp(prec: USize = 128, rnd: RoundingMode = RoundingNearest) =>
    """
    Compute π using the Bailey–Borwein–Plouffe (BBP) formula:
    `π = Σ_{k=0}^∞ (1/16^k) × [4/(8k+1) − 2/(8k+4) − 1/(8k+5) − 1/(8k+6)]`.
    """
    let c = MPFContext(prec, rnd)
    _rep = _MPFAlgo.pi_bbp(c)
    _ctx = c


  //- Internal helpers ---------------------------------------------------------

  fun get_precision(): USize =>
    """
    Get the precision of the `MPFloat` mantissa (significand) in bits.

    Note: This method is kept while the `precision2` return type is not corrected in
    stdlib, for compatibility with GMP implementation.
    """
    _ctx.precision


  fun get_rounding_mode(): RoundingMode =>
    """
    Get the rounding mode that was used when initializing the `MPFloat`.
    """
    _ctx.rounding


  fun _size(): USize =>
    """
    Number of base-256 digits in the internal representation.
    """
    _rep._size()


  fun _trunc(n: USize): MPFloat =>
    """
    Return a new `MPFloat` containing only the `n` most-significant base-256
    digits of `this`. Used to bound intermediate results in Newton iterations
    and series computations so that successive multiplications do not grow the
    array without limit.

    This is a structural truncation (no rounding); the `_ctx` is preserved.
    """
    _from(_rep._trunc(n), _ctx)


  //- Predicates ---------------------------------------------------------------

  fun is_nan(): Bool =>
    """
    `true` if and only if `this` is Not-a-Number.
    """
    _rep.is_nan()


  fun is_infinite(): Bool =>
    """
    `true` if and only if `this` is ±∞.
    """
    _rep.is_infinite()


  fun is_finite(): Bool =>
    """
    `true` if and only if `this` is a finite number (not NaN and not ±∞).
    """
    _rep.is_finite()


  fun is_zero(): Bool =>
    """
    `true` if and only if `this` is ±0 (both +0 and −0 are zero).
    """
    _rep.is_zero()


  fun is_negative(): Bool =>
    """
    `true` if `this` is strictly negative (sign bit set and value ≠ 0).
    NaN → false.
    """
    _rep.is_negative()


  fun is_integer(): Bool =>
    """
    `true` if `this` is a finite integer (no fractional part).
    NaN, ±∞ → false.
    """
    _rep.is_integer()


  fun almost_eq(that: box->MPFloat,
                rel_tol: MPFloat = MPFloat.epsilon().sqrt(),
                abs_tol: MPFloat = MPFloat.epsilon().sqrt())
               : Bool =>
    """
    Return `true` when `this` and `that` are approximately equal.

    Two `MPFloat` values are almost equal when

      `|this − that| ≤ max(rel_tol × max(|this|, |that|), abs_tol)`

    Both tolerance parameters default to `sqrt(MPFloat.epsilon())`,
    using `MPFloat`'s default precision and rounding'. If you use
    `MPFloat` with different precision or round mode, you should
    give values to `rel_tol` and `abs_tol` tuned for your needs.

    The comparison is performed entirely in `MPFloat` arithmetic.

    Special cases (matching IEEE 754 conventions):
    - Either operand is `NaN` → `false`.
    - Both operands are the same infinity (same sign) → `true`.
    - One is infinite, the other finite → `false`.

    This is the MPFloat analogue of `Complex.almost_eq`.
    """
    if is_nan() or that.is_nan() then
      return false
    end
    if is_infinite() then
      return that.is_infinite() and (_rep.sign_bit() == that._rep.sign_bit())
    end
    if that.is_infinite() then
      return false
    end

    // Warn when operands have different precisions: the comparison is valid but
    // the result is only as accurate as the lower-precision value.
    ifdef debug then
      try
        Assert(_size() == that._size(),
              "[MPFloat.almost_eq] Precision mismatch: `this` has " + _size().string() +
              " bytes, `that` has " + that._size().string() +
              " bytes — comparison accuracy limited to the smaller precision")?
      end
    end

    let p: USize = _size().max(that._size())
    let diff = (this - that).abs()
    let abs_this = this.abs()
    let abs_that = that.abs()
    let max_mag = if abs_this < abs_that then abs_that else abs_this end
    let p_bits: USize = p * 8
    let ctx_p = MPFContext(p_bits, _ctx.rounding)
    let rel_part = (_from(MPFRep.from_mpfloat_rep(rel_tol._rep, ctx_p.p_bytes()), ctx_p) * max_mag)
    let abs_part = _from(MPFRep.from_mpfloat_rep(abs_tol._rep, ctx_p.p_bytes()), ctx_p)
    let threshold = if rel_part < abs_part then abs_part else rel_part end
    diff <= threshold


  fun sign(): Compare =>
    """
    Return the sign of this value as a `Compare`:
    - `Less` for strictly negative values
    - `Greater` for strictly positive values
    - `Equal` for zero and NaN (conventional)
    """
    if _rep.is_nan() then
      return Equal
    end
    if is_zero() then
      return Equal
    end
    if _rep.sign_bit() then
      Less
    else
      Greater
    end


  fun exponent(): I64 =>
    """
    Return the base-256 exponent `e` such that `|value| = 0.d₀d₁… × 256^e`.

    The integer part of the value occupies the first `e` bytes of
    `raw_digits()` (when `e > 0`); bytes at position `e` and beyond represent
    the fractional part.  A non-positive exponent means the value is strictly
    in the open interval `(−1, 1)`.
    """
    _rep.exponent()


  fun raw_digits(): Array[U8] val =>
    """
    Return the internal base-256 mantissa as a big-endian `Array[U8]`.

    Together with `exponent()` this determines the magnitude exactly:
    `|value| = 0.d₀d₁… × 256^exponent()`.  The first byte `d₀` is always
    non-zero for finite non-zero values.

    This accessor is intended for low-level conversions such as
    `MPInt.from_mpfloat`; prefer the high-level API for ordinary use.
    """
    _rep.raw_digits()


  fun sign_bit(): Bool =>
    """
    Return `true` when this value is strictly negative or is −0.
    NaN has an undefined sign; this method returns the raw `_sign` field.
    """
    _rep.sign_bit()


  fun rep(): MPFRep val =>
    """
    Return the pure representation of this value as an `MPFRep`.

    The `MPFRep` contains the same sign, exponent, and digit array as `this`,
    without the rounding mode.  Intended for interoperability with `MPFContext`
    and `_MPFAlgo`.
    """
    _rep


  fun ctx(): MPFContext =>
    """
    Return the execution context (precision + rounding mode) of this value
    as an `MPFContext`.
    """
    _ctx


  //- Arithmetic ---------------------------------------------------------------

  fun add(that: MPFloat): MPFloat =>
    """
    Calculate `this + that` with full sign and exponent handling.

    Special value rules:
    - NaN propagates: `NaN + x = NaN`.
    - `+∞ + (−∞) = NaN`.
    - `±∞ + finite = ±∞`.
    - Zero is the additive identity.

    When the operand signs agree, magnitudes are added (`_add_mag`). When
    they differ, the smaller magnitude is subtracted from the larger, and the
    result takes the sign of the larger.

    **Rounding**: correctly rounded — IEEE 754 §5.4. The intermediate sum is
    computed exactly (no cancellation loss at working precision), and a single
    `_round_to` step produces the nearest representable value under
    `rounding_mode()`.
    """
    ifdef debug then
      if _ctx.precision != that._ctx.precision then
        Debug.out("[MPFloat.add] Precision mismatch: this=" + _ctx.precision.string() +
          " bits, that=" + that._ctx.precision.string() + " bits")
      end
    end
    _from(_MPFAlgo.add(_ctx, _rep, that._rep), _ctx)


  fun sub(that: MPFloat): MPFloat =>
    """
    Calculate `this − that`. Delegates to `add` after negating `that`,
    so all sign, exponent, and special-value handling is inherited.

    **Rounding**: correctly rounded — IEEE 754 §5.4 (see `add`).
    """
    _from(_MPFAlgo.sub(_ctx, _rep, that._rep), _ctx)


  fun mul(that: MPFloat): MPFloat =>
    """
    Multiply `this` by `that`.

    **Rounding**: correctly rounded — IEEE 754 §5.4. The full product is
    computed exactly via FFT convolution (all `p_a + p_b` output digits), and
    a single `_round_to` step rounds to the output precision. Note: for
    operands larger than ~10⁶ digits, F64 rounding in the FFT may introduce
    errors; use NTT-based multiplication for extreme precisions.
    """
    ifdef debug then
      if _ctx.precision != that._ctx.precision then
        Debug.out("[MPFloat.mul] Precision mismatch: this=" + _ctx.precision.string() +
          " bits, that=" + that._ctx.precision.string() + " bits")
      end
    end
    _from(_MPFAlgo.mul(_ctx, _rep, that._rep), _ctx)


  fun inv(): MPFloat =>
    """
    Calculate `1 / this` using Newton's method:  `x_{n+1} = x_n × (2 − G × x_n)`
    where `G = 256 × |this|` (the value of `this` reinterpreted with exponent 1).

    `G ∈ [1, 256)`, so the iteration converges to `y = 1/G ∈ (1/256, 1]`.
    The final result exponent is `(res._exponent + 1) − _exponent`.

    Sign propagation: `1/positive = positive`, `1/negative = negative`.
    A safety limit of `size × 4` iterations prevents non-convergence.

    **Rounding**: correctly rounded — IEEE 754 §5.4. The Newton kernel adds
    one extra refinement step at higher working precision to push the error
    below 0.5 ULP before `_round_to` is applied.
    """
    _from(_MPFAlgo.inv(_ctx, _rep), _ctx)


  fun div(that: MPFloat): MPFloat =>
    """
    Calculate `this / that` as `this × (1 / that)`.

    Special cases follow IEEE 754:
    - NaN in either operand → NaN
    - ±∞ / ±∞ → NaN;  finite / ±∞ → ±0;  ±∞ / finite → ±∞
    - finite / 0 → ±∞ (sign = XOR of operand signs);  0 / 0 → NaN

    **Rounding**: correctly rounded — IEEE 754 §5.4. Computed as
    `this × inv(that)`; `inv`'s refinement step ensures the combined error is
    below 0.5 ULP before `_round_to` is applied.
    """
    ifdef debug then
      if _ctx.precision != that._ctx.precision then
        Debug.out("[MPFloat.div] Precision mismatch: this=" + _ctx.precision.string() +
          " bits, that=" + that._ctx.precision.string() + " bits")
      end
    end
    _from(_MPFAlgo.div(_ctx, _rep, that._rep), _ctx)


  fun sqrt(): MPFloat =>
    """
    Calculate the square root of `this` using Newton's method for the
    reciprocal square root: `y_{n+1} = y_n × (3 − H × y_n²) / 2`.

    The iteration converges quadratically to `1 / √H`; the result is then
    multiplied by `H` to obtain `√H`, and the exponent is adjusted to give
    `√this`.

    Parity split: for odd `_exponent` use `H = MPFloat(exp=1, same digits)`
    so `H ∈ [1, 256)`; for even `_exponent` use `H = MPFloat(exp=2, same
    digits)` so `H ∈ [256, 65536)`.  With this choice `_exponent − h_exp`
    is always even, so the result exponent

      `sqrt_h._exponent + (_exponent − h_exp) / 2`

    is always an integer.

    A safety limit of `size × 4` iterations prevents non-convergence.

    **Rounding**: correctly rounded — IEEE 754 §5.4. The Newton kernel adds
    one extra refinement step at higher working precision to push the error
    below 0.5 ULP before `_round_to` is applied.
    Note: `sqrt(9)` returns `2.9999…` rather than exactly 3 — Newton does not
    converge to exact integers for non-dyadic perfect squares at finite
    precision. Use `almost_eq` when testing square root results.
    """
    _from(_MPFAlgo.sqrt(_ctx, _rep), _ctx)


  fun exact_string(base: U8 = 10, rnd: RoundingMode = RoundingNearest)
                  : (String iso^, I64, Bool) =>
    """
    Return `(mantissa, exponent, inexact)` where:

      `|this| = 0.mantissa_digits × 10^exponent`

    The mantissa string begins with `"-"` when `this` is negative, followed
    by decimal digit characters (no decimal point). The number of digits is
    approximately `⌈_size() × log₁₀(256)⌉ + 2`; trailing zeros are NOT
    stripped. The `exponent` is the number of decimal digits to the left of
    the decimal point (so `exponent = 1` means the value is in `[0.1, 1)`
    after dividing by `10^(exponent-1)`).

    Special values:
    - NaN → `("nan", 0, false)`
    - ±∞ → `("±inf", 0, false)`
    - ±0 → `("0", 0, false)`

    Only base 10 is currently implemented; other bases return `("", 0, true)`.
    The `inexact` flag is always `false` (rounding is a future TODO).
    """
    _rep.exact_string(base)


  fun string(): String iso^ =>
    """
    Convert to a human-readable decimal string.

    The number of decimal digits printed is approximately `_size() × log₁₀(256)`.
    The format is `[-]d.dddde±eee` in scientific notation when the exponent
    is large, and plain decimal otherwise.
    """
    _rep.string()


  fun format(spec: String = ""): String =>
    """
    Format this value according to `spec`.
    
    **Type codes**

    - `'e'`/`'E'`: scientific notation `d.dddde±ee`. Precision `p` is digits
      after the decimal point in the mantissa (total `p + 1` significant digits).
      Default: all stored digits.
    - `'f'`/`'F'`: fixed-point notation. Precision `p` is digits after the
      decimal point. Default: all stored digits.
    - `'g'`/`'G'`: general — uses `'e'` when exponent < −4 or ≥ precision,
      `'f'` otherwise. Trailing zeros are stripped unless `#` is set.
      Default precision: all stored digits.
    - No type code: same as `'g'`.
    - `'d'`, `'b'`, `'o'`, `'x'`/`'X'`: truncate toward zero, then format as
      integer in base 10 / 2 / 8 / 16.
    - `'%'`: multiply by 100, format as `'f'`, append `'%'`.

    For special values NaN and ±∞, the type code selects upper/lower case only.
    Width, fill, alignment, sign, `#`, `z`, and grouping follow the standard
    `FormatSpec` grammar.
    """
    _rep.format(spec)


  //- Unsafe operations -------------------------------------------------------

  fun add_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe addition of `this` and `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    _from(_MPFAlgo.add(_ctx, _rep, that._rep), _ctx)


  fun sub_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe substraction of `that` from `this`.If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    _from(_MPFAlgo.sub(_ctx, _rep, that._rep), _ctx)


  fun mul_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe multiplication of `this` and `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.

    TODO: Complete rounding propagation
    """
    _from(_MPFAlgo.mul(_ctx, _rep, that._rep), _ctx)


  fun _trunc_frac(): MPFloat =>
    """
    Internal alias for `trunc()`.  Called by `divrem`, `fld_unsafe`, and
    `divrem_unsafe` to obtain the truncated-toward-zero integer part without
    going through the public API.
    """
    trunc()


  fun divrem(that: MPFloat): (MPFloat, MPFloat) =>
    """
    Truncated division with remainder: returns `(q, r)` such that
    `this = q × that + r`, where `q = trunc(this / that)` is the integer
    nearest to `this / that` in the direction of zero, and `r` has the same
    sign as `this`.

    Special cases follow IEEE 754 conventions for the quotient; the remainder
    is derived as `r = this − q × that`:
    - NaN in either operand → `(NaN, NaN)`.
    - `this` is ±∞ and `that` is finite → `(±∞, NaN)`.
    - `this` is ±∞ and `that` is ±∞ → `(NaN, NaN)`.
    - `that` is ±∞ and `this` is finite → `(0, this)`.
    - `that` is 0 and `this` is 0 → `(NaN, NaN)`.
    - `that` is 0 and `this` is non-zero → `(±∞, NaN)`.
    - `this` is 0 → `(0, 0)`.

    See also [`rem`](#rem), [`fld`](#fld), [`mod`](#mod).
    """
    ifdef debug then
      if _ctx.precision != that._ctx.precision then
        Debug.out("[MPFloat.divrem] Precision mismatch: this=" + _ctx.precision.string() +
          " bits, that=" + that._ctx.precision.string() + " bits")
      end
    end
    (let q, let r) = _MPFAlgo.divrem(_ctx, _rep, that._rep)
    (_from(q, _ctx), _from(r, _ctx))


  fun rem(that: MPFloat): MPFloat =>
    """
    Truncated remainder of `this / that`: the value `r` satisfying
    `this = trunc(this / that) × that + r`.

    `r` has the same sign as `this` (the dividend).  For integer-valued
    operands this is equivalent to `this mod that` in the C / truncation sense.

    Special cases are inherited from `divrem`.

    See also [`divrem`](#divrem), [`fld`](#fld), [`mod`](#mod).
    """
    divrem(that)._2


  fun fld(that: MPFloat): MPFloat =>
    """
    Floored division of `this` by `that`: `floor(this / that)`.

    The result is the largest integer `q` such that `q × that ≤ this`.
    Equivalently, it is the truncated quotient adjusted downward by 1 when the
    truncated remainder is non-zero and `this` and `that` have opposite signs.

    The invariant `this = fld(this, that) × that + mod(this, that)` holds.

    Special cases follow IEEE 754 conventions for division; in addition:
    - When the truncated remainder is NaN (e.g. ±∞ dividend), `fld` returns NaN.

    See also [`divrem`](#divrem), [`rem`](#rem), [`mod`](#mod).
    """
    _from(_MPFAlgo.fld(_ctx, _rep, that._rep), _ctx)


  fun mod(that: MPFloat): MPFloat =>
    """
    Floored remainder of `this / that`: the value `r` satisfying
    `this = fld(this, that) × that + r`.

    `r` has the same sign as `that` (the divisor).  The relationship
    `mod(this, that) = rem(this, that) + that` holds whenever the signs of
    `this` and `that` differ and the truncated remainder is non-zero.

    Special cases are inherited from `rem`; additionally, when `rem` returns
    NaN the result is NaN.

    See also [`divrem`](#divrem), [`rem`](#rem), [`fld`](#fld).
    """
    _from(_MPFAlgo.mod(_ctx, _rep, that._rep), _ctx)


  //- Unsafe arithmetic -------------------------------------------------------

  fun div_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe division of `this` by `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    _from(_MPFAlgo.div(_ctx, _rep, that._rep), _ctx)


  fun inv_unsafe(): MPFloat =>
    """
    Inverse of `this`. If `this` is +/- zero or NaN, the result is undefined.

    Calculate `1 / this` using Newton's method:  `x_{n+1} = x_n × (2 − G × x_n)`
    where `G = 256 × |this|` (the value of `this` reinterpreted with exponent 1).

    `G ∈ [1, 256)`, so the iteration converges to `y = 1/G ∈ (1/256, 1]`.
    The final result exponent is `(res._exponent + 1) − _exponent`.

    Sign propagation: `1/positive = positive`, `1/negative = negative`.
    A safety limit of `size × 4` iterations prevents non-convergence.

    TODO: Rounding management
    """
    _from(_MPFAlgo.inv(_ctx, _rep), _ctx)


  fun fld_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe floored division of `this` by `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    _from(_MPFAlgo.fld(_ctx, _rep, that._rep), _ctx)


  fun divrem_unsafe(that: MPFloat): (MPFloat, MPFloat) =>
    """
    Unsafe truncated division with remainder. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    divrem(that)


  fun rem_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe truncated remainder of `this / that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    rem(that)


  fun mod_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe floored remainder of `this / that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    mod(that)


  fun digit_shl(n: USize = 1): MPFloat =>
    """
    Drop the `n` most-significant base-256 digits and shift the remaining
    digits toward the most-significant position, padding with `n` zero bytes
    at the least-significant end.  The exponent is unchanged.

    When `n ≥ _size()` the result is all-zero (the value becomes exactly zero).
    Used internally by some transcendental algorithms; exposed for testing.
    """
    _from(_rep.digit_shl(n), _ctx)


  fun neg(): MPFloat =>
    """
    Negate `this`: flip the sign bit.
    """
    _from(MPFRep._create(
      not _rep.sign_bit(), _rep.is_nan(), _rep.is_infinite(),
      _rep.exponent(), _rep.raw_digits()), _ctx)


  fun neg_unsafe(): MPFloat =>
    """
    Negate `this` without special-value checking.
    """
    neg()


  fun sqrt_unsafe(): MPFloat =>
    """
    Unsafe square root. If `this` is negative, NaN, or ±∞, the result is undefined.
    """
    _from(_MPFAlgo.sqrt(_ctx, _rep), _ctx)


  //- Comparisons -------------------------------------------------------------

  fun eq(that: MPFloat): Bool =>
    """
    Equality of `this` and `that`.

    * NaN is never equal to anything, including itself
    * Infinities of the same sign are equal; `+∞ ≠ −∞`
    * `finite == ±∞` and `±∞ == finite` are `false`
    * `+0.0` and `-0.0` are equal
    """
    (not _rep.is_nan()) and
    (not that._rep.is_nan()) and
    (
      if _rep.is_infinite() or that._rep.is_infinite() then
        _rep.is_infinite() and that._rep.is_infinite() and
          (_rep.sign_bit() == that._rep.sign_bit())
      else
        eq_unsafe(that)
      end
    )


  fun ne(that: MPFloat): Bool =>
    """
    Difference of `this` and `that`, defined as `not (this == that)`.

    * `-0.0 != +0.0` is `false`
    * `NaN != x` is `true`. NaN is not-equal to everything, including itself.
    """
    not eq(that)


  fun lt(that: MPFloat): Bool =>
    """
    `this` is strictly inferior to `that`.

    * `-0.0 < +0.0` is `false`
    * `NaN < x` is `false`. `NaN` is unordered with everything.
    * `-inf < -inf` is false; `-inf < x` is true. −∞ is less than
    everything except itself and NaN.
    * `x < +inf` is true; `+inf < +inf` is false. +∞ is less than nothing.
    """
    if is_nan() or that.is_nan() then
      return false
    end
    if is_infinite() then
      if _rep.sign_bit() then
        return not (that.is_infinite() and that._rep.sign_bit())
      else
        return false
      end
    end
    if that.is_infinite() then
      if that._rep.sign_bit() then
        return false
      else
        return not(is_infinite() and not _rep.sign_bit())
      end
    end
    lt_unsafe(that)


  fun le(that: MPFloat): Bool =>
    """
    `this` is inferior or equal to `that`.

    * `NaN <= x` or `x <= NaN` are `false`
    * `+0.0 <= -0.0` is `true`
    """
    lt(that) or eq(that)


  fun ge(that: MPFloat): Bool =>
    """
    `this` is greater or equal to `that`.

    * `NaN >= x` or `x >= NaN` are `false`
    * `-0.0 >= +0.0` is `true`
    """
    if is_nan() or that.is_nan() then
      return false
    end
    not lt(that)


  fun gt(that: MPFloat): Bool =>
    """
    `this` is strictly superior to `that`.

    * `NaN > x` or `x > NaN` are `false`
    * `+0.0 > -0.0` is `false`
    """
    if is_nan() or that.is_nan() then
      return false
    end
    not le(that)


  fun compare(that: MPFloat): Compare =>
    """
    Compare `this` to `that`.

    See [`eq`](#eq) and [`lt`](#lt) to understand how `-inf`, `-0.0`, `+0.0`,
    `+inf` and `NaN` compare to each other.
    """
    if eq(that) then
      Equal
    elseif lt(that) then
      Less
    else
      Greater
    end


  //- Unsafe comparisons ------------------------------------------------------

  fun eq_unsafe(that: MPFloat): Bool =>
    """
    Equality of `this` and `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    if is_zero() and that.is_zero() then
      return true
    end
    if _rep.sign_bit() != that._rep.sign_bit() then
      return false
    end
    if _rep.exponent() != that._rep.exponent() then
      return false
    end
    // Compare the common prefix digit by digit.
    let na: USize = _size()
    let nb: USize = that._size()
    let ad: Array[U8] val = _rep.raw_digits()
    let bd: Array[U8] val = that._rep.raw_digits()
    try
      var i: USize = 0
      while i < na.min(nb) do
        if ad(i)? != bd(i)? then
          return false
        end
        i = i + 1
      end
      // If one array is longer, any extra nonzero byte breaks equality.
      if na > nb then
        var j: USize = nb
        while j < na do
          if ad(j)? != 0 then return false end
          j = j + 1
        end
      elseif nb > na then
        var j: USize = na
        while j < nb do
          if bd(j)? != 0 then return false end
          j = j + 1
        end
      end
    end
    true


  fun ne_unsafe(that: MPFloat): Bool =>
    """
    Difference of `this` and `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    not eq_unsafe(that)


  fun lt_unsafe(that: MPFloat): Bool =>
    """
    `this` is strictly inferior to `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    if is_zero() and that.is_zero() then
      return false
    end
    // Zero is less than any positive value, and greater than any negative value.
    if is_zero() then
      return not that._rep.sign_bit()  // 0 < that iff that > 0 (not negative)
    end
    if that.is_zero() then
      return _rep.sign_bit()  // this < 0 iff this is negative
    end
    let this_sign = _rep.sign_bit()
    let that_sign = that._rep.sign_bit()
    let na: USize = _size()
    let nb: USize = that._size()
    let ad: Array[U8] val = _rep.raw_digits()
    let bd: Array[U8] val = that._rep.raw_digits()
    if this_sign then
      if not that_sign then
        // negative < positive → always true
        return true
      end
      // Both negative: this < that ↔ |this| > |that|.
      if _rep.exponent() > that._rep.exponent() then
        return true
      end
      if _rep.exponent() < that._rep.exponent() then
        return false
      end
      try
        var i: USize = 0
        while i < na.min(nb) do
          let a_d: U8 = ad(i)?
          let b_d: U8 = bd(i)?
          if a_d > b_d then return true end
          if a_d < b_d then return false end
          i = i + 1
        end
        if na > nb then
          var j: USize = nb
          while j < na do
            if ad(j)? != 0 then return true end
            j = j + 1
          end
        end
      end
      false
    else
      if that_sign then
        // positive < negative → always false
        return false
      end
      // Both positive: compare exponents then digits.
      if _rep.exponent() < that._rep.exponent() then return true end
      if _rep.exponent() > that._rep.exponent() then return false end
      try
        var i: USize = 0
        while i < na.min(nb) do
          let a_d: U8 = ad(i)?
          let b_d: U8 = bd(i)?
          if a_d < b_d then return true end
          if a_d > b_d then return false end
          i = i + 1
        end
        if nb > na then
          var j: USize = na
          while j < nb do
            if bd(j)? != 0 then return true end
            j = j + 1
          end
        end
      end
      false
    end


  fun le_unsafe(that: MPFloat): Bool =>
    """
    `this` is inferior or equal to `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined as it is defined based on
    [`lt_unsafe`](#lt_unsafe) and [`eq_unsafe`](#eq_unsafe).
    """
    lt_unsafe(that) or eq_unsafe(that)


  fun ge_unsafe(that: MPFloat): Bool =>
    """
    `this` is greater or equal to `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined as it is defined based on
    [`lt_unsafe`](#lt_unsafe).
    """
    not lt_unsafe(that)


  fun gt_unsafe(that: MPFloat): Bool =>
    """
    `this` is strictly greater than `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined as it is defined based on
    [`lt_unsafe`](#lt_unsafe) and [`eq_unsafe`](#eq_unsafe).
    """
    not lt_unsafe(that) and not eq_unsafe(that)


  fun compare_unsafe(that: MPFloat): Compare =>
    """
    Compare `this` and `that` using unsafe comparisons. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined as it is defined based on
    [`lt_unsafe`](#lt_unsafe) and [`eq_unsafe`](#eq_unsafe).
    """
    if eq_unsafe(that) then
      Equal
    elseif lt_unsafe(that) then
      Less
    else
      Greater
    end


  //- Transcendental functions -------------------------------------------------

  fun ln(): MPFloat =>
    """
    Compute the natural logarithm of `this`.

    NaN → NaN. +∞ → +∞. −∞ → NaN. 0 → −∞. Negative → NaN.

    **Rounding**: correctly rounded — the result equals the infinite-precision
    `ln(this)` rounded under `rounding_mode()`. Ziv's iteration automatically
    extends working precision until the rounding direction is unambiguous.
    """
    _from(_MPFAlgo.ln(_ctx, _rep), _ctx)


  fun log(): MPFloat =>
    """
    Natural logarithm alias — same as `ln`. Correctly rounded.
    """
    ln()


  fun log2(): MPFloat =>
    """
    Compute `log₂(this)`.

    NaN → NaN. +∞ → +∞. 0 → −∞. Negative → NaN.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    _from(_MPFAlgo.log2(_ctx, _rep), _ctx)


  fun log10(): MPFloat =>
    """
    Compute `log₁₀(this)`.

    NaN → NaN. +∞ → +∞. 0 → −∞. Negative → NaN.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    _from(_MPFAlgo.log10(_ctx, _rep), _ctx)


  fun logb(base: MPFloat): MPFloat =>
    """
    Compute `log_base(this)`.

    NaN in either operand → NaN. +∞ → +∞. 0 → −∞. Negative → NaN.

    **Rounding**: faithfully rounded (≤ 1 ULP error). Computed as
    `ln(this)/ln(base)`; each `ln` is correctly rounded but the division
    introduces one additional rounding step.
    """
    _from(_MPFAlgo.logb(_ctx, _rep, base._rep), _ctx)


  fun exp(): MPFloat =>
    """
    Compute `e^this`.

    NaN → NaN. +∞ → +∞. −∞ → +0. 0 → 1.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    _from(_MPFAlgo.exp(_ctx, _rep), _ctx)


  fun exp2(): MPFloat =>
    """
    Compute `2^this`.

    NaN → NaN. +∞ → +∞. −∞ → +0. 0 → 1.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    _from(_MPFAlgo.exp2(_ctx, _rep), _ctx)


  fun powi(n: ILong, rnd: RoundingMode = RoundingNearest): MPFloat =>
    """
    Compute `this^n` for integer exponent `n`.

    NaN → NaN. ±∞ → depends on sign and parity of `n`.

    **Rounding**: correctly rounded via Ziv's iteration. The final result is
    the nearest representable value to `this^n` under `rounding_mode()`.
    """
    _from(_MPFAlgo.powi(_ctx, _rep, n), _ctx)


  fun pow(that: MPFloat): MPFloat =>
    """
    Compute `this^that`.

    For positive base: `exp(that × ln(this))`. For negative integer base:
    delegates to `powi`. NaN base or exponent → NaN.
    Negative base with non-integer exponent → NaN.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    _from(_MPFAlgo.pow(_ctx, _rep, that._rep), _ctx)


  fun sin(): MPFloat =>
    """
    Compute `sin(this)` (argument in radians).

    NaN or ±∞ → NaN. 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error) under `rounding_mode()`.
    """
    _from(_MPFAlgo.sin(_ctx, _rep), _ctx)


  fun cos(): MPFloat =>
    """
    Compute `cos(this)` (argument in radians).

    NaN or ±∞ → NaN. 0 → 1.

    **Rounding**: faithfully rounded (≤ 1 ULP error) under `rounding_mode()`.
    """
    _from(_MPFAlgo.cos(_ctx, _rep), _ctx)


  fun tan(): MPFloat =>
    """
    Compute `tan(this)` (argument in radians).

    NaN or ±∞ → NaN. 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error) under `rounding_mode()`.
    """
    _from(_MPFAlgo.tan(_ctx, _rep), _ctx)


  fun csc(): MPFloat =>
    """
    Compute `csc(this) = 1/sin(this)` (argument in radians).

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    sin().inv()


  fun sec(): MPFloat =>
    """
    Compute `sec(this) = 1/cos(this)` (argument in radians).

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    cos().inv()


  fun cot(): MPFloat =>
    """
    Compute `cot(this) = cos(this)/sin(this)` (argument in radians).

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    cos().div(sin())


  fun sinh(): MPFloat =>
    """
    Compute `sinh(this)`.

    NaN → NaN. ±∞ → ±∞. 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error) under `rounding_mode()`.
    """
    _from(_MPFAlgo.sinh(_ctx, _rep), _ctx)


  fun cosh(): MPFloat =>
    """
    Compute `cosh(this)`.

    NaN → NaN. ±∞ → +∞. 0 → 1.

    **Rounding**: faithfully rounded (≤ 1 ULP error) under `rounding_mode()`.
    """
    _from(_MPFAlgo.cosh(_ctx, _rep), _ctx)


  fun tanh(): MPFloat =>
    """
    Compute `tanh(this)`.

    NaN → NaN. +∞ → +1. −∞ → −1. 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error) under `rounding_mode()`.
    """
    _from(_MPFAlgo.tanh(_ctx, _rep), _ctx)


  fun csch(): MPFloat =>
    """
    Compute `csch(this) = 1/sinh(this)`.

    NaN → NaN. ±∞ → +0. 0 → NaN.

    **Rounding**: faithfully rounded (≤ 1 ULP error) under `rounding_mode()`.
    """
    _from(_MPFAlgo.csch(_ctx, _rep), _ctx)


  fun sech(): MPFloat =>
    """
    Compute `sech(this) = 1/cosh(this)`.

    NaN → NaN. ±∞ → +0.

    **Rounding**: faithfully rounded (≤ 1 ULP error) under `rounding_mode()`.
    """
    _from(_MPFAlgo.sech(_ctx, _rep), _ctx)


  fun coth(): MPFloat =>
    """
    Compute `coth(this) = cosh(this)/sinh(this)`.

    NaN → NaN. +∞ → +1. −∞ → −1. 0 → NaN.

    **Rounding**: faithfully rounded (≤ 1 ULP error) under `rounding_mode()`.
    """
    _from(_MPFAlgo.coth(_ctx, _rep), _ctx)


  fun abs(): MPFloat =>
    """
    Absolute value of `this`, when it makes sense. Result is undefined for NaN.
    """
    _from(MPFRep._create(false, _rep.is_nan(), _rep.is_infinite(),
      _rep.exponent(), _rep.raw_digits()), _ctx)


  fun trunc(): MPFloat =>
    """
    Truncation toward zero: the nearest integer to `this` in the direction of
    zero.  Equivalently, the integer part of `this` with the fractional
    base-256 bytes discarded.

    Representation: `value = 0.d₀d₁… × 256^e`.  The first `e` bytes
    (`d₀…d_{e−1}`) are the integer part; bytes at index `e` and beyond are
    fractional and are dropped.

    - NaN → NaN.  ±∞ → ±∞.  ±0 → ±0.
    - Purely fractional (`_exponent ≤ 0`) → +0.
    - All-integer (`_exponent ≥ _size()`) → `this` unchanged.
    - Mixed → keeps the first `_exponent` bytes, drops the rest.

    See also `floor`, `ceil`, `round`.
    """
    _from(_rep._trunc_frac(), _ctx)


  fun floor(): MPFloat =>
    """
    Floor: the largest integer less than or equal to `this`.

    - NaN → NaN.  ±∞ → ±∞.  ±0 → ±0.
    - Non-negative or exact integer → `trunc(this)`.
    - Negative with a non-zero fractional part → `trunc(this) − 1`.

    See also `trunc`, `ceil`, `round`.
    """
    let t = trunc()
    if _rep.sign_bit() and _rep._has_frac() then
      return (t - MPFloat.from[F64](1.0, get_precision(), _ctx.rounding))
    end
    t


  fun ceil(): MPFloat =>
    """
    Ceiling: the smallest integer greater than or equal to `this`.

    - NaN → NaN.  ±∞ → ±∞.  ±0 → ±0.
    - Non-positive or exact integer → `trunc(this)`.
    - Positive with a non-zero fractional part → `trunc(this) + 1`.

    See also `trunc`, `floor`, `round`.
    """
    let t = trunc()
    if (not _rep.sign_bit()) and _rep._has_frac() then
      return (t + MPFloat.from[F64](1.0, get_precision(), _ctx.rounding))
    end
    t


  fun round(): MPFloat =>
    """
    Round to the nearest integer, with ties broken away from zero (i.e.
    `round(0.5) = 1`, `round(−0.5) = −1`).

    - NaN → NaN.  ±∞ → ±∞.  ±0 → ±0.
    - Non-negative: `floor(this + 0.5)`.
    - Negative: `ceil(this − 0.5)`.

    The value 0.5 is exactly representable in base 256 (`[128]`, exponent 0),
    so no precision is lost by the shift before the floor/ceil call.

    See also `trunc`, `floor`, `ceil`.
    """
    let half = MPFloat.from[F64](0.5, get_precision(), _ctx.rounding)
    if _rep.sign_bit() then
      return (this  - half).ceil()
    end
    (this + half).floor()


  //- Hashing -----------------------------------------------------------------

  fun hash(): USize =>
    """
    Calculate a hash of this `MPFloat`.
    """
    _rep.hash()


  fun hash64(): U64 =>
    """
    Calculate a 64-bits hash of the `MPFloat`.
    """
    _rep.hash64()


  //- Limits ------------------------------------------------------------------

  fun \do_not_use\ tag precision2(): U8 =>
    """
    Mantissa precision in bits.

    DO NOT USE BECAUSE OF RISK OF OVERFLOW OF RESULT. RETURN 0
    """
    0

  fun \do_not_use\ tag precision10(): U8 =>
    """
    Mantissa precision in decimal digits.

    DO NOT USE BECAUSE OF RISK OF OVERFLOW OF RESULT. RETURN 0
    """
    0


  fun \do_not_use\ tag min_exp2(): I16 =>
    """
    Minimum exponent value such that (2^exponent) - 1 is representable at full
    precision (ie a normalised number).

    DO NOT USE BECAUSE OF RISK OF OVERFLOW OF RESULT. RETURN 0
    """
    0


  fun \do_not_use\ tag min_exp10(): I16 =>
    """
    Minimum exponent value such that (10^exponent) - 1 is representable at full
    precision (ie a normalised number).

    DO NOT USE BECAUSE OF RISK OF OVERFLOW OF RESULT. RETURN 0
    """
    0


  fun \do_not_use\ tag max_exp2(): I16 =>
    """
    Maximum exponent value such that (2^exponent) - 1 is representable.

    DO NOT USE BECAUSE OF RISK OF OVERFLOW OF RESULT. RETURN 0
    """
    0


  fun \do_not_use\ tag max_exp10(): I16 =>
    """
    Maximum exponent value such that (10^exponent) - 1 is representable.

    DO NOT USE BECAUSE OF RISK OF OVERFLOW OF RESULT. RETURN 0
    """
    0


  fun ldexp(x: MPFloat, e: I32): MPFloat =>
    """
    Multiply `MPFloat` x number by integral power of 2
    """
    MPFloat.create()


  fun frexp(): (MPFloat, U32) =>
    """
    Decomposes a `MPFloat` number into significand and base-2 exponent
    """
    (MPFloat.create(), 0)


//- Conversions ---------------------------------------------------------------

  fun f64(): F64 =>
    """
    Convert the current `MPFloat` to `F64`. Overflows are converted to ±∞.
    Underflow are converted to 0.

    This always hold:
    * If `f` is a `F64`: `MPFloat.from[F64](f).f64() == f`
    * If `mpf` is a `MPFloat`: `MPFloat.from[F64](mpf.f64()) == mpf`
    """
    _rep.f64()


  fun f32(): F32 =>
    """
    Convert the current `MPFloat` to `F32`. Overflows are converted to ±∞.
    Underflow are converted to 0.
    """
    _rep.f32()


  fun i128(): I128 =>
    """
    Convert `this` to an `I128`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `I128.max_value()` or `I128.min_value()`.
    """
    _rep.i128()


  fun i64(): I64 =>
    """
    Convert `this` to an `I64`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `I64.max_value()` or `I64.min_value()`.
    """
    _rep.i64()


  fun i32(): I32 =>
    """
    Convert `this` to an `I32`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `I32.max_value()` or `I32.min_value()`.
    """
    _rep.i32()


  fun i16(): I16 =>
    """
    Convert `this` to an `I16`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `I16.max_value()` or `I16.min_value()`.
    """
    _rep.i16()


  fun i8(): I8 =>
    """
    Convert `this` to an `I8`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `I8.max_value()` or `I8.min_value()`.
    """
    _rep.i8()


  fun ilong(): ILong =>
    """
    Convert `this` to an `ILong`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `ILong.max_value()` or `ILong.min_value()`.
    """
    _rep.ilong()


  fun isize(): ISize =>
    """
    Convert `this` to an `ISize`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `ISize.max_value()` or `ISize.min_value()`.
    """
    _rep.isize()


  fun i128_unsafe(): I128 =>
    """
    Convert `this` to an `I128` without overflow checking. For values outside
    the I128 range the result is undefined.
    """
    _rep.i128_unsafe()


  fun i64_unsafe(): I64 =>
    """
    Convert `this` to an `I64` without overflow checking. For values outside
    the I64 range the result is undefined.
    """
    _rep.i64_unsafe()


  fun i32_unsafe(): I32 =>
    """
    Convert `this` to an `I32` without overflow checking. For values outside
    the I32 range the result is undefined.
    """
    _rep.i32_unsafe()


  fun i16_unsafe(): I16 =>
    """
    Convert `this` to an `I16` without overflow checking. For values outside
    the I16 range the result is undefined.
    """
    _rep.i16_unsafe()


  fun i8_unsafe(): I8 =>
    """
    Convert `this` to an `I8` without overflow checking. For values outside
    the I8 range the result is undefined.
    """
    _rep.i8_unsafe()


  fun ilong_unsafe(): ILong =>
    """
    Convert `this` to an `ILong` without overflow checking. For values outside
    the ILong range the result is undefined.
    """
    _rep.ilong_unsafe()


  fun isize_unsafe(): ISize =>
    """
    Convert `this` to an `ISize` without overflow checking. For values outside
    the ISize range the result is undefined.
    """
    _rep.isize_unsafe()


  fun u128(): U128 =>
    """
    Convert `this` to a `U128`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `U128.max_value()` or 0.
    """
    _rep.u128()


  fun u64(): U64 =>
    """
    Convert `this` to a `U64`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `U64.max_value()` or 0.
    """
    _rep.u64()


  fun u32(): U32 =>
    """
    Convert `this` to a `U32`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `U32.max_value()` or 0.
    """
    _rep.u32()


  fun u16(): U16 =>
    """
    Convert `this` to a `U16`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `U16.max_value()` or 0.
    """
    _rep.u16()


  fun u8(): U8 =>
    """
    Convert `this` to a `U8`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `U8.max_value()` or 0.
    """
    _rep.u8()


  fun ulong(): ULong =>
    """
    Convert `this` to a `ULong`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `ULong.max_value()` or 0.
    """
    _rep.ulong()


  fun usize(): USize =>
    """
    Convert `this` to a `USize`. Fractional digits are truncated toward zero.
    In case of overflow the value is saturated to `USize.max_value()` or 0.
    """
    _rep.usize()
