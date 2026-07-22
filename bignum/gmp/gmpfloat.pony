// MPFloat using GMP/MPFR libraries
//
// The API of this implementation is limited and restricted to provinding
// a compatible interface with native MPFloat implementation. You need to
// rethink the API if you want to get full support for GMP/MPFR, particularly
// to benefit from these libraries memory optimisations.


use "debug"

use "../../assertx"
use "../../formatx"
use bignum = "../"
use "../../mathx"
use "../../pony_testx"

use "lib:gmp"
use "lib:mpfr"

// Size of MPFR types on amd64
// mpfr_int = 4
// mpfr_uint = 4
// mpfr_long = 8
// mpfr_ulong = 8
// mpfr_size_t = 8
// mpfr_prec_t = 8
// mpfr_sign_t = 4
// mpfr_exp_t = 8
// mpfr_uexp_t = 8
// mp_limb_t = 8
// mpfr_t = 32
// mpfr_t* = 8

// With Pony
// USize = 8
// ULong = 8
// ISize = 8
// ILong = 8


use @pony_ctx[Pointer[None]]()
use @pony_alloc[Pointer[U8]](ctx: Pointer[None], size: USize)


class val MPFloat is (bignum.BigReal[MPFloat] & Approximated[MPFloat, F64])
  """
  The implementation of `MPFloat` (Unlimited precision floating point numbers)
  using the [GNU Multiple Precision Arithmetic Library (GMP)](https://gmplib.org/)
  and the [GNU MPFR Library](https://www.mpfr.org/). Refer to
  [MPFR documentation](https://www.mpfr.org/mpfr-current/mpfr.html) to
  complement the explanations provided here.

  This implementation favours ease of use where the developer can translate
  mathematical formulas using literal Pony syntax. It provides automatic memory
  management.

  If a developer looks for raw performance and accuracy instead, like MPFR
  provides, it should turn to `MPF` class that provides a direct binding to
  MPFR, with the consequences that the developer will have to manage variables
  declaration and memory management. It is more difficult to translate
  mathematical formulas with `MPF` but the developer can fine tune roundings,
  reuse variables or use different precisions.

  The precision of a `MPFloat` is set when the instance is created. Contrarily
  to MPFR where the precision of an operation is defined by the result variable,
  it is expected that all `MPFloat`s involved in operations use the same
  precision. If that's not the case, warning messages are printed when compiled
  in debug mode, and results can be not the one expected by the developer. As a
  consequence, the precision of the result of the operation `this op that` or
  `this.op(that, ...)` is set to the precision of `this`.

  The default precision of 112 bits for the mantissa correspond to the size of
  the mantissa of a `F128`. That way, `MPFloat` can be used as a replacement for
  `F128` when the precision is not specified. In fact, there's a one bit
  difference in the mantissa size as `F128` assumes the most significant bit
  to be `1` and could be considered using 113 bits of mantissa.

  Usage of `MPFloat` class assumes that floats have the same precision and it
  prints warnings, when compiled in debug mode, on operations involving
  instances with different precisions. This is because accuracy can degrade
  when mixing operands with different precisions. But nothing prevent you from
  doing it, and there's even a constructor to create a new `MPFloat` from
  another one and changing its precision.
  """

  var _mpfr: SMPFr
    """
    The MPFT structure representing the float.
    """

  let _rnd: RoundingMode
    """
    The default rounding mode used when initializing the `MPFloat`.
    """

  //- Constructors ------------------------------------------------------------

  new val from_string(s: String = "",
                      prec: USize = 112,
                      base: U32 = 10,
                      rnd: RoundingMode = RoundingNearest) ? =>
    """
    Create a new `MPFloat` initialized to `s` value interpreted in base `base`,
    (default decimal base) with precision `prec` (default 112) and rounding mode
    `rnd` (nearest). If `s` is empty, the value of the `MPFloat` is `0`.

    When creating `MPFloat` from literal values, it's better to use that
    function than the `from_f64`, `from_f32`, `from_ulong`, `from_usize` functions
    where conversions occur and precision is lost.

    The base can be in the range `[2..36]`.

    The string format uses the following grammar:

    ```BNF
    mpfloat ::= sign? digit* decpoint? digit* exponent?
    sign ::= '+' | '-'
    digit ::= '0' | '1' | ... | '9' | 'a' | ... | 'z' | 'A' | ... | 'Z'
    decpoint := '.' | LOCALE_DEFINED_DECIMAL_SEPARATOR
    exponent ::= [ 'e' | 'E' | '@' ] sign? digit+
    ```

    The exponent part is always in decimal. For base larger than 10, the exponent
    character is `'@'`. For instance, you must write `1.234@56` instead of
    `1.234e56` if you are using base 16. For bases lower or equal to 10, one can
    use `e`, `E` or `@` to indicate the exponent part of the number.

    Values of `s` `"+inf"`, `"-inf"`, `"-0.0`" or `"+0.0"`, `"nan"` are
    correctly interpreted and give the corresponding `MPFloat`. These string
    values must be given exactly as documented. For instance, `"-0."` or `"-0.00"`
    won't be recognized as `+0.0`. These string literals are those used by Pony.
    MPFR/GMP can use different values (`@NaN@`, `-@Inf@` and `@Inf@`, `-0`) that
    are also accepted in this version of `MPFloat`.

    The whole string `s` must represent a floating point number. Leading or
    trailing white spaces are not accepted. In this implementation based on MPFR,
    digit separator `'_'` is not accepted.

    When the string can't be interpreted as a valid number in base `base`, an
    error is raised.

    If there is an overflow interpreting the value and it can't be coded with
    the desired precision, the number becomes `-inf` or `inf`.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    match s
    | "" => MPF.set_zero(_mpfr, 1)
    | "+inf" => MPF.set_inf(_mpfr, 1)
    | "@Inf@" => MPF.set_inf(_mpfr, 1)
    | "-inf" => MPF.set_inf(_mpfr, -1)
    | "-@Inf@" => MPF.set_inf(_mpfr, -1)
    | "0" => MPF.set_zero(_mpfr, 1)
    | "0.0" => MPF.set_zero(_mpfr, 1)
    | "+0.0" => MPF.set_zero(_mpfr, 1)
    | "-0" => MPF.set_zero(_mpfr, -1)
    | "-0.0" => MPF.set_zero(_mpfr, -1)
    | "nan" => MPF.set_nan(_mpfr)
    | "NaN" => MPF.set_nan(_mpfr)
    | "@NaN@" => MPF.set_nan(_mpfr)
    else
      let r = MPF.set_str(_mpfr, s, base.i32(), rnd)
      if r < 0 then
        error
      end
      r
    end


  new val from[N: (Number & Real[N])](n: N,
                                      prec: USize = 112,
                                      rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` with value `n` and precision `prec` (default 112),
    using rounding mode `rnd` (default nearest).

    When accuracy is important, you should use `create` constructor.

    If there is an overflow interpreting the value and it can't be coded with
    the desired precision, the number becomes `-inf` or `inf`.

    Internally, this implementation uses a conversion to a `String` and then
    conversion to `MPFloat`.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    let tmp = try
        from_string(n.string(), prec, 10, rnd)?
      else
        Debug("MPFloat.from: Conversion failed. Converting from F64 with possible " +
              "lost accuracy")
        from_f64(n.f64(), prec, rnd)
      end
    MPF.set(_mpfr, tmp._mpfr, rnd)


  new val from_f64(d: F64 = 0.0,
                   prec: USize = 112,
                   rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` initialized to `d` value (default 0.0) with precision
    `prec` (default 112) and rounding mode `rnd` (default nearest).

    Values of `d` +inf, -inf, -0.0 or +0.0, NaN are correctly treated and give
    the corresponding `MPFloat`.

    If there is an overflow interpreting the value and it can't be coded with
    the desired precision, the number becomes `-inf` or `inf`.

    When accuracy is important and with precision larger than 53 bits, which is
    the case of the default precision, you should use `create` constructor.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    if d.infinite() then
      if d < 0.0 then
        MPF.set_inf(_mpfr, -1)
      else
        MPF.set_inf(_mpfr, 1)
      end
    elseif (d == 0.0) and (d < 0.0) then
      MPF.set_zero(_mpfr, -1)
    elseif (d == 0.0) and (d > 0.0) then
      MPF.set_zero(_mpfr, 1)
    elseif d.nan() then
      MPF.set_nan(_mpfr)
    else
      MPF.set_d(_mpfr, d, rnd)
    end


  new val from_f32(f: F32 = 0.0,
                   prec: USize = 112,
                   rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` initialized to `f` value (default 0.0) with precision
    `prec` (default 112) and rounding mode `rnd` (default nearest).

    Values of `f` +inf, -inf, -0.0 or +0.0, NaN are correctly treated and give
    the corresponding `MPFloat`.

    If there is an overflow interpreting the value and it can't be coded with
    the desired precision, the number becomes `-inf` or `inf`.

    This function is specific to MPFR implementation.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    if f.infinite() then
      if f < 0.0 then
        MPF.set_inf(_mpfr, -1)
      else
        MPF.set_inf(_mpfr, 1)
      end
    elseif (f == 0.0) and (f < 0.0) then
      MPF.set_zero(_mpfr, -1)
    elseif (f == 0.0) and (f > 0.0) then
      MPF.set_zero(_mpfr, 1)
    elseif f.nan() then
      MPF.set_nan(_mpfr)
    else
      MPF.set_flt(_mpfr, f, rnd)
    end


  new val from_ulong(u: ULong = 0,
                     prec: USize = 112,
                     rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` initialized to `u` value (default 0) with precision
    `prec` (default 112) and rounding mode `rnd` (default nearest).

    If there is an overflow interpreting the value and it can't be coded with
    the desired precision, the number becomes `-inf` or `inf`.

    This function is specific to MPFR implementation.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    MPF.set_ui(_mpfr, u, rnd)


  new val from_usize(u: USize = 0,
                     prec: USize = 112,
                     rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` initialized to `u` value (default 0) with precision
    `prec` (default 112) and rounding mode `rnd` (default nearest).

    If there is an overflow interpreting the value and it can't be coded with
    the desired precision, the number becomes `-inf` or `inf`.

    This function is specific to MPFR implementation.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    MPF.set_uj(_mpfr, u, rnd)


  new val from_ilong(i: ILong = 0,
                     prec: USize = 112,
                     rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` initialized to `i` value (default 0) with precision
    `prec` (default 112) and rounding mode `rnd` (default nearest).

    If there is an overflow interpreting the value and it can't be coded with
    the desired precision, the number becomes `-inf` or `inf`.

    This function is specific to MPFR implementation.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    MPF.set_si(_mpfr, i, rnd)


  new val from_isize(i: ISize = 0,
                     prec: USize = 112,
                     rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` initialized to `i` value (default 0) with precision
    `prec` (default 112) and rounding mode `rnd` (default nearest).

    If there is an overflow interpreting the value and it can't be coded with
    the desired precision, the number becomes `-inf` or `inf`.

    This function is specific to MPFR implementation.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    MPF.set_sj(_mpfr, i, rnd)


  new val from_mpfloat(f: MPFloat,
                       prec: USize = 112,
                       rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` whose value is equal to `f` but with precision `prec`
    (default 112) and using rounding mode `rnd` (default nearest). This
    constructor is useful when you want to change the precision of the
    resulting `MPFloat`.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    MPF.set(_mpfr, f._mpfr, _rnd)


  new val nan_val(prec: USize = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a `MPFloat` representing `NaN` (Not a Number).
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    MPF.set_nan(_mpfr)


  new val inf_val(positive: Bool, prec: USize = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a `MPFloat` representing `+inf` (when `positive` is `true`) or
    `-inf` (when `positive` is `false`).
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    MPF.set_inf(_mpfr, if positive then 1 else -1 end)


  //- Finalizer ---------------------------------------------------------------
  // We need to use a finalizer as the memory in `_mpfr` is managed by the MPFR
  // library.

  fun _final() =>
    """
    Free the memory allocated by MPFR.
    """
    MPF.clear(_mpfr)


  //- Conversions -------------------------------------------------------------

  fun f64(): F64 =>
    """
    Convert the `MPFloat` to a `F64` using the stored rounding mode.
    If overflow occurs during the conversion, the result is a special value like
    `-inf`, `+inf`, `-0.0`, `+0.0` or `nan`.
    """
    MPF.get_d(_mpfr, _rnd)


  fun f32(): F32 =>
    """
    Convert the `MPFloat` to a `F32` using the stored rounding mode.
    If overflow occurs during the conversion, the result is a special value like
    `-inf`, `+inf`, `-0.0`, `+0.0` or `nan`.
    """
    MPF.get_flt(_mpfr, _rnd)


  fun isize(): ISize ? =>
    """
    Convert the `MPFloat` to a `ISize` using the stored rounding mode.
    If the conversion can't occur because the float does not fit into
    a `ISize`, then an error is raised.
    """
    if MPF.fits_intmax_p(_mpfr, _rnd) then
      MPF.get_sj(_mpfr, _rnd)
    else
      error
    end


  fun usize(): USize ? =>
    """
    Convert the `MPFloat` to a `USize` using the stored rounding mode.
    If the conversion can't occur because the float does not fit into
    a `USize`, then an error is raised.
    """
    if MPF.fits_uintmax_p(_mpfr, _rnd) then
      MPF.get_uj(_mpfr, _rnd)
    else
      error
    end


  fun ilong(): ILong ? =>
    """
    Convert the `MPFloat` to a `ILong` using the stored rounding mode.
    If the conversion can't occur because the float does not fit into
    a `ILong`, then an error is raised.
    """
    if MPF.fits_slong_p(_mpfr, _rnd) then
      MPF.get_si(_mpfr, _rnd)
    else
      error
    end


  fun ulong(): ULong ? =>
    """
    Convert the `MPFloat` to a `ULong` using the stored rounding mode.
    If the conversion can't occur because the float does not fit into
    a `ULong`, then an error is raised.
    """
    if MPF.fits_ulong_p(_mpfr, _rnd) then
      MPF.get_ui(_mpfr, _rnd)
    else
      error
    end


  fun string(): String iso^ =>
    """
    Convert the `MPFloat` to a string using decimal representation and with
    all digits of the associated precision, in fixed point or scientific
    format. Converting the value and reading it back does not always yield the
    original value. This string conversion tries to find a representation
    that is short and easy to use for mathematicians. Use `exact_string`
    if you need to control rounding and accuracy.
    """
    let p = ((get_precision().f64() * F64(2.0).log()) / F64(10.0).log()).ceil().u64()
    let fmt = ("%." + p.string() + "RNg").cstring()
    let buf_size = @mpfr_snprintf(Pointer[U8], 0, fmt, _mpfr) + 1
    if buf_size < 0 then
      // We return an error message instead of the string representation
      return "[MPFloat.string] MPFR error".clone()
    end
    let result: String iso = recover iso
      // Memory managed by Pony
      let buffer: Pointer[U8] = @pony_alloc(@pony_ctx(), buf_size.usize())
      @mpfr_snprintf(buffer, buf_size.u64(), fmt, _mpfr)
      String.copy_cstring(buffer)
    end
    consume result


  fun exact_string(base: U8 = 10, rnd: RoundingMode = RoundingNearest)
                  : (String iso^, I64, Bool) =>
    """
    Convert the `MPFloat` to a `String` using base `base` (default 10) and
    rounding mode `rnd` (default to nearest). The result is a tuple where

    * `_1` is the string representation of the mantissa (significand).
    * `_2` is the exponent
    * `_3` is `true` if the conversion is not exact.

    The base can be in the range `[2..36]`.

    The generated string is a fraction, with an implicit radix point immediately
    to the left of the first digit. For example, the number `-3.1416` would
    be returned as `"-31416"` in the String `_1` and `1` in the exponent `_2`.

    If the current float is a special value, like `+inf`, `-inf` or `nan`, the
    value of the exponent is not defined.
    """
    ifdef debug then
      ((2 <= base) and (base <= 36)) or
        Fail(Format("[MPFloat.string] The base ({}) must be in the range [2..36]", base))
    end
    (let mantissa, let exponent) = MPF.get_str(base.i32(), get_precision().ilong(), _mpfr, rnd)

    // Format special values to Pony's format
    mantissa.replace("@Inf@", "inf")
    mantissa.replace("@NaN@", "nan")

    let inexact = MPF.inexflag_p()
    (consume mantissa, exponent, inexact)


  //- Constants ---------------------------------------------------------------

  new val pi(prec: USize = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Calculate Pi ($$\pi = 3.14159...$$) with the requested precision `prec`
    (default 112) and rounding mode `rnd` (default to nearest).

    MPFR caches the value of Pi.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    MPF.const_pi(_mpfr, rnd)


  new val e(prec: USize = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Calculate the base of Neperian logarithm e ($$\e = 2.718...$$) with the
    requested precision `prec` (default 112) and rounding mode `rnd` (default
    to nearest).
    """
    _mpfr = SMPFr
    _rnd = rnd
    let one = from_usize(1, prec, rnd)
    MPF.init2(_mpfr, prec.ilong())
    MPF.exp(_mpfr, one._mpfr, rnd)


  new val const_log2(prec: USize = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Calculate the value of the logarithm of 2 (`log(2)`) with the requested
    precision `prec` (default 112) and rounding mode `rnd` (default nearest).

    MPFR caches the result.
    """
    _mpfr = SMPFr
    _rnd = rnd
    MPF.init2(_mpfr, prec.ilong())
    MPF.const_log2(_mpfr, rnd)


  //- Miscellaneous functions -------------------------------------------------

  fun mpfr_version(): String iso^ =>
    """
    Return the version of MPFR used by `MPFloat` as a string.
    """
    MPF.get_version()


  fun mpfr_patches(): String iso^ =>
    """
    Return the patch level of the MPFR library used by `MPFloat` as a string.
    """
    MPF.get_patches()


  fun get_precision(): USize =>
    """
    Get the precision of the `MPFloat` mantissa (significand) in bits.
    """
    MPF.get_prec(_mpfr).usize()


  fun get_rounding_mode(): RoundingMode =>
    """
    Get the rounding mode that was used when initializing the `MPFloat`.
    """
    _rnd


  fun sign(): Compare =>
    """
    Return the sign of `this`. If `this < 0`, returns `Less`; if `this > 0`,
    returns `Greater`; and returns `Equal` when `this` is `NaN`.
    """
    let sgn = MPF.sgn(_mpfr)
    if sgn < 0 then
      Less
    elseif sgn > 0 then
      Greater
    else
      Equal
    end


  fun get_base(): USize =>
    """
    Get the value of the base used internally for encoding the `MPFloat`. In the
    case of MPFR, the floats are encoded in base 2.
    """
    2


  fun get_exponent(): I64 ? =>
    """
    Get the value of the exponent in base `get_base` used internally for
    encoding the `MPFloat`. This is defined only when `this` is a number.
    An error is raised if that's not the case.
    """
    if number() then
      MPF.get_exp(_mpfr)
    else
      error
    end


  fun get_exponent10(): I64 ? =>
    """
    Get an approximation of the value of the exponent of the current `MPFloat`
    in base 10. This is defined only if `this` is a number. An error is raised
    if that's not the case.

    This is an approximation calculated from `get_exponent` and `get_base` for
    a quick result with evaluating the mantissa. You must use `exact_string`
    in case you need the exact result.
    """
    ((get_exponent()?.f64() / F64(10.0).log()) * get_base().f64()).i64()


  fun clone(): MPFloat =>
    """
    Create a copy of `this`.

    From the view point of MPFR, it assigns the value of `this` to a newly
    created `MPFloat`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.set(result._mpfr, _mpfr, _rnd)
    result


  fun next_toward(that: MPFloat): MPFloat =>
    """
    Create the next `MPFloat` value from `this` toward `that`. If `this` or
    `that` is `NaN`, the result is `NaN`. If `this` and `that` are equal,
    the result is unchanged. Otherwise, if `this` is different from `that`,
    the result is the next floating-point number (with the precision of `this`
    and the current exponent range) in the direction of `that` (the infinite
    values are seen as the smallest and largest floating-point numbers). If
    the result is zero, it keeps the same sign.
    """
    let result = clone()
    MPF.nexttoward(result._mpfr, that._mpfr)
    result


  fun next_below(): MPFloat =>
    """
    Return the next floating point number with the same precision as `this`
    inferior to `this`.
    """
    let result = clone()
    MPF.nextbelow(result._mpfr)
    result


  fun next_above(): MPFloat =>
    """
    Return the next floating point number with the same precision as `this`
    superior to `this`.
    """
    let result = clone()
    MPF.nextabove(result._mpfr)
    result


  fun min(that: MPFloat): MPFloat =>
    """
    Return the minimal value between `this` and `that`. If both floats are
    `NaN`, then the result is `NaN`. If one of the floats is `NaN`, then
    return the numeric value. If `this` and `that` are zeros of different
    signs, then the result is `-0`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.min(result._mpfr, _mpfr, that._mpfr, _rnd)
    result


  fun max(that: MPFloat): MPFloat =>
    """
    Return the maximal value between `this` and `that`. If both floats are
    `NaN`, then the result is `NaN`. If one of the floats is `NaN`, then
    return the numeric value. If `this` and `that` are zeros of different
    signs, then the result is `+0`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.max(result._mpfr, _mpfr, that._mpfr, _rnd)
    result


  //- Comparisons -------------------------------------------------------------

  fun eq(that: MPFloat): Bool =>
    """
    Is `this == that`? Both floats must have the same precision.
    """
    ifdef debug then
      let this_prec = get_precision()
      let that_prec = that.get_precision()
      if this_prec != that_prec then
        Debug("[MPFloat.eq] The two floats have different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). Comparison can be incorrect")
      end
    end
    MPF.equal_p(_mpfr, that._mpfr)


  fun ne(that: MPFloat): Bool =>
    """
    Is `this != that`? Both floats must have the same precision.
    """
    ifdef debug then
      let this_prec = get_precision()
      let that_prec = that.get_precision()
      if this_prec != that_prec then
        Debug("[MPFloat.ne] The two floats have different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). Comparison can be incorrect")
      end
    end
    not MPF.equal_p(_mpfr, that._mpfr)


  fun lt(that: MPFloat): Bool =>
    """
    Is `this < that`? Both floats must have the same precision.
    """
    ifdef debug then
      let this_prec = get_precision()
      let that_prec = that.get_precision()
      if this_prec != that_prec then
        Debug("[MPFloat.lt] The two floats have different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). Comparison can be incorrect")
      end
    end
    MPF.less_p(_mpfr, that._mpfr)


  fun le(that: MPFloat): Bool =>
    """
    Is `this <= that`? Both floats must have the same precision.
    """
    ifdef debug then
      let this_prec = get_precision()
      let that_prec = that.get_precision()
      if this_prec != that_prec then
        Debug("[MPFloat.le] The two floats have different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). Comparison can be incorrect")
      end
    end
    MPF.lessequal_p(_mpfr, that._mpfr)


  fun ge(that: MPFloat): Bool =>
    """
    Is `this >= that`? Both floats must have the same precision.
    """
    ifdef debug then
      let this_prec = get_precision()
      let that_prec = that.get_precision()
      if this_prec != that_prec then
        Debug("[MPFloat.ge] The two floats have different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). Comparison can be incorrect")
      end
    end
    MPF.greaterequal_p(_mpfr, that._mpfr)


  fun gt(that: MPFloat): Bool =>
    """
    Is `this > that`? Both floats must have the same precision.
    """
    ifdef debug then
      let this_prec = get_precision()
      let that_prec = that.get_precision()
      if this_prec != that_prec then
        Debug("[MPFloat.gt] The two floats have different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). Comparison can be incorrect")
      end
    end
    MPF.greater_p(_mpfr, that._mpfr)


  fun finite(): Bool =>
    """
    Is `this` a number, that is not infinite nor NaN?
    """
    MPF.number_p(_mpfr)


  fun number(): Bool =>
    """
    Is `this` a regular number, that is not an infinity nor NaN nor `-0.0` or
    `0.0`.
    """
    MPF.regular_p(_mpfr)


  fun infinite(): Bool =>
    """
    Is `this` infinite, that is equal to `-inf` or `+inf`?
    """
    MPF.inf_p(_mpfr)


  fun nan(): Bool =>
    """
    Is `this` not a number? We can't write `this == nan` as a comparison to NaN
    will always fail. Some cases of NaN:
    * `+inf - +inf` or `+inf + -inf`
    * `0 / 0` or `inf / inf`
    * `log(x)` with `x < 0`

    See [Wikipedia page](https://en.wikipedia.org/wiki/NaN) for more information.
    """
    MPF.nan_p(_mpfr)


  fun zero(): Bool =>
    """
    Is `this == 0.0`?
    """
    MPF.zero_p(_mpfr)


  fun is_nan(): Bool =>
    """Canonical alias for `nan()`."""
    MPF.nan_p(_mpfr)

  fun is_infinite(): Bool =>
    """Canonical alias for `infinite()`."""
    MPF.inf_p(_mpfr)

  fun is_finite(): Bool =>
    """Canonical alias for `finite()`."""
    MPF.number_p(_mpfr)

  fun is_zero(): Bool =>
    """Canonical alias for `zero()`."""
    MPF.zero_p(_mpfr)

  fun is_negative(): Bool =>
    """Is `this < 0`? Returns false for NaN and zero."""
    MPF.sgn(_mpfr) < 0

  fun is_integer(): Bool =>
    """Is `this` a finite integer value?"""
    if not is_finite() then return false end
    MPF.integer_p(_mpfr)

  fun compare(that: MPFloat): Compare =>
    """
    Three-way comparison. Returns `Less`, `Equal`, or `Greater`.
    Consistent with `eq` and `lt`: NaN comparisons follow the same rules.
    """
    if eq(that) then Equal
    elseif lt(that) then Less
    else Greater
    end

  fun almost_eq(that: box->MPFloat,
                rel_tol: F64 = F64.epsilon().sqrt(),
                abs_tol: F64 = F64.epsilon().sqrt()): Bool =>
    """
    Returns true if `|this - that| <= max(rel_tol * max(|this|, |that|), abs_tol)`.
    NaN arguments always return false. Inf returns true only when both are the
    same infinity.
    """
    if is_nan() or that.is_nan() then return false end
    if is_infinite() or that.is_infinite() then
      return is_infinite() and that.is_infinite() and (f64() == that.f64())
    end
    let diff = sub(that).abs().f64()
    let m = (abs().f64().max(that.abs().f64())) * rel_tol
    diff <= m.max(abs_tol)

  fun ln(): MPFloat =>
    """Alias for `log()` — natural logarithm."""
    log()


  //- Arithmetic --------------------------------------------------------------

  fun add(that: MPFloat): MPFloat =>
    """
    Addition `this + that` using the stored rounding mode.
    The result is a new `MPFloat` with a precision equal to the precision of
    `this`.

    If `this` or `that` is infinity, the result is infinity. If `this` or `that`
    is `NaN`, the result is `NaN`.

    If the precision of `that` is less than the one of `this`, a warning
    is printed when compiled in `debug` mode.
    """
    let this_prec = get_precision()
    ifdef debug then
      let that_prec = that.get_precision()
      if that_prec < this_prec then
        Debug("[MPFloat.add] Trying to add two floats with different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). The result with precision " + this_prec.string() + " is inexact")
      end
    end
    let result = MPFloat.from_f64(0.0, this_prec, _rnd)
    MPF.add(result._mpfr, _mpfr, that._mpfr, _rnd)
    result


  fun sub(that: MPFloat): MPFloat =>
    """
    Subtraction `this - that` using the stored rounding mode.
    The result is a new `MPFloat` with a precision equal to the precision of
    `this`.

    If `this` or `that` is infinity, the result is infinity. If `this` or `that`
    is `NaN`, the result is `NaN`.

    If the precision of `that` is less than the one of `this`, a warning
    is printed when compiled in `debug` mode.
    """
    let this_prec = get_precision()
    ifdef debug then
      let that_prec = that.get_precision()
      if that_prec < this_prec then
        Debug("[MPFloat.sub] Trying to subtract two floats with different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). The result with precision " + this_prec.string() + " is inexact")
      end
    end
    let result = MPFloat.from_f64(0.0, this_prec, _rnd)
    MPF.sub(result._mpfr, _mpfr, that._mpfr, _rnd)
    result


  fun mul(that: MPFloat): MPFloat =>
    """
    Multiplication `this * that` using the stored rounding mode.
    The result is a new `MPFloat` with a precision equal to the precision of
    `this`.

    If `this` or `that` is infinity, the result is infinity. If `this` or `that`
    is `NaN`, the result is `NaN`.

    If the precision of `that` is less than the one of `this`, a warning
    is printed when compiled in `debug` mode.
    """
    let this_prec = get_precision()
    ifdef debug then
      let that_prec = that.get_precision()
      if that_prec < this_prec then
        Debug("[MPFloat.mul] Trying to multiply two floats with different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). The result with precision " + this_prec.string() + " can be less precise")
      end
    end
    let result = MPFloat.from_f64(0.0, this_prec, _rnd)
    MPF.mul(result._mpfr, _mpfr, that._mpfr, _rnd)
    result


  fun div(that: MPFloat): MPFloat =>
    """
    Division `this / that` using the stored rounding mode.
    The result is a new `MPFloat` with a precision equal to the precision of
    `this`.

    If `this` and `that` are infinity, the result is `NaN`. If `this` or `that`
    is `NaN`, the result is `NaN`. If `that` is infinity, the result is `0`.
    If `this` and `that` are `0`, the result is `NaN`.

    If the precision of `that` is less than the one of `this`, a warning
    is printed when compiled in `debug` mode.

    Note: This differs from Pony float `div` that is rounded toward zero
    (`RoundingZero`).
    """
    let this_prec = get_precision()
    ifdef debug then
      let that_prec = that.get_precision()
      if that_prec < this_prec then
        Debug("[MPFloat.div] Trying to divide two floats with different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). The result with precision " + this_prec.string() + " can be less precise")
      end
    end
    let result = MPFloat.from_f64(0.0, this_prec, _rnd)
    MPF.div(result._mpfr, _mpfr, that._mpfr, _rnd)
    result


  fun sqrt(): MPFloat =>
    """
    Square root of `this` using the stored rounding mode.
    The result has the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.sqrt(result._mpfr, _mpfr, _rnd)
    result


  fun cbrt(): MPFloat =>
    """
    Cubic root of `this` using the stored rounding mode.
    The result has the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.cbrt(result._mpfr, _mpfr, _rnd)
    result


  fun rootn(n: USize): MPFloat =>
    """
    Calculate the n-th root of `this` using the stored rounding mode.
    The result has the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.rootn_ui(result._mpfr, _mpfr, n.ulong(), _rnd)
    result


  fun neg(): MPFloat =>
    """
    Return `-this`, the negation of the float value. The result has the same
    precision and the same rounding mode.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.neg(result._mpfr, _mpfr, _rnd)
    result


  fun abs(): MPFloat =>
    """
    Return the absolute value of the float. The result has the same precision
    and the same rounding mode.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.abs(result._mpfr, _mpfr, _rnd)
    result


  fun inv(): MPFloat =>
    """
    Return the multiplicative inverse `1/this`. The result has the same
    precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.ui_div(result._mpfr, 1, _mpfr, _rnd)
    result


  fun copysign(sgn: MPFloat): MPFloat =>
    """
    Return a `MPFloat` with the value of `this` but with the sign of `sgn`.
    The result has the same precision and the same rounding mode. The result is
    a new `MPFloat`, even when there is no reason to change the sign.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.copysign(result._mpfr, _mpfr, sgn._mpfr, _rnd)
    result


  //- Transcendental functions ------------------------------------------------

  fun log(): MPFloat =>
    """
    Return natural logarithm of `this` using the stored rounding mode.
    The result has the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.log(result._mpfr, _mpfr, _rnd)
    result


  fun log2(): MPFloat =>
    """
    Return base-2 logarithm of `this` using the stored rounding mode.
    The result has the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.log2(result._mpfr, _mpfr, _rnd)
    result


  fun log10(): MPFloat =>
    """
    Return decimal logarithm of `this` using the stored rounding mode.
    The result has the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.log10(result._mpfr, _mpfr, _rnd)
    result


  fun exp(): MPFloat =>
    """
    Return the exponential of `this` using the stored rounding mode.
    The result has the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.exp(result._mpfr, _mpfr, _rnd)
    result


  fun exp2(): MPFloat =>
    """
    Return 2 to the power of `this` using the stored rounding mode.
    The result has the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.exp2(result._mpfr, _mpfr, _rnd)
    result


  fun pow(that: MPFloat): MPFloat =>
    """
    Return `this` raised to the power of `that` using the stored rounding mode.
    The result is a new `MPFloat` with a precision equal to the precision of
    `this`.

    If the precision of `that` is less than the one of `this`, a warning
    is printed when compiled in `debug` mode.
    """
    let this_prec = get_precision()
    ifdef debug then
      let that_prec = that.get_precision()
      if that_prec != this_prec then
        Debug("[MPFloat.pow] The two MPFloats have different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). The result with precision " + this_prec.string() + " can be less precise")
      end
    end
    let result = MPFloat.from_f64(0.0, this_prec, _rnd)
    MPF.pow(result._mpfr, _mpfr, that._mpfr, _rnd)
    result


  fun powi(n: ILong): MPFloat =>
    """
    Return `this` raised to the integer power `n` using the stored rounding
    mode. The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.pow_si(result._mpfr, _mpfr, n, _rnd)
    result


  fun cos(): MPFloat =>
    """
    Return the cosine of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.cos(result._mpfr, _mpfr, _rnd)
    result


  fun sin(): MPFloat =>
    """
    Return the sine of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.sin(result._mpfr, _mpfr, _rnd)
    result


  fun tan(): MPFloat =>
    """
    Return the tangent of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.tan(result._mpfr, _mpfr, _rnd)
    result


  fun acos(): MPFloat =>
    """
    Return the arc-cosine of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.acos(result._mpfr, _mpfr, _rnd)
    result


  fun asin(): MPFloat =>
    """
    Return the arc-sine of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.asin(result._mpfr, _mpfr, _rnd)
    result


  fun atan(): MPFloat =>
    """
    Return the arc-tangent of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.atan(result._mpfr, _mpfr, _rnd)
    result


  fun cosh(): MPFloat =>
    """
    Return the hyperbolic cosine of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.cosh(result._mpfr, _mpfr, _rnd)
    result


  fun sinh(): MPFloat =>
    """
    Return the hyperbolic sine of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.sinh(result._mpfr, _mpfr, _rnd)
    result


  fun tanh(): MPFloat =>
    """
    Return the hyperbolic tangent of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.tanh(result._mpfr, _mpfr, _rnd)
    result


  fun acosh(): MPFloat =>
    """
    Return the hyperbolic arc-cosine of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.acosh(result._mpfr, _mpfr, _rnd)
    result


  fun asinh(): MPFloat =>
    """
    Return the hyperbolic arc-sine of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.asinh(result._mpfr, _mpfr, _rnd)
    result


  fun atanh(): MPFloat =>
    """
    Return the hyperbolic arc-tangent of `this` using the stored rounding mode.
    The result is a new `MPFloat` with the same precision as `this`.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.atanh(result._mpfr, _mpfr, _rnd)
    result


  //- Format string -----------------------------------------------------------

  fun format(fmt: String): String iso^ =>
    """
    Format the float according to `fmt` specification.
    """
    // TODO
    "TODO".clone()


  //- Integer and remainder related functions ---------------------------------

  fun ceil(): MPFloat =>
    """
    Return a new `MPFloat` equal to `this` rounded to the next higher or equal
    representable integer, with the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.ceil(result._mpfr, _mpfr)
    result


  fun floor(): MPFloat =>
    """
    Return a new `MPFloat` equal to `this` rounded to the next lower or equal
    representable integer, with the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.floor(result._mpfr, _mpfr)
    result


  fun round(): MPFloat =>
    """
    Return a new `MPFloat` equal to `this` rounded to the nearest representable
    integer, rounding halfway cases away from zero, with the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.round(result._mpfr, _mpfr)
    result


  fun trunc(): MPFloat =>
    """
    Return a new `MPFloat` equal to `this` rounded to the next representable
    integer toward zero, with the same precision.
    """
    let result = MPFloat.from_f64(0.0, get_precision(), _rnd)
    MPF.trunc(result._mpfr, _mpfr)
    result


  fun mod(that: MPFloat): MPFloat =>
    """
    Return the remainder of the integral division of `this` by `that`, with the
    same precision, using the stored rounding mode.

    If the precision of `that` is less than the one of `this`, a warning
    is printed when compiled in `debug` mode.
    """
    let this_prec = get_precision()
    ifdef debug then
      let that_prec = that.get_precision()
      if that_prec != this_prec then
        Debug("[MPFloat.mod] The two MPFloats have different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). The results with precision " + this_prec.string() + " can be less precise")
      end
    end
    let result = MPFloat.from_f64(0.0, this_prec, _rnd)
    MPF.fmod(result._mpfr, _mpfr, that._mpfr, _rnd)
    result


  fun divrem(that: MPFloat): (MPFloat, MPFloat) =>
    """
    Return the quotient and the remainder of the integral division of `this`
    by `that` in a tuple, with the same precision, using the stored rounding
    mode.

    If the precision of `that` is less than the one of `this`, a warning
    is printed when compiled in `debug` mode.
    """
    let this_prec = get_precision()
    ifdef debug then
      let that_prec = that.get_precision()
      if that_prec != this_prec then
        Debug("[MPFloat.divrem] The two MPFloats have different precisions (" +
              this_prec.string() + " != " + that_prec.string() +
              "). The results with precision " + this_prec.string() + " can be less precise")
      end
    end
    let remainder = mod(that)
    let quotient = MPFloat.from_f64(0.0, this_prec, _rnd)
    let frac = MPFloat.from_f64(0.0, this_prec, _rnd)
    MPF.modf(quotient._mpfr, frac._mpfr, div(that)._mpfr, _rnd)
    (quotient, remainder)


  //- Precision change --------------------------------------------------------

  fun val change_precision(prec: USize, rnd: RoundingMode = RoundingNearest): MPFloat =>
    """
    Change the precision of the current `MPFloat` to `prec` and round it using
    the `rnd` rounding mode (default is nearest). If `prec` is greater than or
    equal to the current precision, then new space is allocated for the mantissa
    (significand), and it is filled with zeros. Otherwise, the mantissa
    (significand) is rounded to precision `prec` with the given direction `rnd`.
    """
    MPF.prec_round(_mpfr, prec.ilong(), rnd)
    this


  fun min_precision(): ULong =>
    """
    Return the minimum number of bits required to store the mantissa (significand),
    and `0` for special values including `0`.
    """
    MPF.min_prec(_mpfr).ulong()
