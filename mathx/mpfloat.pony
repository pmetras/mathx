// Multi-precision floating point numbers

use "../assertx"

use "collections"
use "debug"


class val MPFloat
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

  The default precision of 112 bits for the mantissa correspond to
  the size of the mantissa of a `F128`. That way, `MPFloat` can be used as a
  replacement for `F128` when the precision is not specified.

  Usage of `MPFloat` class assumes that floats have the same precision and it
  prints warnings, when compiled in debug mode, on operations involving
  instances with different precisions. This is because accuracy can degrade
  when mixing operands with different precisions. But nothing prevent you from
  doing it, and there's even a constructor to create a new `MPFloat` from
  another one and changing its precision.
  """

  let _sign: Bool
    """
    Sign of the number: `true` when the value is strictly negative.
    Ignored when `_nan` is `true`. For −0 both `_sign` is `true` and the
    digit array is all-zero.
    """

  let _nan: Bool
    """
    `true` when the number is Not-a-Number. When set, `_inf`, `_sign` and
    `_digits` are irrelevant.
    """

  let _inf: Bool
    """
    `true` when the number is infinite. Sign is given by `_sign`.
    Ignored when `_nan` is `true`.
    """

  let _exponent: I64
    """
    Base-256 exponent. The value of a finite non-special `MPFloat` is:

      0.d₀d₁… × 256^_exponent

    For the common arithmetic operations, `_exponent = 1` is used throughout
    (meaning the first digit `d₀` is the integer part). Proper exponent
    normalisation for values outside [0, 256) will be added in a later step.
    """

  let _digits: Array[U8] val
    """
    The arbitrary-precision fractional digits encoded in base 256, big-endian
    (index 0 is most significant). For a non-zero finite number, `_digits(0)`
    is guaranteed to be non-zero (in normalised representation). For zero, the
    array may be all-zeros or empty.
    """

  let _base: U32 = 256
    """
    The base of the digits.
    """

  let _base_bits: USize = 8
    """
    The number of bits in a digit
    """


  let _rounding: RoundingMode
    """
    The rounding mode of the `MPFloat`.
    """

  new val _create(
    sgn: Bool,
    nan: Bool,
    inf: Bool,
    expn: I64,
    digits: Array[U8] val,
    rnd: RoundingMode /* TODO TO FIND BUGS = RoundingNearest */)
  =>
    """
    Private canonical constructor. All public constructors and operations
    produce their result through this constructor.
    """
    _sign = sgn
    _nan = nan
    _inf = inf
    _exponent = expn
    _digits = digits
    _rounding = rnd


  new val create(prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a positive zero with `prec` bits of mantissa (significand).

    This is the default constructor: `MPFloat()` gives a size-112 bits positive zero,
    and `MPFloat(n)` gives a size-n bits positive zero (all digits zero).
    """
    _sign = false
    _nan = false
    _inf = false
    _exponent = 0
    let p = (prec.usize() + (_base_bits - 1)) / _base_bits
    _digits = Array[U8].init(0, p)
    _rounding = rnd


  new val nan_val() =>
    """
    Create a Not-a-Number value.
    """
    _sign = false
    _nan = true
    _inf = false
    _exponent = 0
    _digits = Array[U8].create()
    _rounding = RoundingNearest


  new val inf_val(positive: Bool = true) =>
    """
    Create an infinite value. Pass `positive = false` for −∞.
    """
    _sign = not positive
    _nan = false
    _inf = true
    _exponent = 0
    _digits = Array[U8].create()
    _rounding = RoundingNearest


  new val from_f64(f: F64, prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` from the `F64` value `f` with `prec` bits of
    precision (default 112, giving ~33 decimal digits) and a default rounding
    mode to nearest digit.

    Special values (NaN, ±∞, ±0) are preserved. The conversion normalises `f`
    so that `0.d₀d₁… × 256^_exponent` with `d₀ ≠ 0` for non-zero values.

    The rounding mode is not currently used by the `MPFloat`.
    """
    if f.nan() then
      _sign = false
      _nan = true
      _inf = false
      _exponent = 0
      _digits = Array[U8].create()
      _rounding = RoundingNearest
    elseif f.infinite() then
      _sign = f < 0.0
      _nan = false
      _inf = true
      _exponent = 0
      _digits = Array[U8].create()
      _rounding = RoundingNearest
    elseif f == 0.0 then
      // Distinguish −0 from +0 via the sign bit.
      _sign = (f.bits() == 0x8000000000000000)
      _nan = false
      _inf = false
      _exponent = 0
      let p = (prec.usize() + (_base_bits - 1)) / _base_bits
      _digits = Array[U8].init(0, p)
      _rounding = rnd
    else
      _sign = f < 0.0
      _nan = false
      _inf = false

      // Normalise exactly using frexp
      (let m, let e_u32) = f.abs().frexp()
      let e = e_u32.i32().i64()
      
      // value = m * 2^e = 0.d0d1... * 256^expn = frac * 2^(8*expn)
      // frac = m * 2^(e - 8*expn)
      // We want 1/256 <= frac < 1. 
      // Since 0.5 <= m < 1, e - 8*expn must be in [-7, 0].
      let expn = (e.f64() / 8.0).ceil().i64()
      let shift = e - (expn * 8)
      var frac = m.f64() * F64(2).pow(shift.f64())
      _exponent = expn

      // Extract prec bits (mapped to base-256 digits)
      let p: USize = ((prec.usize() + (_base_bits - 1)) / _base_bits).max(1)
      let base = _base.f64()
      _digits = recover
        let d = Array[U8].init(0, p)
        var i: USize = 0
        while i < p do
          frac = frac * base
          let di: U8 = frac.u8()
          try d.update(i, di)? end
          frac = frac - di.f64()
          i = i + 1
        end
        d
      end
      _rounding = rnd
    end


  new val from_f32(f: F32, prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` from the `F32` value `f` with `prec` bits of
    precision (default 112, giving ~33 decimal digits) and a default rounding
    mode to nearest digit.

    Special values (NaN, ±∞, ±0) are preserved. The conversion normalises `f`
    so that `0.d₀d₁… × 256^_exponent` with `d₀ ≠ 0` for non-zero values.

    The rounding mode is not currently used by the `MPFloat`.
    """
    if f.nan() then
      _sign = false
      _nan = true
      _inf = false
      _exponent = 0
      _digits = Array[U8].create()
      _rounding = RoundingNearest
    elseif f.infinite() then
      _sign = f < 0.0
      _nan = false
      _inf = true
      _exponent = 0
      _digits = Array[U8].create()
      _rounding = RoundingNearest
    elseif f == 0.0 then
      // Distinguish −0 from +0 via the sign bit.
      _sign = (f.bits() == 0x80000000)
      _nan = false
      _inf = false
      _exponent = 0
      let p = (prec.usize() + (_base_bits - 1)) / _base_bits
      _digits = Array[U8].init(0, p)
      _rounding = rnd
    else
      _sign = f < 0.0
      _nan = false
      _inf = false

      // Normalise exactly using frexp
      (let m, let e_u32) = f.abs().frexp()
      let e = e_u32.i32().i64()
      
      // value = m * 2^e = 0.d0d1... * 256^expn = frac * 2^(8*expn)
      // frac = m * 2^(e - 8*expn)
      // We want 1/256 <= frac < 1. 
      // Since 0.5 <= m < 1, e - 8*expn must be in [-7, 0].
      let expn = (e.f32() / 8.0).ceil().i64()
      let shift = e - (expn * 8)
      var frac = m * F32(2).pow(shift.f32())
      _exponent = expn

      // Extract prec bits (mapped to base-256 digits)
      let p: USize = ((prec.usize() + (_base_bits - 1)) / _base_bits).max(1)
      let base = _base.f32()
      _digits = recover
        let d = Array[U8].init(0, p)
        var i: USize = 0
        while i < p do
          frac = frac * base
          let di: U8 = frac.u8()
          try d.update(i, di)? end
          frac = frac - di.f32()
          i = i + 1
        end
        d
      end
      _rounding = rnd
    end


  new val from_string(s: String = "",
                      prec: ULong = 112,
                      base: U8 = 10,
                      rnd: RoundingMode = RoundingNearest) ? =>
    """
    Create a new `MPFloat` by parsing the string `s` with `prec` bits
    of precision (default 112, ≈33 significant decimal digits). Raises
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
    // Spaces are not significant.
    let st: String = s.clone() .> strip()

    if base != 10 then
      error  // TODO: implement parsing for bases other than 10
    end

    let p_digits = (prec.usize() + (_base_bits - 1)) / _base_bits

    // Empty string → +0
    if st.size() == 0 then
      _sign = false
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, p_digits)
      _rounding = rnd
      return
    end

    // Special values.
    if (st == "nan") or (st == "NaN") or (st == "@NaN@") then
      _sign = false
      _nan = true
      _inf = false
      _exponent = 0
      _digits = Array[U8].create()
      _rounding = RoundingNearest
      return
    end

    if (st == "+inf") or (st == "@Inf@") or (st == "inf") then
      _sign = false
      _nan = false
      _inf = true
      _exponent = 0
      _digits = Array[U8].create()
      _rounding = RoundingNearest
      return
    end

    if (st == "-inf") or (st == "-@Inf@") then
      _sign = true
      _nan = false
      _inf = true
      _exponent = 0
      _digits = Array[U8].create()
      _rounding = RoundingNearest
      return
    end

    if (st == "0") or (st == "0.0") or (st == "+0") or (st == "+0.0") then
      _sign = false
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, p_digits)
      _rounding = rnd
      return
    end

    if (st == "-0") or (st == "-0.0") then
      _sign = true
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, p_digits)
      _rounding = rnd
      return
    end

    // Multi-precision decimal parsing using exact MPInt arithmetic.
    //
    // Phase 1 — Digit accumulation: collect up to `sig_limit` significant
    // decimal digits into an MPInt N via Horner's method (N = N × 10 + digit).
    // MPInt arithmetic is exact, so no rounding error is introduced here.
    //
    // Phase 2 — Scaling: value = N × 10^dec_exp.
    //   dec_exp ≥ 0 → multiply N by 10^dec_exp (exact MPInt), then from_mpint.
    //   dec_exp < 0 → value = N / (2^n × 5^n), n = −dec_exp.
    //     Dividing by 2^n is a free bit-shift; dividing by 5^n is done via
    //     MPInt long division at (prec + ⌈n / 8⌉ + 4) bytes of working precision.
    //     This is exact when N is divisible by 5^n (e.g. "2.0", "1.5", "0.25")
    //     and correctly rounded to nearest otherwise (e.g. "0.1").
    //
    // log10(256^p_digits) ≈ p_digits × 2.408 decimal digits; add 4 guard digits.
    let sig_limit: USize = ((p_digits.f64() * 2.41).usize() + 4).max(1)

    var pos: USize = 0
    let sz: USize = st.size()
    var fsign: Bool = false

    if st(pos)? == '-' then
      fsign = true
      pos = pos + 1
    elseif st(pos)? == '+' then
      pos = pos + 1
    end

    // Accumulators.
    let ten_int: MPInt = MPInt.from[ILong](10)
    var n_int: MPInt = MPInt.from[ILong](0)
    var sig_count: USize = 0
    // int_exp: extra powers of 10 from integer digits skipped beyond sig_limit.
    var int_exp: I64 = 0
    // frac_count: decimal places consumed (significant or leading-zero).
    var frac_count: I64 = 0
    var has_digit: Bool = false

    // Integer digits.
    var sep_seen: Bool = false
    while (pos < sz) and (((st(pos)? >= '0') and (st(pos)? <= '9')) or (st(pos)? == '_')) do
      // Digits separator must be unique
      if st(pos)? == '_' then
        if not sep_seen then
          sep_seen = true
          pos = pos + 1
          continue
        else
          Assert(false, "[MPFloat.from_string] You can't have multiple consecutive separators '_' to group digits in mantissa.")?
        end
      end
      sep_seen = false

      // Digits
      let d: U8 = st(pos)? - '0'
      has_digit = true
      if (d != 0) or (sig_count > 0) then
        if sig_count < sig_limit then
          // TODO DELETE n_int = n_int.mul(ten_int).add(MPInt.from[ILong](d.ilong()))
          n_int = (n_int * ten_int) + MPInt.from[ILong](d.ilong())
          sig_count = sig_count + 1
        else
          int_exp = int_exp + 1
        end
      end
      pos = pos + 1
    end

    // Fractional digits.
    if (pos < sz) and (st(pos)? == '.') then
      pos = pos + 1

      sep_seen = false
      while (pos < sz) and (((st(pos)? >= '0') and (st(pos)? <= '9')) or (st(pos)? == '_')) do
        // Digits separator must be unique
        if st(pos)? == '_' then
          if not sep_seen then
            sep_seen = true
            pos = pos + 1
            continue
          else
            Assert(false, "[MPFloat.from_string] You can't have multiple consecutive separators '_' to group digits in decimals.")?
          end
        end
        sep_seen = false

        let d: U8 = st(pos)? - '0'
        has_digit = true
        if (d != 0) or (sig_count > 0) then
          if sig_count < sig_limit then
            // TODO DELETE n_int = n_int.mul(ten_int).add(MPInt.from[ILong](d.ilong()))
            n_int = (n_int * ten_int) + MPInt.from[ILong](d.ilong())
            sig_count = sig_count + 1
            frac_count = frac_count + 1
          end
          // Digits beyond sig_limit in the fractional part fall below the
          // requested precision and are dropped.
        else
          // Leading fractional zeros shift the effective decimal point.
          frac_count = frac_count + 1
        end
        pos = pos + 1
      end
    end

    if not has_digit then
      error
    end

    // Decimal exponent ('e', 'E', or '@' for base ≠ 10).
    var str_exp: I64 = 0
    if (pos < sz) and
       ((st(pos)? == 'e') or (st(pos)? == 'E') or (st(pos)? == '@'))
    then
      pos = pos + 1
      var esign: Bool = false
      if (pos < sz) and (st(pos)? == '-') then
        esign = true
        pos = pos + 1
      elseif (pos < sz) and (st(pos)? == '+') then
        pos = pos + 1
      end
      var ehas_digit: Bool = false
      sep_seen = false
      while (pos < sz) and (((st(pos)? >= '0') and (st(pos)? <= '9')) or (st(pos)? == '_')) do
        // Digits separator must be unique
        if st(pos)? == '_' then
          if not sep_seen then
            sep_seen = true
            pos = pos + 1
            continue
          else
            Assert(false, "[MPFloat.from_string] You can't have multiple consecutive separators '_' to group digits in exponent.")?
          end
        end
        sep_seen = false

        str_exp = (str_exp * 10) + (st(pos)? - '0').i64()
        ehas_digit = true
        pos = pos + 1
      end

      if not ehas_digit then
        error
      end

      if esign then
        str_exp = -str_exp
      end
    end

    // Reject trailing garbage.
    if pos != sz then
      error
    end

    // True value = n_int × 10^dec_exp.
    let dec_exp: I64 = (int_exp - frac_count) + str_exp

    // ── Zero check ─────────────────────────────────────────────────────────
    if n_int.is_zero() then
      _sign = fsign
      _nan  = false
      _inf  = false
      _exponent = 0
      _digits = Array[U8].init(0, p_digits)
      _rounding = rnd
      return
    end

    // ── Scale n_int by 10^dec_exp ───────────────────────────────────────────
    // Positive dec_exp: convert n_int exactly to MPFloat, then multiply by
    // 10^dec_exp via truncated binary exponentiation (O(log dec_exp) steps,
    // each bounded to p2 bytes). This is accurate regardless of dec_exp size.
    //
    // Negative dec_exp (n = −dec_exp):
    // For small n (≤ n_max_exact): exact MPInt division by 5^n at
    // sufficient working precision — gives exact results when n_int is
    // divisible by 5^n (e.g. "2.0" → exactly 2, "1.5" → exactly 1.5).
    // For large n: start from exact n_int MPFloat, divide by 10^n via
    // truncated binary exponentiation (MPFloat path has Newton error only in
    // the last few ULPs of the mantissa, which is acceptable for large n).
    //
    // Working precision p2_digits = p_digits + 2 suppresses rounding error in
    // intermediate multiplications.
    let p2_digits: USize = p_digits + 2
    let p2_bits: ULong = (p2_digits * _base_bits).ulong()
    // Threshold for exact MPInt division: 5^n must be a manageable size.
    // For n ≤ n_max_exact, five_pow and n_shifted fit within ~200 bytes.
    let n_max_exact: USize = ((p2_digits.f64() * 14.0).usize() + 50).min(300)

    // Helper: binary-exponentiation scale of ten_mp^sn, truncated to p2_digits.
    let ten_mp = MPFloat.from_f64(10.0, p2_bits, rnd)

    let tmp: MPFloat =
      if dec_exp >= 0 then
        // Start from exact MPInt → MPFloat conversion (no Horner error), then
        // scale up by 10^dec_exp with per-step truncation to p2_digits.
        var n_mp = MPFloat.from_mpint(n_int, p2_bits, rnd)

        if dec_exp > 0 then
          var scale = MPFloat.from_f64(1.0, p2_bits, rnd)
          var sbase: MPFloat = ten_mp
          var sn: I64 = dec_exp
          while sn > 0 do
            if (sn and 1) == 1 then
              scale = (scale * sbase)._trunc(p2_digits)
            end
            sbase = (sbase * sbase)._trunc(p2_digits)
            sn = sn / 2
          end
          n_mp = (n_mp * scale)._trunc(p2_digits)
        end
        n_mp._trunc(p_digits)
      else
        // dec_exp < 0: value = n_int / (2^n × 5^n), n = −dec_exp.
        let n: USize = (-dec_exp).usize()
        if n <= n_max_exact then
          // Exact path: compute q = ⌊ n_int × 2^{_base_bits × extra − n} / 5^n ⌋
          // = ⌊ value × 256^extra ⌋.
          //
          // extra is chosen so that q has at least p_digits significant bytes:
          //   extra ≥ n / log₁₀(256) + p_digits + guard
          //        ≈ n × 5/12 + p_digits + 4.
          // shift_bits = extra × _base_bits − n is always positive within the threshold.
          let extra: USize = (((n * 5) + 11) / 12) + p_digits + 4
          let shift_bits: USize = (extra * _base_bits) - n
          let five_pow: MPInt = MPInt.from[ILong](5).pow(MPInt.from[ILong](n.ilong()))
          let n_shifted: MPInt = n_int.shl(MPInt.from[ILong](shift_bits.ilong()))
          (let q, let r) = n_shifted.divrem(five_pow)
          // Round to nearest: increment q if 2×r ≥ five_pow.
          let q_rounded: MPInt =
            // TODO DELETE if r.add(r).ge(five_pow) then
            if (r + r) >= five_pow then
              // TODO DELETE q.add(MPInt.from[ILong](1))
              q + MPInt.from[ILong](1)
            else
              q
            end

          // q_rounded encodes value in units of 1/256^extra.
          // Big-endian bytes → MPFloat: e = byte_count − extra.
          let qbytes: Array[U8] val = q_rounded.raw_digits()
          let qlen: USize = qbytes.size()
          let new_exp: I64 = qlen.i64() - extra.i64()
          let keep: USize = p_digits.min(qlen)

          let new_digits: Array[U8] val = recover
            let d = Array[U8].create(keep)
            var i: USize = 0
            while i < keep do
              try d.push(qbytes(i)?) end
              i = i + 1
            end
            d
          end
          MPFloat._create(false, false, false, new_exp, new_digits, rnd)
        else
          // Large n: exact MPInt → MPFloat start, divide by 10^n via
          // truncated binary exponentiation.
          var n_mp = MPFloat.from_mpint(n_int, p2_bits, rnd)
          var scale = MPFloat.from_f64(1.0, p2_bits, rnd)
          var sbase: MPFloat = ten_mp
          var sn: I64 = (-dec_exp)
          while sn > 0 do
            if (sn and 1) == 1 then
              scale = (scale * sbase)._trunc(p2_digits)
            end
            sbase = (sbase * sbase)._trunc(p2_digits)
            sn = sn / 2
          end
          (n_mp / scale)._trunc(p_digits)
        end
      end

    _sign = fsign
    _nan  = tmp._nan
    _inf  = tmp._inf
    _exponent = tmp._exponent
    _digits = tmp._digits
    _rounding = rnd


  new val from_mpint(n: MPInt, prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` from the `MPInt` value `n` with `prec` bits of
    precision (default 112, giving ~33 decimal digits) and rounding mode `rnd`.

    The conversion is exact up to the requested precision: the magnitude of
    `n` is represented without error as long as `n` fits within `prec` bits
    (i.e. `|n| < 2^prec`); larger values are truncated to the `prec` bits.
    As a consquence, if you want to keep all the digits of `n` in the resulting
    `MPFloat`, select a `prec` value that is at least 3.32 × number of decimal
    digits of `n`.

    Special cases:
    - Zero → `+0` (positive zero regardless of any sign on the MPInt zero).
    - Sign is preserved: a negative `MPInt` produces a negative `MPFloat`.

    Algorithm: call `MPInt.raw_digits()` to get the absolute value as
    a big-endian `Array[U8]` (each base-65536 word split into two bytes, MSW
    first, leading zeros stripped).  The result maps directly onto the MPFloat
    `_digits` layout:
    - `_exponent` = total byte count of the full magnitude (before truncation).
    - `_digits` = the first `prec` bits (truncated to byte boundary).

    This is O(n) in the word count.

    The rounding mode `rnd` is accepted for API compatibility; rounding support
    is a future TODO.
    """
    _rounding = rnd

    let p_digits = (prec.usize() + (_base_bits - 1)) / _base_bits

    if n.is_zero() then
      _sign = false
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, p_digits)
      return
    end

    _sign = n.is_negative()
    _nan  = false
    _inf  = false

    // Get the magnitude as big-endian bytes (leading zeros already stripped).
    let mag = n.raw_digits()
    let total: USize = mag.size()

    // _exponent = number of bytes in the full magnitude, so that
    //   value = 0.d₀d₁… × 256^_exponent  with d₀ = mag(0) ≠ 0.
    _exponent = total.i64()

    // Keep only the `p_digits` most-significant bytes.
    _digits = recover
      let keep: USize = total.min(p_digits)
      let d = Array[U8].create(keep)
      mag.copy_to(d, 0, 0, keep)
      d
    end


  new val from_mpfloat(f: MPFloat,
                       prec: ULong = 112,
                       rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` whose value is equal to `f` but with precision `prec`
    bits (default 112) and using rounding mode `rnd` (default nearest). This
    constructor is useful when you want to change the precision and the rounding
    mode of the resulting `MPFloat`.

    TODO: Manage rounding. Presently, the result is truncated to fit with
    requested precision.
    """
    _sign = f._sign
    _nan = f._nan
    _inf = f._inf
    _exponent = f._exponent
    _rounding = rnd

    let p_digits = (prec.usize() + (_base_bits - 1)) / _base_bits

    _digits = recover
      let size = p_digits.min(f._size())
      let d = Array[U8].init(0, p_digits)
      f._digits.copy_to(d, 0, 0, size)
      d
    end


  new val from_ulong(n: ULong, prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` whose value is equal to `n` with precision `prec` bits (default 112)
    and rounding mode `rnd` (default nearest).

    TODO: Implement rounding
    """
    _sign = false
    _nan = false
    _inf = false
    _rounding = rnd

    let p_digits = (prec.usize() + (_base_bits - 1)) / _base_bits

    let base: ULong = _base.ulong()
    var q: ULong = n

    _digits = recover
      let d: Array[U8] = Array[U8].create(p_digits)
      var i: USize = 0
      while q >= base do
        (q, let r) = q.divrem(base)
        d.push(r.u8())
        i = i + 1
      end
      d.push(q.u8())
      d.reverse_in_place()
      d
    end
    // _exponent = number of bytes in the full magnitude, so that
    // n = 0.d₀d₁… × 256^_exponent  with d₀ = mag(0) ≠ 0.
    _exponent = _digits.size().i64()



  new min_normalized(prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    The smallest normalized floating point number.

    As don't have the notion of normalized number in the sense of the binary
    floating point representation, the smallest normalized `MPFloat` has been
    set to `0.0`.
    """
    _sign = false
    _nan = false
    _inf = false
    _exponent = 0
    let p = (prec.usize() + (_base_bits - 1)) / _base_bits
    _digits = Array[U8].init(0, p)
    _rounding = rnd


  new epsilon(prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Create the *machine epsilon* for the given precision: the smallest positive
    `MPFloat` ε such that `1 + ε ≠ 1` in base 256, i.e. `ε = 256^(1 − p_digits)`.

    For the default precision of 112 bits (14 bytes, ≈ 33 significant decimal digits):
    `ε ≈ 4.93 × 10^{−32}`.

    This is the base-256 analogue of `F64.epsilon() ≈ 2.22 × 10^{−16}`.
    In general `MPFloat.epsilon(p)` decreases as `p` grows — the user controls
    precision and thus controls the machine epsilon.

    Example use in tests: compare with absolute tolerance `2 × epsilon`:
    ```
    let eps = MPFloat.epsilon(p)
    h.assert_true(got.almost_eq(expected, MPFloat.from_f64(1e-31), MPFloat.from_f64(1e-31)))
    ```
    """
    _sign = false
    _nan = false
    _inf = false
    let p_digits: USize = ((prec.usize() + (_base_bits - 1)) / _base_bits).max(1)
    // ε = 0.1_{256} × 256^{2−p_digits} = (1/256) × 256^{2−p_digits} = 256^{1−p_digits}
    _exponent = 2 - p_digits.i64()
    _digits = recover
      let a = Array[U8].init(0, p_digits)
      try a(0)? = 1 end
      a
    end
    _rounding = rnd


  new val min_value() =>
    """
    The minimum value is −∞ (`-inf`).
    """
    _sign = true   // negative: −∞
    _nan = false
    _inf = true
    _exponent = 0
    _digits = Array[U8].create()
    _rounding = RoundingNearest



  new val max_value() =>
    """
    The maximum value is +∞ (`+inf`).
    """
    _sign = false  // positive: +∞
    _nan = false
    _inf = true
    _exponent = 0
    _digits = Array[U8].create()
    _rounding = RoundingNearest


  new val pi(prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    The pi constant, calculated with the specified `prec` accuracy (number of
    bits, default 112 giving ~33 decimal digits). The rounding
    mode `rnd` is not used yet (TODO).

    Uses Machin's formula:
      π = 16·arctan(1/5) − 4·arctan(1/239)
    where each arctan is computed by the Taylor series
      arctan(x) = x − x³/3 + x⁵/5 − x⁷/7 + ⋯
    Both arguments are small (< 0.25), giving fast geometric convergence.
    """
    let p_digits: USize = (prec.usize() + (_base_bits - 1)) / _base_bits
    let p: USize = p_digits + 4   // guard digits for intermediate results
    let p_bits: ULong = (p * _base_bits).ulong()

    // Compute 1/5 and 1/239 to full p-digit precision via MPFloat division.
    let one = MPFloat.from_f64(1.0, p_bits, rnd)
    let x5 = MPFloat.from_f64(5.0, p_bits, rnd).inv()._trunc(p) //one.div(MPFloat.from_f64(5.0,   p, rnd))._trunc(p)
    let x239 = MPFloat.from_f64(239.0, p_bits, rnd).inv()._trunc(p) // one.div(MPFloat.from_f64(239.0, p, rnd))._trunc(p)

    // arctan(1/5) via Taylor series: arctan(x) = x - x³/3 + x⁵/5 - ...
    // Keep `pow5` as the raw odd power x^(2k+1); divide separately for each term.
    let neg_x52: MPFloat = (-x5 * x5)._trunc(p) //x5.mul(x5).neg()._trunc(p)  // −(1/5)²
    var pow5: MPFloat = x5
    var atan5: MPFloat = x5
    var k5: USize = 1
    var iters5: USize = 0
    let max_iters: USize = p * 30 // TODO: Why 30?
    while iters5 < max_iters do
      pow5 = pow5.mul(neg_x52)._trunc(p)   // next odd power: −x^(2k+1)
      let d5: USize = (2 * k5) + 1
      let term5: MPFloat = if d5 <= 255 then
          (let q5, _) = pow5._short_div(d5.u8())
          q5
        else
          pow5.div(MPFloat.from_f64(d5.f64(), p_bits, rnd))._trunc(p)
        end
      let new_a5: MPFloat = atan5.add(term5)._trunc(p)
      if new_a5.eq(atan5) then break end
      atan5 = new_a5
      k5 = k5 + 1
      iters5 = iters5 + 1
    end

    // arctan(1/239) via Taylor series.
    let neg_x2392: MPFloat = x239.mul(x239).neg()._trunc(p) // −(1/239)²
    var pow239: MPFloat = x239
    var atan239: MPFloat = x239
    var k239: USize = 1
    var iters239: USize = 0
    while iters239 < max_iters do
      pow239 = pow239.mul(neg_x2392)._trunc(p)
      let d239: USize = (2 * k239) + 1
      let term239: MPFloat = if d239 <= 255 then
          (let q239, _) = pow239._short_div(d239.u8())
          q239
        else
          pow239.div(MPFloat.from_f64(d239.f64(), p_bits, rnd))._trunc(p)
        end
      let new_a239: MPFloat = atan239.add(term239)._trunc(p)
      if new_a239.eq(atan239) then break end
      atan239 = new_a239
      k239 = k239 + 1
      iters239 = iters239 + 1
    end

    // π = 16·arctan(1/5) − 4·arctan(1/239)
    let sixteen = MPFloat.from_f64(16.0, p_bits, rnd)
    let four = MPFloat.from_f64(4.0, p_bits, rnd)
    let pi_val: MPFloat = sixteen.mul(atan5).sub(four.mul(atan239))._trunc(p_digits)

    _sign     = false
    _nan      = false
    _inf      = false
    _exponent = pi_val._exponent
    _digits   = pi_val._digits
    _rounding = rnd


new val pi_chudnovsky(prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
  """
  Chudnovsky’s Algorithm: https://en.wikipedia.org/wiki/Chudnovsky_algorithm

  This formula converge very quickly but because of numerous factorial calculations
  in numerator and denominator, loss of precision occurs also very quickly.
  ⚠ Don't use it as is!
  """
    let p_digits: USize = (prec.usize() + (_base_bits - 1)) / _base_bits
    let p: USize = p_digits + 4   // guard digits for intermediate results
    let p_bits: ULong = (p * _base_bits).ulong()
    let one = MPFloat.from_ulong(1, p_bits, rnd)
    let minus_one = -one

    let k_3 = MPFloat.from_ulong(3, p_bits, rnd)
    let k_1_5 = k_3 / MPFloat.from_ulong(2, p_bits, rnd)
    let k_640320 = MPFloat.from_ulong(640320, p_bits, rnd)
    let k_13591409 = MPFloat.from_ulong(13591409, p_bits, rnd)
    let k_545140134 = MPFloat.from_ulong(545140134, p_bits, rnd)

    var result = MPFloat(p_bits, rnd)
    var prev_res = result
    var k: ULong = 0
    var term = MPFloat(p_bits, rnd)
    repeat
      (let fact_k, let fact_3k, let fact_6k) = _fact_chudnovsky(k, p_bits, rnd)

      let sgn = if (k %% 2) == 0 then one else minus_one end 
      //let num = sgn * fact_6k * (k_13591409 + (k_545140134 * MPFloat.from_ulong(k)))
      //let den = fact_3k * (fact_k.powi(3)) * (k_640320.pow((k_3 * MPFloat.from_ulong(k)) + k_1_5))
      //let term = num / den
      term = ((sgn * fact_6k) / fact_3k) * ((k_13591409 + (k_545140134 * MPFloat.from_ulong(k, p_bits))) / (fact_k.powi(3) * (k_640320.pow((k_3 * MPFloat.from_ulong(k, p_bits)) + k_1_5))))
      prev_res = result
      result = result + term
      k = k + 1
    until (result.almost_eq(prev_res)) or (k > p.ulong()) end
    let pi_val = (MPFloat.from_ulong(12, p_bits, rnd) * result).inv()._trunc(p_digits)

    _sign     = false
    _nan      = false
    _inf      = false
    _exponent = pi_val._exponent
    _digits   = pi_val._digits
    _rounding = rnd


  fun tag _fact_chudnovsky(n: ULong, prec: ULong, rnd: RoundingMode): (MPFloat, MPFloat, MPFloat) =>
    """
    Calculate the factorials used in Chudnovsky's algorithm:
    `n!, (3 × n)!, (6 × n)!`

    Calculations are done in `MPInt` and the results are converted to `MPFloat`
    """
    var fact_n = MPInt.from[ULong](1)
    for i in Range[ULong](1, n + 1) do
      fact_n = fact_n * MPInt.from[ULong](i)
    end

    var fact_3n = fact_n
    for i in Range[ULong](n + 1, (3 * n) + 1) do
      fact_3n = fact_3n * MPInt.from[ULong](i)
    end

    var fact_6n = fact_3n
    for i in Range[ULong]((3 * n) + 1, (6 * n) + 1) do
      fact_6n = fact_6n * MPInt.from[ULong](i)
    end

    (MPFloat.from_mpint(fact_n, prec, rnd), MPFloat.from_mpint(fact_3n, prec, rnd), MPFloat.from_mpint(fact_6n, prec, rnd))


  new val pi_bbp(prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    The Beiley-Borwein-Plouffe formula is known to allow calculation of any decimal
    of π, in 16-base, without requiring to know the previous decimals. It is also
    a very efficient formula, in O(n ln(n)).

    See https://en.wikipedia.org/wiki/Bailey%E2%80%93Borwein%E2%80%93Plouffe_formula
    See https://en.wikipedia.org/wiki/Approximations_of_pi#Efficient_methods
    """
    let p_digits: USize = (prec.usize() + (_base_bits - 1)) / _base_bits
    let p: USize = p_digits + 4   // guard digits for intermediate results
    let p_bits: ULong = (p * _base_bits).ulong()
    let k_p = p.ulong()

    let k_1 = MPFloat.from_ulong(1, p_bits, rnd)
    let k_2 = MPFloat.from_ulong(2, p_bits, rnd)
    let k_4 = MPFloat.from_ulong(4, p_bits, rnd)
    let k_5 = MPFloat.from_ulong(5, p_bits, rnd)
    let k_6 = MPFloat.from_ulong(6, p_bits, rnd)
    let k_8 = MPFloat.from_ulong(8, p_bits, rnd)
    let k_16 = MPFloat.from_ulong(16, p_bits, rnd)

    var k: ULong = 0
    var result = MPFloat(p_bits, rnd)
    var term = MPFloat(p_bits, rnd)
    var resultf: F64 = 0.0
    var termf: F64 = 0.0
    var prev_res = result
    repeat
      let t0 = k_8 * MPFloat.from_ulong(k, p_bits)
      let t0f: F64 = 8.0 * k.f64()
      Debug("Pi_bpp t0=" + t0.string() + "/ t0f=" + t0f.string())
      let t1 = k_4 / (t0 + k_1)
      let t2 = k_2 / (t0 + k_4)
      let t3 = k_1 / (t0 + k_5)
      let t4 = k_1 / (t0 + k_6)
      let t1f = 4.0 / (t0f + 1.0)
      let t2f = 2.0 / (t0f + 4.0)
      let t3f = 1.0 / (t0f + 5.0)
      let t4f = 1.0 / (t0f + 6.0)
      Debug("Pi_bpp t1=" + t1.string() + "/ t1f=" + t1f.string() +
            "; t2=" + t2.string() + "/ t2f=" + t2f.string() +
            "; t3=" + t3.string() + "/ t3f=" + t3f.string() +
            "; t4=" + t4.string() + "/ t4f=" + t4f.string())
      let t5 = k_16.powi(k.ilong(), RoundingNearest).inv()._trunc(p)
      let t5f': F64 = 16.0
      let t5f: F64 = 1.0 / t5f'.powi(k.i32())
      Debug("Pi_bbp: t5=" + t5.string() + "/ t5f=" + t5f.string())
      term = (t5 * (t1 - (t2 + t3 + t4)))._trunc(p)
      termf = (t5f * (t1f - (t2f + t3f +t4f)))
      Debug("Pi_bbp: term=" + term.string() + "/ termf=" + termf.string())
      prev_res = result
      result = (result + term)._trunc(p)
      resultf = (resultf + termf)
      Debug("Pi_bbp [" + k.string() + "]: result=" + result.string() + "/ resultf=" + resultf.string())
      k = k + 1
    until (result.almost_eq(prev_res)) or (k > k_p) end
    result

    _sign     = false
    _nan      = false
    _inf      = false
    _exponent = result._exponent
    _digits   = result._digits
    _rounding = rnd


  //- Internal helpers ---------------------------------------------------------

  fun get_precision(): ULong =>
    """
    Get the precision of the `MPFloat` mantissa (significand) in bits.

    Note: This method is kept while the `precision2` return type is not corrected in
    stdlib, for compatibility with GMP implementation.
    """
    (_digits.size() * _base_bits).ulong()


  fun get_rounding_mode(): RoundingMode =>
    """
    Get the rounding mode that was used when initializing the `MPFloat`.
    """
    _rounding


  fun _size(): USize =>
    """
    Number of base-256 digits in the internal representation.
    """
    _digits.size()


  fun _lowb(a: U16): U8 =>
    """
    Lower byte of `a`.
    """
    a.u8()


  fun _highb(a: U16): U8 =>
    """
    Upper byte of `a`.
    """
    a.shr(8).u8()


  fun _addc(a: U8, b: U8, c: U16): (U8, U16) =>
    """
    Single-digit addition `a + b + c` returning `(digit, carry)`.
    """
    let r: U16 = a.u16() + b.u16() + c
    (r.u8(), r.shr(8))


  fun _subc(a: U8, b: U8, c: U16): (U8, U16) =>
    """
    Single-digit subtraction `a − b` with borrow `c`, returning
    `(digit, new_borrow)`. Uses unsigned-complement arithmetic.
    """
    let r: U16 = ((a.max_value().u16() + a.u16()) - b.u16()) + c.shr(8)
    (r.u8(), r.shr(8))


  fun _inc_first(): MPFloat =>
    """
    Return a new `MPFloat` identical to `this` but with the most-significant
    digit (`_digits(0)`) incremented by 1. No carry propagation is performed;
    the caller must ensure `_digits(0) < 255`. Used internally for Newton
    iteration updates of the form `x ← x + 1`.
    """
    let size = _size()
    let d: Array[U8] val = recover
      let a = Array[U8].init(0, size)
      _digits.copy_to(a, 0, 0, size)
      try a.update(0, (try a(0)? else 0 end) + 1)? end
      a
    end
    MPFloat._create(_sign, _nan, _inf, _exponent, d, _rounding)


  fun _trunc(n: USize): MPFloat =>
    """
    Return a new `MPFloat` containing only the `n` most-significant base-256
    digits of `this`. Used to bound intermediate results in Newton iterations
    so that successive multiplications do not grow the array without limit.

    TODO: Take into account rounding
    """
    let s: USize = _size().min(n)
    let d: Array[U8] val = recover
      let a = Array[U8].init(0, s)
      _digits.copy_to(a, 0, 0, s)
      a
    end
    MPFloat._create(_sign, _nan, _inf, _exponent, d, _rounding)


  fun _neg_comp(): MPFloat =>
    """
    Return the byte-level one's/two's-complement negation of `_digits`. This
    is NOT a mathematical float negation — it is a fixed-point unsigned
    complement used internally by `inv()` and `sqrt()` Newton iterations
    to compute `2 − x` in unsigned arithmetic. Use `neg()` for the
    mathematical sign flip.
    """
    let size = _size()
    let d: Array[U8] val = recover
      let a = Array[U8].init(0, size)
      try
        let max: U16 = U8.max_value().u16()
        var carry: U16 = max + 1
        var k: USize = size
        repeat
          k = k - 1
          carry = (max - _digits(k)?.u16()) + _highb(carry).u16()
          a.update(k, _lowb(carry))?
        until k == 0 end
      end
      a
    end
    MPFloat._create(_sign, _nan, _inf, _exponent, d, _rounding)


  //- Predicates ---------------------------------------------------------------

  fun is_nan(): Bool =>
    """
    Return `true` if this value is Not-a-Number.
    """
    _nan


  fun is_infinite(): Bool =>
    """
    Return `true` if this value is ±∞.
    """
    _inf and not _nan


  fun is_finite(): Bool =>
    """
    Return `true` if this value is finite (not NaN and not ±∞).
    """
    not _nan and not _inf


  fun is_zero(): Bool =>
    """
    Return `true` if this value is zero (either +0 or −0).
    """
    if not is_finite() then
      return false
    end
    for d in _digits.values() do
      if d != 0 then
        return false
      end
    end
    true


  fun is_negative(): Bool =>
    """
    Return `true` if the sign bit is set (value is strictly negative or −0).
    Always `false` for NaN.
    """
    _sign and not _nan


  fun is_integer(): Bool =>
    """
    Return `true` if `this` is a finite value with no fractional part, i.e.
    `this == trunc(this)`.

    - NaN → `false`.  ±∞ → `false`.
    - ±0 → `true` (zero is an integer).
    - Any value whose base-256 representation has no non-zero fractional
      bytes (bytes at index ≥ `_exponent`) → `true`.
    - Any value with at least one non-zero fractional byte → `false`.

    Equivalent to `is_finite() and not _has_frac()`.
    """
    is_finite() and not _has_frac()


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
      return that.is_infinite() and (_sign == that._sign)
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
    let rel_part = (MPFloat.from_mpfloat(rel_tol, (p * _base_bits).ulong(), _rounding) * max_mag)._trunc(p)
    let abs_part = MPFloat.from_mpfloat(abs_tol, (p * _base_bits).ulong(), _rounding)
    let threshold = if rel_part < abs_part then abs_part else rel_part end
    diff <= threshold


  fun sign(): Compare =>
    """
    Return the sign of this value as a `Compare`:
    - `Less` for strictly negative values
    - `Greater` for strictly positive values
    - `Equal` for zero and NaN (conventional)
    """
    if _nan then
      return Equal
    end
    if is_zero() then
      return Equal
    end
    if _sign then
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
    _exponent


  fun raw_digits(): Array[U8] val =>
    """
    Return the internal base-256 mantissa as a big-endian `Array[U8]`.

    Together with `exponent()` this determines the magnitude exactly:
    `|value| = 0.d₀d₁… × 256^exponent()`.  The first byte `d₀` is always
    non-zero for finite non-zero values.

    This accessor is intended for low-level conversions such as
    `MPInt.from_mpfloat`; prefer the high-level API for ordinary use.
    """
    _digits


  fun sign_bit(): Bool =>
    """
    Return `true` when this value is strictly negative or is −0.
    NaN has an undefined sign; this method returns the raw `_sign` field.
    """
    _sign


  fun rep(): MPFRep =>
    """
    Return the pure representation of this value as an `MPFRep`.

    The `MPFRep` contains the same sign, exponent, and digit array as `this`,
    without the rounding mode.  Intended for interoperability with `MPFContext`
    and `_MPFAlgo`.
    """
    MPFRep._create(_sign, _nan, _inf, _exponent, _digits)


  fun ctx(): MPFContext =>
    """
    Return the execution context (precision + rounding mode) of this value
    as an `MPFContext`.

    The precision is derived from `_digits.size() × 8` bits; the rounding
    mode is `_rounding`.
    """
    MPFContext((_size() * 8).ulong(), _rounding)


  //- Arithmetic ---------------------------------------------------------------

  fun _cmp_mag(that: MPFloat): Compare =>
    """
    Compare the magnitudes `|this|` and `|that|`. Returns `Greater` if
    `|this| > |that|`, `Less` if `|this| < |that|`, `Equal` otherwise.
    Both operands must be finite and non-zero; callers must check.
    """
    if _exponent > that._exponent then
      return Greater
    end
    if _exponent < that._exponent then
      return Less
    end
    let na: USize = _size()
    let nb: USize = that._size()
    let n: USize = na.min(nb)
    var i: USize = 0
    try
      while i < n do
        let ai = _digits(i)?
        let bi = that._digits(i)?
        if ai > bi then
          return Greater
        end
        if ai < bi then
          return Less
        end
        i = i + 1
      end
    end
    // After matching the common prefix, check whether the extra bytes of
    // the longer operand are all zero.  If so the magnitudes are equal.
    if na > nb then
      var all_zero: Bool = true
      try
        while i < na do
          if _digits(i)? != 0 then all_zero = false; break end
          i = i + 1
        end
      end
      if all_zero then Equal else Greater end
    elseif na < nb then
      var all_zero: Bool = true
      try
        while i < nb do
          if that._digits(i)? != 0 then all_zero = false; break end
          i = i + 1
        end
      end
      if all_zero then Equal else Less end
    else
      Equal
    end


  fun _add_mag(that: MPFloat, sgn: Bool): MPFloat =>
    """
    Add the magnitudes `|this|` and `|that|` and return the sum with the
    given `sgn`. Both operands must be finite. The result exponent equals
    `max(this._exponent, that._exponent)`, incremented by one when a carry
    propagates out of the most-significant column.

    Operands are aligned on their most-significant digit before addition.
    The result precision equals `max(|this|_digits, shift + |that|_digits)`
    where `shift` is the exponent difference, so no precision is lost.

    TODO: Complete rounding propagation
    """
    // Align both operands
    (let ea, let eb, let ad, let bd, let na, let nb) =
      if _exponent >= that._exponent then
        (_exponent, that._exponent, _digits, that._digits, _size(), that._size())
      else
        (that._exponent, _exponent, that._digits, _digits, that._size(), _size())
      end
    let shift: USize = (ea - eb).usize()

    // When b falls entirely below a's precision, the sum equals a.
    if shift >= na then
      return MPFloat._create(sgn, false, false, ea, ad, _rounding)
    end

    // One guard digit at index 0 absorbs a possible carry out of column 0.
    let prec: USize = na.max(shift + nb)
    let result_size: USize = prec + 1
    let raw: Array[U8] val = recover
      let res = Array[U8].init(0, result_size)
      try
        var carry: U16 = 0
        var col: USize = prec
        repeat
          col = col - 1
          let ai: U8 = if col < na then ad(col)? else 0 end
          let bi: U8 = if (col >= shift) and ((col - shift) < nb) then
              bd(col - shift)?
            else
              0
            end
          (let sum, let c2) = _addc(ai, bi, carry)
          carry = c2
          res.update(col + 1, sum)?
        until col == 0 end
        res.update(0, _lowb(carry))?
      end
      res
    end

    let carry_digit: U8 = try raw(0)? else 0 end
    if carry_digit == 0 then
      // No carry: strip the guard digit, exponent stays ea.
      let d: Array[U8] val = recover
        let a2 = Array[U8].init(0, prec)
        raw.copy_to(a2, 1, 0, prec)
        a2
      end
      MPFloat._create(sgn, false, false, ea, d, RoundingNearest)
    else
      // Carry out: one extra leading digit, exponent becomes ea + 1.
      MPFloat._create(sgn, false, false, ea + 1, raw, RoundingNearest)
    end


  fun _sub_mag(that: MPFloat, sgn: Bool): MPFloat =>
    """
    Subtract the magnitude `|that|` from `|this|`, where `|this| >= |that|`
    (enforced by the caller via `_cmp_mag`). Returns the difference with the
    given `sgn`. Both operands must be finite. The result is normalised:
    leading zero bytes are stripped and the exponent adjusted accordingly.

    TODO: Complete rounding propagation
    """
    let ea: I64 = _exponent
    let eb: I64 = that._exponent
    let shift: USize = (ea - eb).usize()
    let na: USize = _size()
    let nb: USize = that._size()
    let result_size: USize = na.max(shift + nb)
    // carry.shr(8) == 1 means no borrow; == 0 means borrow — same convention
    // as _neg_comp. carry is maintained as a full U16 to preserve the bit.
    let max_b: U16 = U8.max_value().u16()
    let raw: Array[U8] val = recover
      let res = Array[U8].init(0, result_size)
      try
        var carry: U16 = max_b + 1
        var col: USize = result_size
        repeat
          col = col - 1
          let ai: U16 = if col < na then _digits(col)?.u16() else 0 end
          let bi: U16 = if (col >= shift) and ((col - shift) < nb) then
              that._digits(col - shift)?.u16()
            else
              0
            end
          carry = ((max_b + ai) - bi) + _highb(carry).u16()
          res.update(col, _lowb(carry))?
        until col == 0 end
      end
      res
    end
    
    // Normalise: find the first non-zero byte and adjust the exponent.
    let raw_size: USize = raw.size()
    var leading: USize = 0
    try
      while (leading < raw_size) and (raw(leading)? == 0) do
        leading = leading + 1
      end
    end
    if leading == raw_size then
      return MPFloat.create((na * _base_bits).ulong())
    end

    let new_size: USize = raw_size - leading
    let new_exp: I64 = ea - leading.i64()
    let d: Array[U8] val = recover
      let a2 = Array[U8].init(0, new_size)
      raw.copy_to(a2, leading, 0, new_size)
      a2
    end
    MPFloat._create(sgn, false, false, new_exp, d, RoundingNearest)


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

    TODO: Check rounding propagation
    """
    if _nan or that._nan then
      return MPFloat.nan_val()
    end
    if _inf and that._inf then
      if _sign == that._sign then
        return MPFloat._create(_sign, false, true, 0, Array[U8].create(), RoundingNearest)
      end
      return MPFloat.nan_val()
    end
    if _inf then
      return MPFloat._create(_sign, false, true, 0, Array[U8].create(), RoundingNearest)
    end
    if that._inf then
      return MPFloat._create(that._sign, false, true, 0, Array[U8].create(), RoundingNearest)
    end
    if is_zero() then
      return MPFloat._create(that._sign, false, false, that._exponent, that._digits, that._rounding)
    end
    if that.is_zero() then
      return MPFloat._create(_sign, false, false, _exponent, _digits, _rounding)
    end
    if _sign == that._sign then
      return _add_mag(that, _sign)
    end
    match _cmp_mag(that)
    | Greater => _sub_mag(that, _sign)
    | Less    => that._sub_mag(MPFloat._create(_sign, _nan, _inf, _exponent, _digits, _rounding), that._sign)
    else
      MPFloat.create((_digits.size() * _base_bits).ulong(), _rounding)
    end


  fun sub(that: MPFloat): MPFloat =>
    """
    Calculate `this − that`. Delegates to `add` after negating `that`,
    so all sign, exponent, and special-value handling is inherited.
    """
    add(that.neg())


  fun _short_add(b: U8): MPFloat =>
    """
    Short addition: add the single byte `b` to the least-significant digit of
    `this`, propagating carry toward the most-significant digit.
    """
    let size: USize = _size()
    var err_i: USize = 0
    let d: Array[U8] val = recover
      let res = Array[U8].init(0, size)
      var i: USize = size
      try
        var carry: U16 = b.u16().shl(8)
        repeat
          i = i - 1
          err_i = i
          carry = _digits(i)?.u16() + _highb(carry).u16()
          if (i + 1) < size then
            res.update(i + 1, _lowb(carry))?
          end
        until i == 0 end
        res.update(0, _lowb(carry))?
      else
        try
          Assert(false, "[MPFloat._short_add] index out of bounds at i=" + err_i.string(),
            true)?
        end
      end
      res
    end
    MPFloat._create(_sign, _nan, _inf, _exponent, d, _rounding)


  fun _short_mul(b: U8): MPFloat =>
    """
    Short multiplication: multiply `this` by the single byte `b`.
    The result has the same digit count as `this`. Any overflow beyond the
    most-significant digit is captured in `_digits(0)` via the carry.
    """
    let size: USize = _size()
    var err_i: USize = 0
    let d: Array[U8] val = recover
      let res = Array[U8].init(0, size)
      var i: USize = size
      try
        var carry: U16 = 0
        repeat
          i = i - 1
          err_i = i
          carry = (_digits(i)?.u16() * b.u16()) + _highb(carry).u16()
          if i < (size - 1) then
            res.update(i + 1, _lowb(carry))?
          end
        until i == 0 end
        res.update(0, _highb(carry))?
      else
        try
          Assert(false, "[MPFloat._short_mul] index out of bounds at i=" + err_i.string(),
            true)?
        end
      end
      res
    end
    MPFloat._create(_sign, _nan, _inf, _exponent, d, _rounding)


  fun _short_div(b: U8): (MPFloat, U8) =>
    """
    Short division: divide `this` by the single byte `b`. Returns
    `(quotient, remainder)`. The quotient has the same digit count as `this`.
    """
    let size: USize = _size()
    var remain: U16 = 0
    var err_i: USize = 0
    let d: Array[U8] val = recover
      let res = Array[U8].init(0, size)
      var i: USize = 0
      try
        while i < size do
          err_i = i
          (let quotient, let r2) = (remain.shl(8) + _digits(i)?.u16()).divrem(b.u16())
          remain = r2
          res.update(i, _lowb(quotient))?
          i = i + 1
        end
      else
        try
          Assert(false, "[MPFloat._short_div] index out of bounds at i=" + err_i.string(),
            true)?
        end
      end
      res
    end
    (MPFloat._create(_sign, _nan, _inf, _exponent, d, _rounding), remain.u8())


  fun neg(): MPFloat =>
    """
    Return the arithmetic negation: a new `MPFloat` identical to `this` but
    with the sign flipped. NaN is returned unchanged.
    """
    if _nan then
      return MPFloat.nan_val()
    end
    MPFloat._create(not _sign, _nan, _inf, _exponent, _digits, _rounding)


  fun digit_shl(n: USize = 1): MPFloat =>
    """
    Left-shift by `n` base-256 digit positions: drop the `n` most-significant
    digits and zero-pad on the right. Equivalent to multiplying by `256^n`
    modulo the array size. The exponent field is not adjusted.

    Used internally to strip the integer part after an arithmetic step.
    """
    let size: USize = _size()
    let s: USize = if size > n then size - n else 0 end
    let d: Array[U8] val = recover
      let a = Array[U8].init(0, size)
      var i: USize = 0
      try
        while i < s do
          a.update(i, _digits(i + n)?)?
          i = i + 1
        end
      else
        try
          Assert(false, "[MPFloat.digit_shl] index out of bounds at i=" + i.string(), true)?
        end
      end
      a
    end
    MPFloat._create(_sign, _nan, _inf, _exponent, d, _rounding)


  fun mul(that: MPFloat): MPFloat =>
    """
    Fast multiplication of `this` and `that` using FFT-based convolution.

    Uses real FFT (Cooley-Tukey); limited to arrays of ~10⁶ base-256 digits.
    For larger sizes, a Number-Theoretic Transform would be required.

    TODO: Complete rounding propagation
    """
    // Special values
    if _nan or that._nan then
      return MPFloat.nan_val()
    end
    if _inf then
      if that.is_zero() then
        return MPFloat.nan_val()
      end
      return MPFloat._create((_sign != that._sign), false, true, 0, Array[U8].create(), RoundingNearest)
    end
    if that._inf then
      if is_zero() then
        return MPFloat.nan_val() end
      return MPFloat._create((_sign != that._sign), false, true, 0, Array[U8].create(), RoundingNearest)
    end
    // TODO: better rounding estimation
    let new_rnd: RoundingMode = _rounding
    if is_zero() or that.is_zero() then
      return MPFloat.create(0, new_rnd)
    end
    mul_unsafe(that)


  fun inv(): MPFloat =>
    """
    Calculate `1 / this` using Newton's method:  `x_{n+1} = x_n × (2 − G × x_n)`
    where `G = 256 × |this|` (the value of `this` reinterpreted with exponent 1).

    `G ∈ [1, 256)`, so the iteration converges to `y = 1/G ∈ (1/256, 1]`.
    The final result exponent is `(res._exponent + 1) − _exponent`.

    Sign propagation: `1/positive = positive`, `1/negative = negative`.
    A safety limit of `size × 4` iterations prevents non-convergence.

    TODO: Rounding management
    """
    // The special values
    if _nan then
      return MPFloat.nan_val()
    end
    if _inf then
      return MPFloat.create(0, _rounding)
    end
    if is_zero() then
      return MPFloat.inf_val(not _sign)
    end
    inv_unsafe()


  fun div(that: MPFloat): MPFloat =>
    """
    Calculate `this / that` as `this × (1 / that)`.

    Special cases follow IEEE 754:
    - NaN in either operand → NaN
    - ±∞ / ±∞ → NaN;  finite / ±∞ → ±0;  ±∞ / finite → ±∞
    - finite / 0 → ±∞ (sign = XOR of operand signs);  0 / 0 → NaN

    TODO: Complete rounding propagation
    """
    // Special values
    if _nan or that._nan then
      return MPFloat.nan_val()
    end
    if _inf and that._inf then
      return MPFloat.nan_val()
    end
    if _inf then
      return MPFloat._create((_sign != that._sign), false, true, 0, Array[U8].create(), RoundingNearest)
    end
    // TODO: Better rounding evaluation
    let new_rnd = RoundingNearest
    if that._inf then
      return MPFloat.create(0, new_rnd)
    end
    if that.is_zero() then
      if is_zero() then
        return MPFloat.nan_val()
      end
      return MPFloat._create((_sign != that._sign), false, true, 0, Array[U8].create(), new_rnd)
    end
    if is_zero() then
      return MPFloat.create(0, _rounding)
    end
    mul_unsafe(that.inv_unsafe())


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

    Special cases follow IEEE 754:
    - NaN → NaN
    - +∞ → +∞; −∞ → NaN
    - ±0 → ±0 (sign preserved)
    - negative non-zero → NaN

    A safety limit of `size × 4` iterations prevents non-convergence.

    TODO: Complete rounding propagation
    """
    // Special values
    if _nan then
      return MPFloat.nan_val()
    end
    if _inf then
      if _sign then
        return MPFloat.nan_val()
      end
      return MPFloat.inf_val(true)
    end
    if is_zero() then
      return MPFloat._create(_sign, false, false, 0, _digits, _rounding)
    end
    if _sign then
      return MPFloat.nan_val()
    end
    sqrt_unsafe()


  //- String conversion -------------------------------------------------------

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
    - NaN → `("nan", 0, false)`, ±∞ → `("±inf", 0, false)`, ±0 → `("0", 0, false)`

    Only base 10 is currently implemented; other bases return `("", 0, true)`.
    The `inexact` flag is always `false` (rounding is a future TODO).

    Algorithm (exact, no Newton approximation):
    Let `N` be `_digits` interpreted as a base-256 integer and
    `k = _exponent − _size()`.  Then `|this| = N × 256^k`.

    - Integer case (`k ≥ 0`): `|this| = N × 2^{8k}`.
      Compute `N << (8k)` as an `MPInt` and call `.string()`.
    - Fractional case (`k < 0`): `|this| = N × 5^b / 10^b` where `b = 8|k|`.
      Compute `(N × 5^b).string()` via `MPInt` and set
      `dec_exp = len(str) − b`.

    Both paths are always exact — no scaling error, no Newton undershooting.
    """
    // Special values
    if _nan then
      return ("nan".clone(), 0, false)
    end
    if _inf then
      let s: String iso = if _sign then "-inf".clone() else "+inf".clone() end
      return (consume s, 0, false)
    end
    if is_zero() then
      return ("0".clone(), 0, false)
    end
    if base != 10 then
      return ("".clone(), 0, true)
    end

    // Algorithm: value = N × 256^{k_bytes} = N × 2^b_bits  (or N / 2^b_bits)
    // where N = _digits interpreted as a base-256 integer (all bytes integer),
    // and k_bytes = _exponent − _size().
    //
    // Integer case (k_bytes ≥ 0): value = N × 2^{8 × k_bytes}.
    //   Decimal string = (N << (8 × k_bytes)).string(), dec_exp = len(str).
    //
    // Fractional case (k_bytes < 0): value = N / 2^b where b = 8 × |k_bytes|.
    //   = N × 5^b / 10^b.
    //   Decimal mantissa = (N × 5^b).string(), dec_exp = len − b.
    //
    // Both paths are exact (no Newton approximation).  MPInt handles the
    // big-integer arithmetic; _short_div / scaling tricks are no longer needed.

    let prec: USize = _size()
    let n_dec: USize = ((prec.f64() * 2.41).usize() + 2).max(1)
    let sign_prefix: String = if _sign then "-" else "" end

    // Extract N as an MPInt: a synthetic MPFloat with _exponent = prec treats
    // all bytes as integer bytes.  MPInt.from_mpfloat truncates toward zero,
    // which for a pure-integer synthetic value gives exactly N.
    let n_float = MPFloat._create(false, false, false, prec.i64(), _digits, _rounding)
    let n_int = try MPInt.from_mpfloat(n_float)? else MPInt.from[ILong](0) end

    let k_bytes: I64 = _exponent - prec.i64()

    if (k_bytes >= 0) and (k_bytes < 150) then
      // Integer case: value = N × 2^{8 × k_bytes}
      let b: USize = k_bytes.usize() * _base_bits
      let shifted: MPInt = if b > 0 then
          n_int.shl(MPInt.from[ILong](b.ilong()))
        else
          n_int
        end
      let dec_str_iso: String iso = shifted.string()
      let dec_exp: I64 = dec_str_iso.size().i64()
      let total: USize = n_dec.max(dec_str_iso.size())
      let result: String iso = recover
        let s = String.create(sign_prefix.size() + total)
        s.append(sign_prefix)
        s.append(consume dec_str_iso)
        var i: USize = dec_exp.usize()
        while i < n_dec do
          s.push('0')
          i = i + 1
        end
        s
      end
      (consume result, dec_exp, false)
    elseif (k_bytes < 0) and (k_bytes > -150) then
      // Fractional case: value = N / 2^b = N × 5^b / 10^b
      let b: USize = ((-k_bytes) * 8).usize()
      let five_pow: MPInt = MPInt.from[ILong](5).pow(MPInt.from[ILong](b.ilong()))
      let numerator: MPInt = n_int.mul(five_pow)
      let num_str_iso: String iso = numerator.string()
      let num_len: USize = num_str_iso.size()
      let dec_exp: I64 = num_len.i64() - b.i64()
      let total: USize = num_len.max(n_dec)
      let result: String iso = recover
        let s = String.create(sign_prefix.size() + total)
        s.append(sign_prefix)
        s.append(consume num_str_iso)
        var i: USize = num_len
        while i < n_dec do
          s.push('0')
          i = i + 1
        end
        s
      end
      (consume result, dec_exp, false)
    else
      // Extreme values: use approximate path via F64 to avoid OOM in MPInt/String
      ifdef debug then
        Debug("[MPFloat.exact_string] Extreme value encountered (k_bytes=" + k_bytes.string() + "). Exact decimal conversion range exceeded; using approximate F64 path.")
      end

      let f = f64()
      if f.nan() then
        return ("nan".clone(), 0, true)
      end
      if f.infinite() then
        return ("1".clone(), 1000, true)
      end
      if f == 0 then
        return ("0".clone(), -1000, true)
      end
      
      let s_f = f.abs().string()
      
      // Basic parsing of F64 string to satisfy the contract
      // e.g. "1.234e+45" or "1.234"
      try
        let e_pos = s_f.find("e")?
        let m_str: String val = s_f.substring(0, e_pos)
        let e_str: String val = s_f.substring(e_pos + 1)
        
        var m_raw = recover String end
        var dot_seen = false
        var dot_pos: USize = 0
        var idx: USize = 0
        for char in m_str.values() do
          if char == '.' then
            dot_seen = true
            dot_pos = idx
          else
            m_raw.push(char)
          end
          idx = idx + 1
        end
        
        let exp_val = e_str.i64()? + (if dot_seen then dot_pos.i64() else m_str.size().i64() end)
        (consume m_raw, exp_val, true)
      else
        // No 'e', just a simple float string like "123.456"
        let s_f_val: String val = consume s_f
        var m_raw = recover String end
        var dot_seen = false
        var dot_pos: USize = 0
        var idx: USize = 0
        for char in s_f_val.values() do
          if char == '.' then
            dot_seen = true
            dot_pos = idx
          else
            m_raw.push(char)
          end
          idx = idx + 1
        end
        let exp_val = if dot_seen then dot_pos.i64() else s_f_val.size().i64() end
        (consume m_raw, exp_val, true)
      end
    end


  fun string(): String iso^ =>
    """
    Return a decimal string representation of this `MPFloat`, formatted in
    the style of C's `printf` `%g` with `⌈_size() × log₁₀(256)⌉ + 2` significant
    digits. Trailing zeros are stripped; at least one digit after the decimal
    point is always shown.

    Special values: `"nan"`, `"+inf"`, `"-inf"`, `"0.0"`, `"-0.0"`.

    Scientific notation (`d.dddde±N`) is used when the decimal exponent is
    less than `−4` or greater than or equal to the number of significant
    digits; fixed notation (`ddd.ddd`) is used otherwise.
    """
    // The special cases
    if _nan then
      return "nan".clone()
    end
    if _inf then
      return (if _sign then "-inf" else "+inf" end).clone()
    end
    if is_zero() then
      return (if _sign then "-0.0" else "0.0" end).clone()
    end

    (let raw_mant, let dec_exp, _) = exact_string(10, _rounding)
    let raw: String = consume raw_mant

    // Separate the optional leading '-' from the digit characters.
    var dig_start: USize = 0
    let sgn: String =
      if (raw.size() > 0) and (try raw(0)? == '-' else false end) then
        dig_start = 1
        "-"
      else
        ""
      end

    // Strip trailing zeros but keep at least one significant digit.
    let all_digits: String = raw.trim(dig_start)
    var sig_end: USize = all_digits.size()
    while (sig_end > 1) and (try all_digits(sig_end - 1)? == '0' else false end) do
      sig_end = sig_end - 1
    end
    let digits: String = all_digits.trim(0, sig_end)
    let nd: USize = digits.size()

    // %g rule: scientific if dec_exp < -4 or dec_exp >= full precision digit count.
    let use_sci: Bool = (dec_exp < -4) or (dec_exp >= all_digits.size().i64())

    recover
      let s = String.create(sgn.size() + nd + 8)
      s.append(sgn)
      if use_sci then
        // d.dddde±EE
        s.push(try digits(0)? else '0' end)
        if nd > 1 then
          s.append(".")
          s.append(digits.trim(1))
        end
        s.append("e")
        let exp_str: String = (dec_exp - 1).string()
        if dec_exp >= 1 then s.append("+") end
        s.append(exp_str)
      else
        let e: ISize = dec_exp.isize()
        let ni: ISize = nd.isize()
        if e <= 0 then
          // 0.000digits
          s.append("0.")
          var zi: ISize = 0
          while zi < -e do s.push('0'); zi = zi + 1 end
          s.append(digits)
        elseif e >= ni then
          // digits followed by trailing zeros, then ".0"
          s.append(digits)
          var zi: ISize = e - ni
          while zi > 0 do s.push('0'); zi = zi - 1 end
          s.append(".0")
        else
          // ddd.ddd
          s.append(digits.trim(0, e.usize()))
          s.append(".")
          s.append(digits.trim(e.usize()))
        end
      end
      s
    end


  //- Unsafe arithmetic -------------------------------------------------------

  fun add_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe addition of `this` and `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    if _sign == that._sign then
      return _add_mag(that, _sign)
    end
    match _cmp_mag(that)
    | Greater => _sub_mag(that, _sign)
    | Less    => that._sub_mag(MPFloat._create(_sign, _nan, _inf, _exponent, _digits, _rounding), that._sign)
    else
      MPFloat.create(get_precision(), _rounding)
    end

  
  fun sub_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe substraction of `that` from `this`.If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    add_unsafe(that.neg())


  fun mul_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe multiplication of `this` and `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.

    TODO: Complete rounding propagation
    """

    let result_sign: Bool = _sign != that._sign
    let this_size: USize = _size()
    let that_size: USize = that._size()
    // res_size = na + nb holds carry digit + all na+nb-1 convolution digits.
    let res_size: USize = this_size + that_size
    // pow2 must be a power of 2 >= res_size (no circular aliasing) and >= 4
    // because fourier_real accesses index (size/2)+1 which requires size >= 4.
    let pow2: USize = res_size.next_pow2().max(4)

    var err_i: USize = 0
    let d: Array[U8] val = recover
      let res = Array[U8].init(0, res_size)
      var i: USize = 0
      try
        let a = Array[F64].init(0.0, pow2)
        i = 0
        while i < this_size do
          a.update(i, _digits(i)?.f64())?
          i = i + 1
        end

        let b = Array[F64].init(0.0, pow2)
        i = 0
        while i < that_size do
          b.update(i, that._digits(i)?.f64())?
          i = i + 1
        end

        FFT.fourier_real(a)
        FFT.fourier_real(b)

        // Pointwise complex multiply in the frequency domain.
        b.update(0, b(0)? * a(0)?)?
        b.update(1, b(1)? * a(1)?)?
        i = 2
        while i < pow2 do
          let temp = b(i)?
          b.update(i, (temp * a(i)?) - (b(i + 1)? * a(i + 1)?))?
          b.update(i + 1, (temp * a(i + 1)?) + (b(i + 1)? * a(i)?))?
          i = i + 2
        end

        // fourier_real inverse already normalises: no extra s2 factor needed.
        FFT.fourier_real(b, true)

        // Convert back to U8 digits with carry propagation (high to low).
        var carry: F64 = 0.0
        let base = _base.f64()
        i = pow2
        repeat
          i = i - 1
          err_i = i
          let temp = b(i)? + carry + 0.5
          carry = (temp / base).u64().f64()
          b.update(i, temp - (carry * base))?
        until i == 0 end

        // res[0] = carry (overflow); res[1..res_size-1] = b[0..res_size-2].
        res.update(0, carry.u8())?
        i = 1
        while i < res_size do
          res.update(i, b(i - 1)?.u8())?
          i = i + 1
        end
      else
        try
          Assert(false, "[MPFloat.mul] index out of bounds at i=" + err_i.string(), true)?
        end
      end
      res
    end

    // Strip leading zeros; each stripped zero shifts the exponent down by one.
    var leading: USize = 0
    try
      while (leading < d.size()) and (d(leading)? == 0) do
        leading = leading + 1
      end
    end
    if leading == d.size() then
      return MPFloat.create(0)
    end
    let new_exp: I64 = (_exponent + that._exponent) - leading.i64()
    let new_size: USize = d.size() - leading
    let norm_d: Array[U8] val = recover
      let a2 = Array[U8].init(0, new_size)
      d.copy_to(a2, leading, 0, new_size)
      a2
    end
    // TODO: better rounding estimation
    let new_rnd: RoundingMode = _rounding
    MPFloat._create(result_sign, false, false, new_exp, norm_d, new_rnd)


  fun _trunc_frac(): MPFloat =>
    """
    Internal alias for `trunc()`.  Called by `divrem`, `fld_unsafe`, and
    `divrem_unsafe` to obtain the truncated-toward-zero integer part without
    going through the public API.
    """
    trunc()


  fun _has_frac(): Bool =>
    """
    Return `true` if `this` has at least one non-zero fractional base-256 byte,
    i.e. `this ≠ trunc(this)`.

    - NaN or ±∞ → `false` (no fractional part is defined for non-finite values).
    - Zero → `false`.
    - Purely fractional (`_exponent ≤ 0`) → `true` for any non-zero value.
    - All-integer (`_exponent ≥ _size()`) → `false`.
    - Mixed → scan bytes at index `_exponent` and beyond for any non-zero.
    """
    if not is_finite() then
      return false
    end
    if is_zero() then
      return false
    end
    if _exponent <= 0 then
      return true
    end
    let e: USize = _exponent.usize()
    if e >= _size() then
      return false
    end
    var i: USize = e
    try
      while i < _size() do
        if _digits(i)? != 0 then
          return true
        end
        i = i + 1
      end
    end
    false


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
    let nan = MPFloat.nan_val()
    if _nan or that._nan then
      return (nan, nan)
    end
    if _inf and that._inf then
      return (nan, nan)
    end
    if _inf then
      return (MPFloat._create(_sign != that._sign, false, true, 0, Array[U8].create(), _rounding), nan)
    end
    if that._inf then
      return (MPFloat.create(get_precision(), _rounding), MPFloat._create(_sign, false, false, _exponent, _digits, _rounding))
    end
    if that.is_zero() then
      if is_zero() then
        return (nan, nan)
      end
      return (MPFloat._create(_sign != that._sign, false, true, 0, Array[U8].create(), _rounding), nan)
    end
    if is_zero() then
      let z = MPFloat.create(get_precision(), _rounding)
      return (z, z)
    end
    var q: MPFloat = div(that)._trunc_frac()
    var r: MPFloat = sub(q.mul(that))
    // Post-correction: Newton's inv() always undershoots, so q may be 1 too
    // small when this/that is an exact integer (e.g. 10/5 → q≈1.999 → trunc=1).
    // If |r| ≥ |that|, increment |q| by 1 and adjust r.  At most one step is
    // needed because the Newton error is bounded by one base-256 ULP.
    if (not r.is_zero()) and (not r.abs().lt(that.abs())) then
      let one = MPFloat.from_f64(1.0, get_precision(), _rounding)
      if _sign == that._sign then
        // Positive quotient: increase q toward +∞.
        q = q.add(one)
      else
        // Negative quotient: decrease q toward −∞.
        q = q.sub(one)
      end
      r = sub(q.mul(that))
    end
    (q, r)


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
    (let q, let r) = divrem(that)
    if r.is_nan() then
      return MPFloat.nan_val()
    end
    if r._has_frac() or (not r.is_zero()) then
      if _sign != that._sign then
        return q.sub(MPFloat.from_f64(1.0, get_precision(), _rounding))
      end
    end
    q


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
    let r: MPFloat = rem(that)
    if r.is_nan() then
      return MPFloat.nan_val()
    end
    if (not r.is_zero()) and (_sign != that._sign) then
      return r.add(that)
    end
    r


  //- Unsafe operations -------------------------------------------------------

  fun div_unsafe(that: MPFloat): MPFloat =>
    """
    Unsafe division of `this` by `that`. If any input or output of the operation
    is +/- infinity or NaN, the result is undefined.
    """
    MPFloat._create(_sign, _nan, _inf, _exponent, _digits, _rounding).div(that)


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
    let prec: USize = _size()
    // One extra digit guards Newton rounding.
    let size: USize = prec + 1

    // G = same digits as |this| but with exponent=1, so its integer part is
    // d[0] ∈ [1, 256).  Newton converges to y = 1/G ∈ (1/256, 1].
    let g = MPFloat._create(false, false, false, 1, _digits, _rounding)

    // F64 initial estimate: fg = d[0] + d[1]/256 + … ≈ G.
    let ng: USize = prec.min(4)
    var fg: F64 = 0.0
    try
      var i: USize = ng
      repeat
        i = i - 1
        fg = (fg / 256.0) + _digits(i)?.f64()
      until i == 0 end
    end
    var res = MPFloat.from_f64(1.0 / fg, (size * _base_bits).ulong(), _rounding)

    // Newton refinement: y_{n+1} = y_n × (2 − G × y_n).
    let two = MPFloat.from_f64(2.0, (size * _base_bits).ulong(), _rounding)
    var iters: USize = 0
    let max_iters: USize = size * 4
    while iters < max_iters do
      iters = iters + 1
      let gy: MPFloat = g.mul(res)._trunc(size)
      let delta: MPFloat = two.sub(gy)._trunc(size)
      let new_res: MPFloat = res.mul(delta)._trunc(size)

      // Converged when all leading `prec` digits agree with previous iterate.
      var changed: Bool = false
      try
        let ns: USize = new_res._size().min(res._size()).min(prec)
        var ci: USize = 0
        while ci < ns do
          if new_res._digits(ci)? != res._digits(ci)? then
            changed = true
            break
          end
          ci = ci + 1
        end
      end
      res = new_res
      if not changed then
        break
      end
    end

    // TODO: Evaluate rounding
    let new_rnd = _rounding
    // 1/|A| = res × 256^{1−e}.  res._exponent accounts for any F64-init shift.
    MPFloat._create(_sign, false, false, (res._exponent + 1) - _exponent, res._digits, new_rnd)


  fun fld_unsafe(that: MPFloat): MPFloat =>
    """
    Floored division of `this` by `that`. If any input or output of the
    operation is +/- infinity or NaN, the result is undefined.

    Computes `trunc(this / that)` without special-case guards, then adjusts
    downward by 1 if the signs of `this` and `that` differ and the remainder
    is non-zero.
    """
    let q: MPFloat = div_unsafe(that)._trunc_frac()
    let r: MPFloat = sub_unsafe(q.mul_unsafe(that))
    if (not r.is_zero()) and (_sign != that._sign) then
      q.sub_unsafe(MPFloat.from_f64(1.0, get_precision(), _rounding))
    else
      q
    end


  fun divrem_unsafe(that: MPFloat): (MPFloat, MPFloat) =>
    """
    Truncated division with remainder of `this` by `that`. If any input or
    output of the operation is +/- infinity or NaN, the result is undefined.

    Returns `(q, r)` where `q = trunc(this /~ that)` and `r = this -~ q *~ that`.
    """
    let q: MPFloat = div_unsafe(that)._trunc_frac()
    let r: MPFloat = sub_unsafe(q.mul_unsafe(that))
    (q, r)


  fun rem_unsafe(that: MPFloat): MPFloat =>
    """
    Truncated remainder of the division of `this` by `that`. If any input or
    output of the operation is +/- infinity or NaN, the result is undefined.

    Returns `this -~ trunc(this /~ that) *~ that`.
    """
    divrem_unsafe(that)._2


  fun mod_unsafe(that: MPFloat): MPFloat =>
    """
    Floored remainder of `this` by `that`. If any input or output of the
    operation is +/- infinity or NaN, the result is undefined.

    Returns `this -~ fld_unsafe(this, that) *~ that`.  Equivalent to `rem_unsafe` when
    `this` and `that` have the same sign; adds `that` to the truncated remainder
    when the signs differ and the remainder is non-zero.
    """
    let r: MPFloat = rem_unsafe(that)
    if (not r.is_zero()) and (_sign != that._sign) then
      r.add_unsafe(that)
    else
      r
    end


  fun neg_unsafe(): MPFloat =>
    """
    Change the sign of `this`. If `this` is NaN, the result is undefined.
    """
    MPFloat._create(not _sign, _nan, _inf, _exponent, _digits, _rounding)
  

  fun sqrt_unsafe(): MPFloat =>
    """
    Calculate the square root of `this` using Newton's method for the
    reciprocal square root: `y_{n+1} = y_n × (3 − H × y_n²) / 2`. If `this`
    if +/- infinity or NaN, the result is undefined.

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

    TODO: Complete rounding propagation
    """
    let prec: USize = _size()
    // One extra digit guards Newton rounding.
    let size: USize = prec + 1

    // Choose h_exp so that _exponent − h_exp is always even.
    // Odd _exponent → h_exp=1, H ∈ [1, 256).
    // Even _exponent → h_exp=2, H ∈ [256, 65536).
    let h_exp: I64 = if (_exponent and 1) != 0 then 1 else 2 end
    let h = MPFloat._create(false, false, false, h_exp, _digits, _rounding)

    // F64 initial estimate: fg ≈ d0 + d1/256 + … ≈ H / 256^(h_exp−1).
    let ng: USize = prec.min(4)
    var fg: F64 = 0.0
    let base = _base.f64()
    try
      var k: USize = ng
      repeat
        k = k - 1
        fg = (fg / base) + _digits(k)?.f64()
      until k == 0 end
    end
    let fg_h: F64 = if h_exp == 1 then fg else fg * base end
    var res = MPFloat.from_f64(1.0 / fg_h.sqrt(), (size * _base_bits).ulong(), _rounding)

    // Newton refinement: y_{n+1} = y_n × (3 − H × y_n²) / 2.
    let three = MPFloat.from_f64(3.0, (size * _base_bits).ulong(), _rounding)
    var iters: USize = 0
    let max_iters: USize = size * 4
    while iters < max_iters do
      iters = iters + 1
      let y2: MPFloat = res.mul(res)._trunc(size)
      let hy2: MPFloat = h.mul(y2)._trunc(size)
      let delta: MPFloat = three.sub(hy2)._trunc(size)
      let new_res_full: MPFloat = res.mul(delta)._trunc(size)
      (let new_res_halved, _) = new_res_full._short_div(2)
      let nr: MPFloat = new_res_halved._trunc(size)

      // Converged when all leading `prec` digits agree with previous iterate.
      var changed: Bool = false
      try
        let ns: USize = nr._size().min(res._size()).min(prec)
        var i: USize = 0
        while i < ns do
          if nr._digits(i)? != res._digits(i)? then
            changed = true
            break
          end
          i = i + 1
        end
      end
      res = nr
      if not changed then
        break
      end
    end

    // √H = H × (1/√H); exponent adjusted for the parity split.
    let sqrt_h: MPFloat = h.mul(res)._trunc(size)
    let result_exp: I64 = sqrt_h._exponent + ((_exponent - h_exp) / 2)
    // TODO: Evaluate rounding
    let new_rnd = _rounding
    MPFloat._create(false, false, false, result_exp, sqrt_h._digits, new_rnd)


  //- Comparisons -------------------------------------------------------------

  fun eq(that: MPFloat): Bool =>
    """
    Equality of `this` and `that`.

    * NaN is never equal to anything, including itself
    * Infinities of the same sign are equal; `+∞ ≠ −∞`
    * `finite == ±∞` and `±∞ == finite` are `false`
    * `+0.0` and `-0.0` are equal
    """
    if is_nan() or that.is_nan() then
      return false
    end
    if _inf or that._inf then
      // Both must be infinite and share the same sign; one finite → false.
      return _inf and that._inf and (_sign == that._sign)
    end
    eq_unsafe(that)


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
      if _sign then
        return not (that.is_infinite() and that._sign)
      else
        return false
      end
    end
    if that.is_infinite() then
      if that._sign then
        return false
      else
        return not(is_infinite() and not _sign)
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
    if _sign != that._sign then
      return false
    end
    if _exponent != that._exponent then
      return false
    end
    // Compare the common prefix digit by digit.
    let na: USize = _size()
    let nb: USize = that._size()
    try
      var i: USize = 0
      while i < na.min(nb) do
        if _digits(i)? != that._digits(i)? then
          return false
        end
        i = i + 1
      end
      // If one array is longer, any extra nonzero byte breaks equality.
      if na > nb then
        var j: USize = nb
        while j < na do
          if _digits(j)? != 0 then return false end
          j = j + 1
        end
      elseif nb > na then
        var j: USize = na
        while j < nb do
          if that._digits(j)? != 0 then return false end
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
      return not that._sign  // 0 < that iff that > 0 (not negative)
    end
    if that.is_zero() then
      return _sign  // this < 0 iff this is negative
    end
    if _sign then
      if not that._sign then
        // negative < positive → always true
        return true
      end
      // Both negative: this < that ↔ |this| > |that|.
      // Compare exponents first (larger exponent = larger magnitude).
      if _exponent > that._exponent then
        return true
      end
      if _exponent < that._exponent then
        return false
      end

      // Same exponent: compare digits most-significant first.
      let na: USize = _size()
      let nb: USize = that._size()
      try
        var i: USize = 0
        while i < na.min(nb) do
          let ad: U8 = _digits(i)?
          let bd: U8 = that._digits(i)?
          if ad > bd then return true end
          if ad < bd then return false end
          i = i + 1
        end
        // Common prefix equal: extra nonzero digits in this → larger magnitude.
        if na > nb then
          var j: USize = nb
          while j < na do
            if _digits(j)? != 0 then return true end
            j = j + 1
          end
        end
      end
      false  // equal magnitude → not less
    else
      if that._sign then
        // positive < negative → always false
        return false
      end
      // Both positive: compare exponents then digits.
      if _exponent < that._exponent then return true end
      if _exponent > that._exponent then return false end
      let na: USize = _size()
      let nb: USize = that._size()
      try
        var i: USize = 0
        while i < na.min(nb) do
          let ad: U8 = _digits(i)?
          let bd: U8 = that._digits(i)?
          if ad < bd then return true end
          if ad > bd then return false end
          i = i + 1
        end
        // Common prefix equal: extra nonzero digits in that → that is larger.
        if nb > na then
          var j: USize = na
          while j < nb do
            if that._digits(j)? != 0 then return true end
            j = j + 1
          end
        end
      end
      false  // equal → not less
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
  //
  // Private helpers
  // ───────────────
  //   _atanh_series(x, p)  arctanh(x) via Σ x^{2k+1}/(2k+1)
  //   _ln2_const(p)        ln(2) = 2 × arctanh(1/3)
  //   _exp_taylor(r, p)    e^r   via Σ r^k/k!
  //
  // Public functions (all return a new MPFloat at the same precision as this)
  // ───────────────────────────────────────────────────────────────────────────
  //   ln / log             natural logarithm
  //   log2 / log10         logarithm base 2 and 10
  //   logb(base)           logarithm in an arbitrary base
  //   exp                  e^this
  //   exp2                 2^this
  //   powi(n)              this^n  for integer n (exact repeated squaring)
  //   pow(that)            this^that (general real power, positive base)

  fun _atanh_series(x: MPFloat, p: USize): MPFloat =>
    """
    Compute `arctanh(x)` via the Taylor series
    `Σ_{k=0}^∞ x^{2k+1} / (2k+1)`.

    The caller must ensure `|x| < 1` for convergence.  `p` is the working
    precision in bytes; the loop stops when the current term no longer
    changes the running sum at that precision.

    Uses `_short_div` (exact U8 division, no Newton error) for divisors
    up to 255; falls back to `div` for larger denominators (which only arise
    when `p` exceeds ~35 bytes for the 2/3 seed used by `log10`).
    """
    let p2: USize = p + 2
    let x2: MPFloat = x.mul(x)._trunc(p2)         // x^2 — advances the power
    var term: MPFloat = x._trunc(p2)               // running odd power x^{2k+1}
    var sum: MPFloat  = term                       // accumulate here
    var k: USize = 3                               // denominator (1,3,5,…)
    var iters: USize = 0
    let max_iters: USize = p2 * 30                 // safety: >10× what's needed
    while iters < max_iters do
      term = term.mul(x2)._trunc(p2)              // term = x^k
      // Divide term by k to get the series addend x^k / k.
      let addend: MPFloat =
        if k <= 255 then
          (let q, _) = term._short_div(k.u8()) ; q
        else
          term.div(MPFloat.from_f64(k.f64(), (p2 * _base_bits).ulong(), _rounding))._trunc(p2)
        end
      let new_sum: MPFloat = sum.add(addend)._trunc(p2)
      if new_sum.eq(sum) then break end            // converged
      sum = new_sum
      k = k + 2
      iters = iters + 1
    end
    sum._trunc(p)


  fun _ln2_const(p: USize): MPFloat =>
    """
    Compute `ln(2)` to `p` bytes of working precision using the identity
    `ln(2) = 2 × arctanh(1/3)`.

    The rational `1/3` has the exact repeating base-256 representation
    `[0x55, 0x55, …]` (since `85/255 = 1/3`).  A `p+4`-byte approximation
    introduces an error below `256^{-(p+4)}/3`, well within the guard digits.
    """
    let p2: USize = p + 4
    let d: Array[U8] val = recover Array[U8].init(0x55, p2) end
    let third = MPFloat._create(false, false, false, 0, d, _rounding)
    let two = MPFloat.from_f64(2.0, (p2 * _base_bits).ulong(), _rounding)
    _atanh_series(third, p2).mul(two)._trunc(p)


  fun _exp_taylor(r: MPFloat, p: USize): MPFloat =>
    """
    Compute `e^r` via the Taylor series `Σ_{k=0}^∞ r^k / k!`.

    Works best for `|r| ≤ ln(2)/2 ≈ 0.347` where roughly `2.4 × p` terms
    suffice.  Uses the recurrence `term_k = term_{k-1} × r / k` to avoid
    computing factorials.
    """
    let p2: USize = p + 2
    let one = MPFloat.from_f64(1.0, (p2 * _base_bits).ulong(), _rounding)
    var term: MPFloat = one                        // term_0 = 1
    var sum:  MPFloat = one                        // sum starts at 1
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term = term.mul(r)._trunc(p2)               // term_{k-1} × r  (not yet /k)
      // Divide by k: term becomes r^k/k! via the recurrence term_k = term_{k-1} × r / k.
      let divided: MPFloat =
        if k <= 255 then
          (let q, _) = term._short_div(k.u8()) ; q
        else
          term.div(MPFloat.from_f64(k.f64(), (p2 * _base_bits).ulong(), _rounding))._trunc(p2)
        end
      term = divided                               // update: term now holds r^k/k!
      let new_sum: MPFloat = sum.add(divided)._trunc(p2)
      if new_sum.eq(sum) then break end
      sum = new_sum
      k = k + 1
      iters = iters + 1
    end
    sum._trunc(p)


  fun ln(): MPFloat =>
    """
    Natural logarithm `ln(this)`.

    - `NaN`     → `NaN`
    - `+∞`      → `+∞`
    - `x ≤ 0`   → `NaN`  (ln is not real for non-positive arguments)
    - `0`       → `−∞`

    **Algorithm** — two levels of argument reduction, then an arctanh series:

    1. Write `|this| = 0.d₀d₁… × 256^e`.  Set `k = e − 1` (the base-256
       exponent) and `m = 0.d₀d₁… × 256 ∈ [1, 256)`, so
       `ln(this) = k × ln(256) + ln(m)`.

    2. Find `n = ⌊log₂(d₀)⌋ ∈ {0,…,7}` and set `u = m / 2^n ∈ [1, 2)`.
       Then `ln(m) = n × ln(2) + ln(u)`.

    3. Let `t = (u − 1)/(u + 1) ∈ [0, 1/3)`.  Then
       `ln(u) = 2 × arctanh(t)` — the series converges as `(1/9)^k`.

    Combined: `ln(this) = (8k + n) × ln(2) + 2 × arctanh(t)`.
    """
    if _nan then
      return MPFloat.nan_val()
    end
    if _inf and (not _sign) then
      return MPFloat._create(false, false, true, 0, Array[U8].create(), _rounding)
    end
    // ln(negative) = NaN
    if _sign then
      return MPFloat.nan_val()
    end
    if is_zero() then
      return MPFloat._create(true, false, true, 0, Array[U8].create(), _rounding) // -inf
    end

    let p:  USize = _size()
    let p2: USize = p + 6 // guard bytes for cancellation

    let ln2: MPFloat = _ln2_const(p2)

    // Step 1: extract base-256 exponent k and mantissa m ∈ [1, 256).
    // this = 0.d₀d₁… × 256^_exponent, so m = 0.d₀d₁… × 256, k = _exponent − 1.
    let k: I64 = _exponent - 1

    // Build m at working precision p2: zero-pad _digits to p2 bytes.
    let m_digits: Array[U8] val = recover
      let a = Array[U8].init(0, p2)
      var i: USize = 0
      while i < _size().min(p2) do
        try a(i)? = _digits(i)? end
        i = i + 1
      end
      a
    end
    let m = MPFloat._create(false, false, false, 1, m_digits, _rounding)

    // Step 2: reduce m to [1, 2) by dividing by 2^n where n = floor(log2(d₀)).
    let d0: U8 = try _digits(0)? else 1 end
    var n_inner: I64 = 0
    var tmp: U8 = d0
    while tmp >= 2 do
      tmp = tmp >> 1
      n_inner = n_inner + 1
    end
    // u = m / 2^n_inner ∈ [1, 2).  _short_div is an exact integer division.
    let divisor: U8 = U8(1).shl(n_inner.u8())     // 2^n_inner ∈ {1,…,128}
    (let u, _) = m._short_div(divisor)

    // Step 3: ln(u) = 2 × arctanh((u−1)/(u+1)), t ∈ [0, 1/3).
    let one = MPFloat.from_f64(1.0, (p2 * _base_bits).ulong(), _rounding)
    let two = MPFloat.from_f64(2.0, (p2 * _base_bits).ulong(), _rounding)
    let t: MPFloat = (u - one)._trunc(p2) / (u + one)._trunc(p2)
    let ln_u: MPFloat = (_atanh_series(t, p2) * two)._trunc(p2)

    // Combine: ln(this) = (8k + n_inner) × ln(2) + ln(u).
    let ln2_factor: I64 = (8 * k) + n_inner
    let correction: MPFloat = if ln2_factor == 0 then
        MPFloat.create((p2 * _base_bits).ulong(), _rounding)             // zero
      else
        let fac = MPFloat.from_f64(ln2_factor.f64(), (p2 * _base_bits).ulong(), _rounding)
        (fac * ln2)._trunc(p2)
      end
    (correction + ln_u)._trunc(p)


  fun log(): MPFloat =>
    """
    Natural logarithm of `this`.  Alias for [`ln`](#ln).
    """
    ln()


  fun log2(): MPFloat =>
    """
    Base-2 logarithm of `this`: `log₂(this) = ln(this) / ln(2)`.

    Special cases follow those of [`ln`](#ln); additionally `log₂(2) = 1`
    exactly (both are computed to full precision, so the ratio is 1.0 within
    the last digit).
    """
    let p:  USize = _size()
    let p2: USize = p + 4
    let ln_x:  MPFloat = ln()
    if ln_x.is_nan() or ln_x.is_infinite() then
      return ln_x
    end
    (ln_x / _ln2_const(p2))._trunc(p)


  fun log10(): MPFloat =>
    """
    Base-10 logarithm of `this`: `log₁₀(this) = ln(this) / ln(10)`.

    `ln(10) = ln(2) + ln(5) = 2×arctanh(1/3) + 2×arctanh(2/3)`.
    `2/3` has the exact repeating base-256 representation `[0xAA, 0xAA, …]`
    (since `170/255 = 2/3`).
    """
    let p:  USize = _size()
    let p2: USize = p + 4
    let ln_x: MPFloat = ln()
    if ln_x.is_nan() or ln_x.is_infinite() then
      return ln_x
    end
    // ln(10) = 2×arctanh(1/3) + 2×arctanh(2/3)
    let p3: USize = p2 + 2
    let d3:  Array[U8] val = recover Array[U8].init(0x55, p3) end  // 1/3
    let d23: Array[U8] val = recover Array[U8].init(0xAA, p3) end  // 2/3
    let third = MPFloat._create(false, false, false, 0, d3,  _rounding)
    let two_third = MPFloat._create(false, false, false, 0, d23, _rounding)
    let two = MPFloat.from_f64(2.0, (p3 * _base_bits).ulong(), _rounding)
    let ln2 = (_atanh_series(third, p3) * two)._trunc(p3)
    let ln5 = (_atanh_series(two_third, p3) * two)._trunc(p3)
    let ln10: MPFloat = (ln2 + ln5)._trunc(p2)
    (ln_x / ln10)._trunc(p)


  fun logb(base: MPFloat): MPFloat =>
    """
    Logarithm of `this` in the given `base`: `logb(this, b) = ln(this) / ln(b)`.

    `base` must be positive and not equal to 1; otherwise the result is NaN
    or ±∞ following the same conventions as [`ln`](#ln) and `div`.
    """
    let p = _size()
    let p2 = p + 4
    (ln() / base.ln())._trunc(p)


  fun exp(): MPFloat =>
    """
    Natural exponential `e^this`.

    - `NaN`  → `NaN`
    - `+∞`   → `+∞`
    - `−∞`   → `+0`
    - `0`    → `1`

    **Algorithm** — two-level argument reduction, then a Taylor series:

    1. Write `this = n₂₅₆ × ln(256) + r₁` where `n₂₅₆ = round(this / ln(256))`
       and `|r₁| ≤ ln(256)/2 ≈ 2.77`.  Then `e^this = 256^n₂₅₆ × e^r₁`,
       and scaling by `256^n₂₅₆` is exact (adjust `_exponent`).

    2. Write `r₁ = n₂ × ln(2) + r` where `n₂ = round(r₁ / ln(2)) ∈ {-4,…,4}`
       and `|r| ≤ ln(2)/2 ≈ 0.347`.  Then `e^r₁ = 2^n₂ × e^r`, and
       `2^n₂` is an exactly-representable 1-byte MPFloat.

    3. Compute `e^r` via Taylor series (≈ `2.4 × p` terms for `p` bytes).
    """
    if _nan then
      return MPFloat.nan_val()
    end
    if _inf and (not _sign) then
      return MPFloat._create(false, false, true, 0, Array[U8].create(), _rounding)
    end
    if _inf then // exp(-inf) = +0
      return MPFloat.create(get_precision(), _rounding)
    end
    if is_zero() then
      return MPFloat.from_f64(1.0, get_precision(), _rounding)
    end

    let p:  USize = _size()
    let p2: USize = p + 8 // extra guard for cancellation

    let ln2:  MPFloat = _ln2_const(p2)
    // ln(256) = 8 × ln(2).  Multiply by 8 using a 1-byte MPFloat for 8:
    // [8] with _exponent=1 represents 8/256 × 256^1 = 8.
    let eight = MPFloat._create(false, false, false, 1, recover [U8(8)] end, _rounding)
    let ln256: MPFloat = (ln2 * eight)._trunc(p2)

    // x as a full-precision MPFloat at p2 bytes (zero-padded to p2 digits).
    let x_digits: Array[U8] val = recover
      let a = Array[U8].init(0, p2)
      var i: USize = 0
      while i < _size().min(p2) do
        try a(i)? = _digits(i)? end
        i = i + 1
      end
      a
    end
    let x = MPFloat._create(_sign, false, false, _exponent, x_digits, _rounding)

    // Step 1: n₂₅₆ = round(x / ln256); r₁ = x − n₂₅₆ × ln(256).
    let n256_f: MPFloat = (x / ln256)._trunc(p2).round()
    let n256: ILong = (try MPInt.from_mpfloat(n256_f)? else MPInt.from[ILong](0) end).ilong()
    let n256_mpf = MPFloat.from_f64(n256.f64(), (p2 * _base_bits).ulong(), _rounding)
    let r1: MPFloat = x.sub(n256_mpf.mul(ln256)._trunc(p2))._trunc(p2)

    // Step 2: n₂ = round(r₁ / ln2) ∈ {−4,…,4}; r = r₁ − n₂ × ln(2).
    let n2_f: MPFloat = (r1 / ln2)._trunc(p2).round()
    let n2: ILong = (try MPInt.from_mpfloat(n2_f)? else MPInt.from[ILong](0) end).ilong()
    let n2_mpf = MPFloat.from_f64(n2.f64(), (p2 * _base_bits).ulong(), _rounding)
    let r: MPFloat = (r1 - n2_mpf.mul(ln2)._trunc(p2))._trunc(p2)

    // Step 3: Taylor series for e^r, |r| ≤ ln(2)/2.
    var result: MPFloat = _exp_taylor(r, p2)

    // Scale by 2^n₂: construct an exact 1-byte MPFloat for 2^|n₂|.
    if n2 > 0 then
      // 2^n₂ = [2^n₂] with _exponent=1  (digit ∈ {2,4,8,16}, fits in U8).
      let pow2 = MPFloat._create(false, false, false, 1, recover [U8(1).shl(n2.u8())] end, _rounding)
      result = (result * pow2)._trunc(p2)
    elseif n2 < 0 then
      // 2^n₂ = [128 >> (|n₂|-1)] with _exponent=0  (digit ∈ {64,32,16,8}).
      let pow2 = MPFloat._create(false, false, false, 0, recover [U8(128).shr((-n2 - 1).u8())] end, _rounding)
      result = (result * pow2)._trunc(p2)
    end

    // Scale by 256^n₂₅₆: just adjust the exponent field.
    result = MPFloat._create(result._sign, false, false, result._exponent + n256.i64(), result._digits, _rounding)

    result._trunc(p)


  fun exp2(): MPFloat =>
    """
    Base-2 exponential `2^this = exp(this × ln(2))`.

    - `NaN`  → `NaN`.  `+∞` → `+∞`.  `−∞` → `+0`.  `0` → `1`.

    For integer `this`, the result is exact to full precision (the argument
    reduction in `exp` will produce `r = 0` exactly).
    """
    if _nan then
      return MPFloat.nan_val()
    end
    if _inf and (not _sign) then
      return MPFloat._create(false, false, true, 0, Array[U8].create(), _rounding)
    end
    if _inf then
      return MPFloat.create(get_precision(), _rounding)
    end
    if is_zero() then
      return MPFloat.from_f64(1.0, get_precision(), _rounding)
    end
    let p2: USize = _size() + 4
    let ln2: MPFloat = _ln2_const(p2)
    (this * ln2)._trunc(p2).exp()._trunc(_size())


  fun powi(n: ILong, rnd: RoundingMode = RoundingNearest): MPFloat =>
    """
    Raise `this` to the integer power `n` using binary (square-and-multiply)
    exponentiation.  Exact for all finite inputs (no Newton approximation).

    - `x^0 = 1` for any finite `x` (including `0^0 = 1` by convention).
    - `0^n = 0` for `n > 0`.
    - `x^n` for `n < 0` is computed as `(1/x)^|n|` via `inv()`.
    - NaN → NaN.  ±∞ → ±∞ (for `n > 0`) or `0` (for `n < 0`).
    """
    if _nan then
      return MPFloat.nan_val()
    end
    
    let p: USize = _size()
    if n == 0 then
      return MPFloat.from_f64(1.0, (p * _base_bits).ulong(), _rounding)
    end
    
    // For negative exponent: compute inv() first, then positive power.
    var base = if n < 0 then inv() else MPFloat._create(_sign, _nan, _inf, _exponent, _digits, _rounding) end
    var exp_n: ILong = if n < 0 then -n else n end
    var result = MPFloat.from_f64(1.0, (p * _base_bits).ulong(), _rounding)
    while exp_n > 0 do
      if (exp_n and 1) == 1 then
        result = (result * base)._trunc(p)
      end
      base = (base * base)._trunc(p)
      exp_n = exp_n / 2
    end
    result


  fun pow(that: MPFloat): MPFloat =>
    """
    General real power `this^that = exp(that × ln(this))`.

    - `this > 0`: computed via `exp(that × ln(this))`.
    - `this = 1`: returns `1` regardless of `that`.
    - `this = 0` and `that > 0`: returns `0`.
    - `this = 0` and `that ≤ 0`: returns `NaN`.
    - `this < 0` and `that` is an integer: uses `powi`.
    - `this < 0` and `that` is not an integer: returns `NaN`.
    - NaN in either operand → `NaN`.
    """
    let p: USize = _size()
    let p2: USize = p + 4

    if _nan or that._nan then
      return MPFloat.nan_val()
    end
    // x^0 = 1
    if that.is_zero() then
      return MPFloat.from_f64(1.0, (p * _base_bits).ulong(), _rounding)
    end
    if is_zero() then
      if that.is_negative() then
        return MPFloat.nan_val()
      end
      return MPFloat.create((p * _base_bits).ulong(), _rounding)                                  // 0^pos = 0
    end

    // Positive base: exact path.
    if not _sign then
      let ln_x: MPFloat = ln()
      if ln_x.is_nan() then
        return MPFloat.nan_val()
      end
      return (that * ln_x)._trunc(p2).exp()._trunc(p)
    end

    // Negative base: only valid for integer exponents.
    if that.is_integer() then
      let ni: ILong = (try MPInt.from_mpfloat(that)? else MPInt.from[ILong](0) end).ilong()
      return powi(ni)
    end
    MPFloat.nan_val()


  //- Trigonometric and hyperbolic functions -----------------------------------

  fun _sin_taylor(r: MPFloat, p: USize): MPFloat =>
    """
    Compute `sin(r)` via the Taylor series `r − r³/3! + r⁵/5! − r⁷/7! + …`.

    The recurrence `term_{k+1} = term_k × (−r²) / ((2k)(2k+1))` advances two
    factorial positions at once.  Converges for any finite `r`; works best for
    `|r| ≤ π/4 ≈ 0.785` where the argument reduction in `sin()` leaves `r`.
    """
    let p2: USize = p + 2
    let neg_r2: MPFloat = (r * r)._trunc(p2).neg()  // −r²
    var term: MPFloat = r._trunc(p2)                 // r¹/1!
    var sum:  MPFloat = term
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30 // TODO Check constant
    while iters < max_iters do
      term = (term * neg_r2)._trunc(p2)
      let d1: USize = 2 * k
      let d2: USize = (2 * k) + 1
      term = if d1 <= 255 then
          (let q, _) = term._short_div(d1.u8())
          q
        else
          (term / MPFloat.from_f64(d1.f64(), (p2 * _base_bits).ulong(), _rounding))._trunc(p2)
        end
      term = if d2 <= 255 then
          (let q, _) = term._short_div(d2.u8())
          q
        else
          term.div(MPFloat.from_f64(d2.f64(), (p2 * _base_bits).ulong(), _rounding))._trunc(p2)
        end
      let new_sum = (sum + term)._trunc(p2)

      if new_sum == sum then
        break
      end
      sum = new_sum
      k = k + 1
      iters = iters + 1
    end
    sum._trunc(p)


  fun _cos_taylor(r: MPFloat, p: USize): MPFloat =>
    """
    Compute `cos(r)` via the Taylor series `1 − r²/2! + r⁴/4! − r⁶/6! + …`.

    The recurrence `term_{k+1} = term_k × (−r²) / ((2k−1)(2k))` advances two
    factorial positions at once.  Converges for any finite `r`; works best for
    `|r| ≤ π/4 ≈ 0.785`.
    """
    let p2: USize = p + 2
    let neg_r2 = (r * r)._trunc(p2).neg() // −r²
    let one = MPFloat.from_f64(1.0, (p2 * _base_bits).ulong(), _rounding)
    var term = one // r⁰/0! = 1
    var sum = term
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30 // TODO Check constant
    while iters < max_iters do
      term = (term * neg_r2)._trunc(p2)
      let d1: USize = (2 * k) - 1
      let d2: USize = 2 * k
      term = if d1 <= 255 then
          (let q, _) = term._short_div(d1.u8())
          q
        else
          (term / MPFloat.from_f64(d1.f64(), (p2 * _base_bits).ulong(), _rounding))._trunc(p2)
        end
      term = if d2 <= 255 then
          (let q, _) = term._short_div(d2.u8())
          q
        else
          (term / MPFloat.from_f64(d2.f64(), (p2 * _base_bits).ulong(), _rounding))._trunc(p2)
        end
      let new_sum = (sum + term)._trunc(p2)

      if new_sum == sum then
        break
      end
      sum = new_sum
      k = k + 1
      iters = iters + 1
    end
    sum._trunc(p)


  fun sin(): MPFloat =>
    """
    Sine of `this` (argument in radians).

    - `NaN` or `±∞` → `NaN`  (sine is not defined at infinity)
    - `0` → `0`

    **Algorithm** — argument reduction then Taylor series:

    1. Compute `n = round(this / (π/2))` and `r = this − n × π/2`,
       so `|r| ≤ π/4`.
    2. Based on `k = n mod 4`:
       - `k=0`: `sin = sin(r)`
       - `k=1`: `sin = cos(r)`
       - `k=2`: `sin = −sin(r)`
       - `k=3`: `sin = −cos(r)`

    **Precision note**: argument reduction subtracts two nearly-equal numbers
    when `|this|` is much larger than `π/2`.  For reliable results keep
    `|this| ≲ 256^(_size() − 1)`.
    """
    if _nan or _inf then
      return MPFloat.nan_val()
    end
    if is_zero() then
      return MPFloat.create(get_precision(), _rounding)
    end

    let p: USize = _size()
    let p2: USize = p + 8 // extra guard for argument-reduction cancellation

    // Extend x to p2 bytes.
    let x_digits: Array[U8] val = recover
      let a = Array[U8].init(0, p2)
      var i: USize = 0
      /* TODO DELETE
      while i < _size().min(p2) do
        try a(i)? = _digits(i)? end
        i = i + 1
      end
      */
      _digits.copy_to(a, 0, 0, _size().min(p2))
      a
    end
    let x = MPFloat._create(_sign, false, false, _exponent, x_digits, _rounding)

    let pi_val = MPFloat.pi((p2 * _base_bits).ulong(), _rounding)
    let two = MPFloat.from_f64(2.0, (p2 * _base_bits).ulong(), _rounding)
    let pi_half = (pi_val / two)._trunc(p2)

    // n = round(x / (π/2)), r = x − n × (π/2).
    let xdivph = (x / pi_half)._trunc(p2)
    let n_f = xdivph.round()
    let n: ILong = (try MPInt.from_mpfloat(n_f)? else MPInt.from[ILong](0) end).ilong()
    let n_mpf = MPFloat.from_mpint(MPInt.from[ILong](n), (p2 * _base_bits).ulong())
    let r = (x - n_mpf.mul(pi_half)._trunc(p2))._trunc(p2)

    // Octant k = n mod 4, normalised to [0, 3].
    let k: ILong = ((n % 4) + 4) % 4

    // Only one Taylor series is needed per octant.
    let result: MPFloat = if (k == 0) or (k == 2) then
        let sr: MPFloat = _sin_taylor(r, p2)
        if k == 0 then
          sr
        else
          sr.neg()
        end
      else
        let cr: MPFloat = _cos_taylor(r, p2)
        if k == 1 then
          cr
        else
          cr.neg()
        end
      end
    result._trunc(p)


  fun cos(): MPFloat =>
    """
    Cosine of `this` (argument in radians).

    - `NaN` or `±∞` → `NaN`
    - `0` → `1`

    Uses the same argument reduction as `sin()`:
    - `k=0`: `cos = cos(r)`
    - `k=1`: `cos = −sin(r)`
    - `k=2`: `cos = −cos(r)`
    - `k=3`: `cos = sin(r)`
    """
    if _nan or _inf then
      return MPFloat.nan_val()
    end
    if is_zero() then
      return MPFloat.from_f64(1.0, get_precision(), _rounding)
    end

    let p: USize = _size()
    let p2: USize = p + 8

    let x_digits: Array[U8] val = recover
      let a = Array[U8].init(0, p2)
      /* TODO DELETE
      var i: USize = 0
      while i < _size().min(p2) do try a(i)? = _digits(i)? end ; i = i + 1 end
      */
      _digits.copy_to(a, 0, 0, _size().min(p2))
      a
    end
    let x = MPFloat._create(_sign, false, false, _exponent, x_digits, _rounding)

    let pi_val = MPFloat.pi((p2 * _base_bits).ulong(), _rounding)
    let two = MPFloat.from_f64(2.0, (p2 * _base_bits).ulong(), _rounding)
    let pi_half = (pi_val / two)._trunc(p2)

    let n_f = (x / pi_half)._trunc(p2).round()
    let n: ILong = (try MPInt.from_mpfloat(n_f)? else MPInt.from[ILong](0) end).ilong()
    let n_mpf = MPFloat.from_mpint(MPInt.from[ILong](n), (p2 * _base_bits).ulong())
    let r = (x - n_mpf.mul(pi_half)._trunc(p2))._trunc(p2)

    let k: ILong = ((n % 4) + 4) % 4

    let result: MPFloat = if (k == 0) or (k == 2) then
        let cr = _cos_taylor(r, p2)
        if k == 0 then
          cr
        else
          cr.neg()
        end
      else
        let sr = _sin_taylor(r, p2)
        if k == 3 then
          sr
        else
          sr.neg()
        end
      end
    result._trunc(p)


  fun tan(): MPFloat =>
    """
    Tangent of `this` (argument in radians): `sin(this) / cos(this)`.

    - `NaN` or `±∞` → `NaN`
    - `0` → `0`

    When `cos(this)` underflows to zero at finite precision (near odd multiples
    of `π/2`), `div` returns `±∞`.
    """
    if _nan or _inf then
      return MPFloat.nan_val()
    end
    if is_zero() then
      return MPFloat.create(get_precision(), _rounding)
    end
    let c: MPFloat = cos()
    if c.is_zero() then
      return MPFloat.nan_val()
    end
    (sin() / c)._trunc(_size())


  fun csc(): MPFloat =>
    """
    Cosecant of `this`: `1 / sin(this)`.

    - `NaN`, `±∞`, or `sin(this) = 0` → `NaN`
    """
    if _nan or _inf then
      return MPFloat.nan_val()
    end
    let s: MPFloat = sin()
    if s.is_zero() then
      return MPFloat.nan_val()
    end
    s.inv()._trunc(_size())


  fun sec(): MPFloat =>
    """
    Secant of `this`: `1 / cos(this)`.

    - `NaN`, `±∞`, or `cos(this) = 0` → `NaN`
    """
    if _nan or _inf then
      return MPFloat.nan_val()
    end
    let c: MPFloat = cos()
    if c.is_zero() then
      return MPFloat.nan_val()
    end
    c.inv()._trunc(_size())


  fun cot(): MPFloat =>
    """
    Cotangent of `this`: `cos(this) / sin(this)`.

    - `NaN`, `±∞`, or `sin(this) = 0` → `NaN`
    """
    if _nan or _inf then
      return MPFloat.nan_val()
    end
    let s: MPFloat = sin()
    if s.is_zero() then
      return MPFloat.nan_val()
    end
    (cos() / s)._trunc(_size())


  fun sinh(): MPFloat =>
    """
    Hyperbolic sine: `(e^this − e^{−this}) / 2`.

    - `NaN` → `NaN`
    - `±∞` → `±∞`
    - `0` → `0`

    **Precision note**: for very small `|this|` the subtraction
    `e^x − e^{−x}` cancels leading digits.  Use higher precision `p` when
    `|this| ≪ 256^{−(p−1)/2}`.
    """
    if _nan then
      return MPFloat.nan_val()
    end
    if _inf  then
      return MPFloat._create(_sign, false, true, 0, Array[U8].create(), _rounding)
    end
    if is_zero() then
      return MPFloat.create(get_precision(), _rounding)
    end

    let p: USize = _size()
    let p2: USize = p + 4
    let ex:  MPFloat = exp()
    let emx: MPFloat = neg().exp()
    let two = MPFloat.from_f64(2.0, (p2 * _base_bits).ulong(), _rounding)
    ((ex - emx)._trunc(p2) / two)._trunc(p)


  fun cosh(): MPFloat =>
    """
    Hyperbolic cosine: `(e^this + e^{−this}) / 2`.

    - `NaN` → `NaN`
    - `±∞` → `+∞`  (cosh is always non-negative)
    - `0` → `1`
    """
    if _nan then
      return MPFloat.nan_val()
    end
    if _inf  then
      return MPFloat._create(false, false, true, 0, Array[U8].create(), _rounding)
    end
    if is_zero() then
      return MPFloat.from_f64(1.0, get_precision(), _rounding)
    end
    let p: USize = _size()
    let p2: USize = p + 4
    let ex:  MPFloat = exp()
    let emx: MPFloat = neg().exp()
    let two = MPFloat.from_f64(2.0, (p2 * _base_bits).ulong(), _rounding)
    ((ex + emx)._trunc(p2) / two)._trunc(p)


  fun tanh(): MPFloat =>
    """
    Hyperbolic tangent: `(e^this − e^{−this}) / (e^this + e^{−this})`.

    - `NaN` → `NaN`
    - `+∞` → `+1`,  `−∞` → `−1`
    - `0` → `0`
    """
    if _nan then
      return MPFloat.nan_val()
    end
    if _inf  then
      return MPFloat.from_f64(if _sign then -1.0 else 1.0 end, get_precision(), _rounding)
    end
    if is_zero() then
      return MPFloat.create(get_precision(), _rounding)
    end
    let p: USize = _size()
    let p2: USize = p + 4
    let ex:  MPFloat = exp()
    let emx: MPFloat = neg().exp()
    ((ex - emx)._trunc(p2) / (ex + emx))._trunc(p)


  fun csch(): MPFloat =>
    """
    Hyperbolic cosecant: `1 / sinh(this)`.

    - `NaN`, `±∞`, or `sinh(this) = 0` (i.e. `this = 0`) → `NaN`
    """
    if _nan or _inf then
      return MPFloat.nan_val()
    end
    if is_zero() then
      return MPFloat.nan_val()
    end
    sinh().inv()._trunc(_size())


  fun sech(): MPFloat =>
    """
    Hyperbolic secant: `1 / cosh(this)`.

    - `NaN` → `NaN`
    - `±∞` → `0`  (cosh(±∞) = +∞)
    """
    if _nan then
      return MPFloat.nan_val()
    end
    if _inf  then
      return MPFloat.create(get_precision(), _rounding)
    end
    cosh().inv()._trunc(_size())


  fun coth(): MPFloat =>
    """
    Hyperbolic cotangent: `cosh(this) / sinh(this)`.

    - `NaN`, `±∞`, or `this = 0` → `NaN`
    """
    if _nan or _inf then
      return MPFloat.nan_val()
    end
    if is_zero() then
      return MPFloat.nan_val()
    end
    (cosh() / sinh())._trunc(_size())


  //- Rounding and truncation -------------------------------------------------

  fun abs(): MPFloat =>
    """
    Absolute value of `this`, when it makes sense. Result is undefined for NaN.
    """
    MPFloat._create(false, _nan, _inf, _exponent, _digits, _rounding)


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
    if _nan or _inf then
      return MPFloat._create(_sign, _nan, _inf, _exponent, _digits, _rounding)
    end
    if is_zero() or (_exponent <= 0) then
      return MPFloat.create(get_precision(), _rounding)
    end
    let e: USize = _exponent.usize()
    if e >= _size() then
      return MPFloat._create(_sign, false, false, _exponent, _digits, _rounding)
    end
    let new_digits: Array[U8] val = recover
      let d = Array[U8].create(e)
      /* TODO DELETE
      var i: USize = 0
      while i < e do
        try d.push(_digits(i)?) end
        i = i + 1
      end
      */
      _digits.copy_to(d, 0, 0, e)
      d
    end
    MPFloat._create(_sign, false, false, _exponent, new_digits, _rounding)


  fun floor(): MPFloat =>
    """
    Floor: the largest integer less than or equal to `this`.

    - NaN → NaN.  ±∞ → ±∞.  ±0 → ±0.
    - Non-negative or exact integer → `trunc(this)`.
    - Negative with a non-zero fractional part → `trunc(this) − 1`.

    See also `trunc`, `ceil`, `round`.
    """
    let t: MPFloat = trunc()
    if _sign and _has_frac() then
      return (t - MPFloat.from_f64(1.0, get_precision(), _rounding))
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
    let t: MPFloat = trunc()
    if (not _sign) and _has_frac() then
      return (t + MPFloat.from_f64(1.0, get_precision(), _rounding))
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
    let half = MPFloat.from_f64(0.5, get_precision(), _rounding)
    if _sign then
      return (this  - half).ceil()
    end
    (this + half).floor()


  //- Hashing -----------------------------------------------------------------

  fun hash(): USize =>
    """
    Calculate a hash of this `MPFloat`.
    """
    0


  fun hash64(): U64 =>
    """
    Calculate a 64-bits hash of the `MPFloat`.
    """
    0

  
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
// Subnormal numbers (e.g., 2.22507e-308 for F64 and 1.17549e-38 for F32) can cause
// intermediate underflow during the reconstruction of the floating-point value in
// the f64() and f32() methods.
//
// Reasons for Failure
//
// 1. Underflow of the Scaling Factor: In f64(), the code calculates
// `result * (base_f64.pow(_exponent.f64() - size.f64()))`. For very small numbers,
// the exponent `(_exponent - size)` can be -130 or lower. Since 256^-135 = 2^-1080,
// and the minimum representable `F64` (including subnormals) is 2^-1074, the `pow()`
// call returns `0.0`. This makes the entire result 0.0, even if the final value
// should have been a non-zero subnormal.
// 2. Order of Operations: By calculating the scaling factor in isolation, it hits the
// floating-point limits prematurely. Multiplying the mantissa by the factor in stages
// would keep the intermediate results within the representable range.
// 3. Potential Overflow: Conversely, for high-precision `MPFloat` instances, the integer
// result (accumulated mantissa) can overflow to infinity before the scaling factor is
// applied if too many digits are summed.
//
// How to Solve It
//
// 1. Limit Mantissa Accumulation: Only sum enough digits to fill the destination's precision
// (e.g., 16 base-256 digits for `F64` provide 128 bits, which is more than enough for its
// 53-bit mantissa). This prevents intermediate overflow.
// 2. Conditional Scaling Split: Split the `pow()` calculation only when the exponent is
// outside the "safe" range (where `256^e` is guaranteed not to underflow/overflow). For
// `F64`, a safe range is approximately [-125, 125].
// 3. Preserve Accuracy: Use integer-valued exponents for `pow()` to ensure that for normal
// ranges, the calculation remains bit-for-bit identical to the original implementation,
// avoiding 1-ulp regressions.

  fun f64(): F64 =>
    """
    Convert the current `MPFloat` to `F64`. Overflows are converted to ±∞.
    Underflow are converted to 0.

    This always hold:
    * If `f` is a `F64`: `MPFloat.from_f64(f).f64() == f`
    * If `mpf` is a `MPFloat`: `MPFloat.from_f64(mpf.f64()) == mpf`

    The algorithm limits mantissa accumulation to the necessary precision (16 digits) to
    prevent intermediate overflow. It splits the exponent scaling (`base.pow(e)`) into two
    parts when the exponent is large, avoiding intermediate underflow for near-subnormal
    and subnormal values.
    """
    if _nan then
      return F64.from_bits(0x7FF8000000000000)
    end
    if _inf then
      if _sign then
        return F64.from_bits(0xFFF0000000000000)
      else
        return F64.from_bits(0x7FF0000000000000)
      end
    end
    if is_zero() then
      if _sign then
        return F64.from_bits(0x8000000000000000)
      else
        return 0.0
      end
    end

    // Reconstruct the F64 value
    var result: F64 = 0.0
    let size = _size()
    let p = size.min(16) // 16 digits = 128 bits, covers F64 precision
    let base_f64 = _base.f64()
    let inv_base = 1.0 / base_f64

    // Horner's method from least to most significant digit
    try
      var i: USize = p
      while i > 0 do
        i = i - 1
        result = (result + _digits(i)?.f64()) * inv_base
      end
    end

    // Adjust by the exponent avoiding intermediate underflow/overflow
    let final_value = if (_exponent >= -125) and (_exponent <= 125) then
        result * base_f64.pow(_exponent.f64())
      else
        let e1 = _exponent / 2
        let e2 = _exponent - e1
        (result * base_f64.pow(e1.f64())) * base_f64.pow(e2.f64())
      end

    if _sign then
      return -final_value
    else
      return final_value
    end

  
  fun f32(): F32 =>
    """
    Convert the current `MPFloat` to `F32`. Overflows are converted to ±∞.
    Underflow are converted to 0.

    This always hold:
    * If `f` is a `F32`: `MPFloat.from_f32(f).f32() == f`
    * If `mpf` is a `MPFloat`: `MPFloat.from_f32(mpf.f32()) == mpf`

    The algorithm limits mantissa accumulation to the necessary precision (8 digits) to
    prevent intermediate overflow. It splits the exponent scaling (`base.pow(e)`) into two
    parts when the exponent is large, avoiding intermediate underflow for near-subnormal
    and subnormal values.
    """
    if _nan then
      return F32.from_bits(0x7FC00000)
    end
    if _inf then
      if _sign then
        return F32.from_bits(0xFF800000)
      else
        return F32.from_bits(0x7F800000)
      end
    end
    if is_zero() then
      if _sign then
        return F32.from_bits(0x80000000)
      else
        return 0.0
      end
    end

    // Reconstruct the F32 value
    var result: F64 = 0.0
    let size = _size()
    let p = size.min(16)
    let base_f64 = _base.f64()
    let inv_base = 1.0 / base_f64

    // Horner's method from least to most significant digit
    try
      var i: USize = p
      while i > 0 do
        i = i - 1
        result = (result + _digits(i)?.f64()) * inv_base
      end
    end

    // Adjust by the exponent using F64 to maintain precision
    let final_value_f64 = if (_exponent >= -125) and (_exponent <= 125) then
        (result * base_f64.pow(_exponent.f64()))
      else
        let e1 = _exponent / 2
        let e2 = _exponent - e1
        (result * base_f64.pow(e1.f64())) * base_f64.pow(e2.f64())
      end
    
    let final_value = final_value_f64.f32()

    if _sign then
      return -final_value
    else
      return final_value
    end
  

  fun val i128(): I128 =>
    """
    Convert `this` to an `I128` value. If `this` has decimals, they
    are truncated. In case of overflow, the value is saturated to
    `I128.max_value()` or `I128.min_value()`.

    * NaN ⟶ 0
    * -∞ ⟶ `I128.min_value()`
    * +∞ ⟶ `I128.max_value()`
    * ]-1, +1[ ⟶ 0
    * Other values are saturated to the `I128` range.
    """
    if _nan then
      return 0
    end
    if _inf then
      if _sign then
        return I128.min_value()
      else
        return I128.max_value()
      end
    end
    if is_zero() or (_exponent <= 0) then
      return 0
    end

    // Check for overflow before conversion to avoid huge allocations
    // and to implement saturation.
    // I128 range is [-2^127, 2^127 - 1].
    // 2^127 = 128 * 256^15.
    // _exponent is the number of base-256 integer digits.
    if (_exponent > 16) or ((_exponent == 16) and (try _digits(0)? >= 128 else false end)) then
      if _sign then
        return I128.min_value()
      else
        return I128.max_value()
      end
    end

    // It fits in I128 (magnitude < 2^127, or negative and magnitude == 2^127)
    var res: U128 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << _base_bits.u128()) or _digits(i)?.u128()
      end
    end

    if _exponent.usize() > n then
      res = res << (_base_bits.usize() * (_exponent.usize() - n)).u128()
    end

    if _sign then
      (-res).i128()
    else
      res.i128()
    end
  
  
  fun i64(): I64 =>
      """
      Convert `this` to an `I64` value. If `this` has decimals, they
      are truncated. In case of overflow, the value is saturated to
      `I64.max_value()` or `I64.min_value()`.
  
      * NaN ⟶ 0
      * -∞ ⟶ `I64.min_value()`
      * +∞ ⟶ `I64.max_value()`
      * ]-1, +1[ ⟶ 0
      * Other values are saturated to the `I64` range.
      """
      if _nan then
        return 0
      end
      if _inf then
        if _sign then
          return I64.min_value()
        else
          return I64.max_value()
        end
      end
      if is_zero() or (_exponent <= 0) then
        return 0
      end
  
      // Check for overflow before conversion to MPInt to avoid huge allocations
      // and to implement saturation.
      // I64 range is [-2^63, 2^63 - 1].
      // 2^63 = 128 * 256^7.
      // _exponent is the number of base-256 integer digits.
      if (_exponent > 8) or ((_exponent == 8) and (try _digits(0)? >= 128 else false end)) then
        if _sign then
          return I64.min_value()
        else
          return I64.max_value()
        end
      end
  
      // It fits in I64 (magnitude < 2^63, or negative and magnitude == 2^63)
    var res: U64 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << _base_bits.u64()) or _digits(i)?.u64()
      end
    end

    if _exponent.usize() > n then
      res = res << (_base_bits.usize() * (_exponent.usize() - n)).u64()
    end

    if _sign then
      (-res).i64()
    else
      res.i64()
    end
    
    
  fun i32(): I32 =>
    """
    Convert `this` to an `I32` value. If `this` has decimals, they
    are truncated. In case of overflow, the value is saturated to
    `I32.max_value()` or `I32.min_value()`.

    * NaN ⟶ 0
    * -∞ ⟶ `I32.min_value()`
    * +∞ ⟶ `I32.max_value()`
    * ]-1, +1[ ⟶ 0
    * Other values are saturated to the `I32` range.
    """
    if _nan then
      return 0
    end
    if _inf then
      if _sign then
        return I32.min_value()
      else
        return I32.max_value()
      end
    end
    if is_zero() or (_exponent <= 0) then
      return 0
    end

    // Check for overflow before conversion to MPInt to avoid huge allocations
    // and to implement saturation.
    // I32 range is [-2^31, 2^31 - 1].
    // 2^31 = 128 * 256^3.
    // _exponent is the number of base-256 integer digits.
    if (_exponent > 4) or ((_exponent == 4) and (try _digits(0)? >= 128 else false end)) then
      if _sign then
        return I32.min_value()
      else
        return I32.max_value()
      end
    end

    // It fits in I32 (magnitude < 2^31, or negative and magnitude == 2^31)
    var res: U32 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << _base_bits.u32()) or _digits(i)?.u32()
      end
    end

    if _exponent.usize() > n then
      res = res << (_base_bits.usize() * (_exponent.usize() - n)).u32()
    end

    if _sign then
      (-res).i32()
    else
      res.i32()
    end


  fun i16(): I16 =>
    """
    Convert `this` to an `I16` value. If `this` has decimals, they
    are truncated. In case of overflow, the value is saturated to
    `I16.max_value()` or `I16.min_value()`.

    * NaN ⟶ 0
    * -∞ ⟶ `I16.min_value()`
    * +∞ ⟶ `I16.max_value()`
    * ]-1, +1[ ⟶ 0
    * Other values are saturated to the `I16` range.
    """
    if _nan then
      return 0
    end
    if _inf then
      if _sign then
        return I16.min_value()
      else
        return I16.max_value()
      end
    end
    if is_zero() or (_exponent <= 0) then
      return 0
    end

    // Check for overflow before conversion to MPInt to avoid huge allocations
    // and to implement saturation.
    // I16 range is [-2^15, 2^15 - 1].
    // 2^15 = 128 * 256^1.
    // _exponent is the number of base-256 integer digits.
    if (_exponent > 2) or ((_exponent == 2) and (try _digits(0)? >= 128 else false end)) then
      if _sign then
        return I16.min_value()
      else
        return I16.max_value()
      end
    end

    // It fits in I16 (magnitude < 2^15, or negative and magnitude == 2^15)
    //TODO: Optimize by unrolling the loop
    var res: U16 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << _base_bits.u16()) or _digits(i)?.u16()
      end
    end

    if _exponent.usize() > n then
      res = res << (_base_bits.usize() * (_exponent.usize() - n)).u16()
    end

    if _sign then
      (-res).i16()
    else
      res.i16()
    end
    
    
  fun i8(): I8 =>
    """
    Convert `this` to an `I8` value. If `this` has decimals, they
    are truncated. In case of overflow, the value is saturated to
    `I8.max_value()` or `I8.min_value()`.

    * NaN ⟶ 0
    * -∞ ⟶ `I8.min_value()`
    * +∞ ⟶ `I8.max_value()`
    * ]-1, +1[ ⟶ 0
    * Other values are saturated to the `I8` range.
    """
    if _nan then
      return 0
    end
    if _inf then
      if _sign then
        return I8.min_value()
      else
        return I8.max_value()
      end
    end
    if is_zero() or (_exponent <= 0) then
      return 0
    end

    // Check for overflow before conversion to MPInt to avoid huge allocations
    // and to implement saturation.
    // I8 range is [-2^7, 2^7 - 1].
    // 2^7 = 128.
    // _exponent is the number of base-256 integer digits.
    if (_exponent > 1) or ((_exponent == 1) and (try _digits(0)? >= 128 else false end)) then
      if _sign then
        return I8.min_value()
      else
        return I8.max_value()
      end
    end

    // It fits in I8 (magnitude < 2^7, or negative and magnitude == 2^7)
    var res = try _digits(0)?.u8() else 0 end

    if _exponent > 1 then
      res = res << (_base_bits.u8() * (_exponent.u8() - 1))
    end

    if _sign then
      (-res).i8()
    else
      res.i8()
    end
    
    
  fun ilong(): ILong =>
    """
    Convert `this` to a `ILong` value. On LLP64 and ILP32, that's a conversion
    to `I32`, and to `I64` for all other platforms.

    * NaN ⟶ 0
    * -∞ ⟶ `ILong.min_value()`
    * +∞ ⟶ `ILong.max_value()`
    * ]-1, +1[ ⟶ 0
    * Other values are saturated to the `ILong` range.
    """
    ifdef lp64 then
      i64().ilong()
    elseif llp64 then
      i32().ilong()
    elseif ilp32 then
      i32().ilong()
    else
      i64().ilong()
    end


  fun isize(): ISize =>
    """
    Convert `this` to a `ISize` value. On ILP32, that's a conversion to `I32`,
    and to `I64` for all other platforms.

    * NaN ⟶ 0
    * -∞ ⟶ `ISize.min_value()`
    * +∞ ⟶ `ISize.max_value()`
    * ]-1, +1[ ⟶ 0
    * Other values are saturated to the `ISize` range.
    """
    ifdef lp64 then
      i64().isize()
    elseif llp64 then
      i64().isize()
    elseif ilp32 then
      i32().isize()
    else
      i64().isize()
    end
  

  fun i128_unsafe(): I128 =>
    """
    Convert `this` to an `I128` value. If `this` has decimals, they
    are truncated. In case of overflow or special value, the result is undefined.
    """
    let e = _exponent.usize()
    if e <= 0 then
      return 0
    end

    var res: U128 = 0
    try
      res = _digits(e - 1)?.u128()
      res = res or (_digits(e - 2)?.u128() << 8)
      res = res or (_digits(e - 3)?.u128() << 16)
      res = res or (_digits(e - 4)?.u128() << 24)
      res = res or (_digits(e - 5)?.u128() << 32)
      res = res or (_digits(e - 6)?.u128() << 40)
      res = res or (_digits(e - 7)?.u128() << 48)
      res = res or (_digits(e - 8)?.u128() << 56)
      res = res or (_digits(e - 9)?.u128() << 64)
      res = res or (_digits(e - 10)?.u128() << 72)
      res = res or (_digits(e - 11)?.u128() << 80)
      res = res or (_digits(e - 12)?.u128() << 88)
      res = res or (_digits(e - 13)?.u128() << 96)
      res = res or (_digits(e - 14)?.u128() << 104)
      res = res or (_digits(e - 15)?.u128() << 112)
      res = res or (_digits(e - 16)?.u128() << 120)
    end

    if _sign then
      (-res).i128_unsafe()
    else
      res.i128_unsafe()
    end


  fun i64_unsafe(): I64 =>
    """
    Convert `this` to an `I64` value. If `this` has decimals, they
    are truncated. In case of overflow or special value, the result is undefined.
    """
    let e = _exponent.usize()
    if e <= 0 then
      return 0
    end

    var res: U64 = 0
    try
      res = _digits(e - 1)?.u64()
      res = res or (_digits(e - 2)?.u64() << 8)
      res = res or (_digits(e - 3)?.u64() << 16)
      res = res or (_digits(e - 4)?.u64() << 24)
      res = res or (_digits(e - 5)?.u64() << 32)
      res = res or (_digits(e - 6)?.u64() << 40)
      res = res or (_digits(e - 7)?.u64() << 48)
      res = res or (_digits(e - 8)?.u64() << 56)
    end

    if _sign then
      (-res).i64_unsafe()
    else
      res.i64_unsafe()
    end


  fun i32_unsafe(): I32 =>
    """
    Convert `this` to an `I32` value. If `this` has decimals, they
    are truncated. In case of overflow or special value, the result is undefined.
    """
    let e = _exponent.usize()
    if e <= 0 then
      return 0
    end

    var res: U32 = 0
    try
      res = _digits(e - 1)?.u32()
      res = res or (_digits(e - 2)?.u32() << 8)
      res = res or (_digits(e - 3)?.u32() << 16)
      res = res or (_digits(e - 4)?.u32() << 24)
    end

    if _sign then
      (-res).i32_unsafe()
    else
      res.i32_unsafe()
    end


  fun i16_unsafe(): I16 =>
    """
    Convert `this` to an `I16` value. If `this` has decimals, they
    are truncated. In case of overflow or special value, the result is undefined.
    """
    let e = _exponent.usize()
    if e <= 0 then
      return 0
    end

    var res: U16 = 0
    try
      res = _digits(e - 1)?.u16()
      res = res or (_digits(e - 2)?.u16() << 8)
    end

    if _sign then
      (-res).i16_unsafe()
    else
      res.i16_unsafe()
    end


  fun i8_unsafe(): I8 =>
    """
    Convert `this` to an `I8` value. If `this` has decimals, they
    are truncated. In case of overflow or special value, the result is undefined.
    """
    let e = _exponent.usize()
    if e <= 0 then
      return 0
    end

    var res: U8 = 0
    try res = _digits(e - 1)? end

    if _sign then
      (-res).i8()
    else
      res.i8()
    end


  fun ilong_unsafe(): ILong =>
    """
    Convert `this` to a `ILong` value without checking for overflow or
    special values. On LLP64 and ILP32, that's a conversion to `I32`, and
    to `I64` for all other platforms.
    """
    ifdef lp64 then
      i64_unsafe().ilong_unsafe()
    elseif llp64 then
      i32_unsafe().ilong_unsafe()
    elseif ilp32 then
      i32_unsafe().ilong_unsafe()
    else
      i64_unsafe().ilong_unsafe()
    end


  fun isize_unsafe(): ISize =>
    """
    Convert `this` to a `ISize` value without checking for overflow or
    special values. On ILP32, that's a conversion to `I32`, and to `I64`
    for all other platforms.
    """
    ifdef lp64 then
      i64_unsafe().isize_unsafe()
    elseif llp64 then
      i64_unsafe().isize_unsafe()
    elseif ilp32 then
      i32_unsafe().isize_unsafe()
    else
      i64_unsafe().isize_unsafe()
    end


    