// Multi-precision floating point numbers

use "debug"
use "../assertx"


type RoundingMode is (RoundingNearest | RoundingNegInf | RoundingPosInf |
                      RoundingZero | RoundingAwayZ | RoundingFaithful )
  """
  The rounding mode that must be applied to operations.

  * MPFloat_RNDN: round to nearest, with the even rounding rule (roundTiesToEven in IEEE 754).
  * MPFloat_RNDD: round toward negative infinity (roundTowardNegative in IEEE 754).
  * MPFloat_RNDU: round toward positive infinity (roundTowardPositive in IEEE 754).
  * MPFloat_RNDZ: round toward zero (roundTowardZero in IEEE 754).
  * MPFloat_RNDA: round away from zero.
  * MPFloat_RNDF: faithful rounding. This feature is currently experimental.

  For full definition and background: https://www.mpfr.org/mpfr-current/mpfr.html#Rounding
  and https://en.wikipedia.org/wiki/Rounding
  """


primitive RoundingNearest
    """
    Round to the nearest, with the even rounding rule.
    """
  fun apply(): I32 =>
    0

  fun string(): String =>
    "MPFloat_RNDN"


primitive RoundingNegInf
    """
    Round toward negative infinity.
    """
  fun apply(): I32 =>
    3

  fun string(): String =>
    "MPFloat_RNDD"


primitive RoundingPosInf
    """
    Round toward positive infinity.
    """
  fun apply(): I32 =>
    2

  fun string(): String =>
    "MPFloat_RNDU"


primitive RoundingZero
    """
    Round toward zero.
    """
  fun apply(): I32 =>
    1

  fun string(): String =>
    "MPFloat_RNDZ"


primitive RoundingAwayZ
    """
    Round away from zero.
    """
  fun apply(): I32 =>
    4

  fun string(): String =>
    "MPFloat_RNDA"


primitive RoundingFaithful
    """
    Faithfull rounding.
    """
  fun apply(): I32 =>
    5

  fun string(): String =>
    "MPFloat_RNDF"


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
  arithmetic; precision is limited only by `prec` bytes, not by F64.
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
    is guaranteed to be non-zero (normalised representation). For zero, the
    array may be all-zeros or empty.
    """

  let _base: U32 = 256
    """
    The base of the digits.
    """

  let _rounding: RoundingMode
    """
    The rounding mode of the `MPFloat`
    """

  new val _create(
    sgn: Bool,
    nan: Bool,
    inf: Bool,
    exp: I64,
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
    _exponent = exp
    _digits = digits
    _rounding = rnd


  new val create(prec: USize = 0, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a positive zero with `prec` bytes of digit capacity.

    This is the default constructor: `MPFloat()` gives a size-0 positive zero,
    and `MPFloat(n)` gives a size-n positive zero (all digits zero).
    """
    _sign = false
    _nan = false
    _inf = false
    _exponent = 0
    _digits = Array[U8].init(0, prec)
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


  new val from_f64(f: F64, prec: USize = 8, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` from the `F64` value `f` with `prec` bytes of
    precision (default 8, giving ~19 decimal digits) and a default rounding
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
      _digits = Array[U8].init(0, prec)
      _rounding = rnd
    else
      _sign = f < 0.0
      _nan = false
      _inf = false
      // Normalise: find frac in [1/256, 1) and exp such that |f| = frac×256^exp.
      let base = F64.from[U32](_base)
      var frac: F64 = f.abs()
      var exp: I64 = 0
      if frac >= 1.0 then
        while frac >= 1.0 do
          frac = frac / base
          exp = exp + 1
        end
      elseif frac < (1.0 / base) then
        while frac < (1.0 / base) do
          frac = frac * base
          exp = exp - 1
        end
      end
      _exponent = exp
      // Extract prec base-256 digits: d_i = floor(frac × 256), then remove it.
      let p: USize = prec.max(1)
      _digits = recover
        let d = Array[U8].init(0, p)
        var i: USize = 0
        var fr: F64 = frac
        while i < p do
          fr = fr * base
          let di: U8 = fr.u8()
          try d.update(i, di)? end
          fr = fr - di.f64()
          i = i + 1
        end
        d
      end
      _rounding = rnd
    end


  new val from_string(s: String = "",
                      prec: USize = 8,
                      base: U8 = 10,
                      rnd: RoundingMode = RoundingNearest) ? =>
    """
    Create a new `MPFloat` by parsing the string `s` with `prec` base-256
    digits of precision (default 8, ≈19 significant decimal digits). Raises
    an error if `s` is not a recognised floating-point representation.

    The `base` parameter (default 10) selects the numeral base; currently
    only base 10 is implemented. The `rnd` parameter (default nearest) is
    accepted for API compatibility with GMP; rounding support is a future TODO.

    Accepted formats:
    - `""`, `"0"`, `"0.0"`, `"+0"`, `"+0.0"` → +0
    - `"-0"`, `"-0.0"` → −0
    - `"nan"`, `"NaN"`, `"@NaN@"` → NaN
    - `"+inf"`, `"@Inf@"` → +∞
    - `"-inf"`, `"-@Inf@"` → −∞
    - `[+-]d+[.d*][e|E|@[+-]d+]` → finite decimal

    Leading and trailing whitespace is stripped. Unlike the `F64` path, the
    decimal mantissa is parsed directly into multi-precision arithmetic, so
    the full `prec` bytes of precision are exploited regardless of how many
    significant digits the string contains.
    """
    // Spaces are not significant.
    let st: String = s.clone() .> strip()

    if base != 10 then
      error  // TODO: implement parsing for bases other than 10
    end

    // Empty string → +0
    if st.size() == 0 then
      _sign = false
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, prec)
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

    if (st == "+inf") or (st == "@Inf@") then
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
      _digits = Array[U8].init(0, prec)
      _rounding = rnd
      return
    end

    if (st == "-0") or (st == "-0.0") then
      _sign = true
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, prec)
      _rounding = rnd
      return
    end

    // Multi-precision decimal parsing.
    //
    // We collect at most `sig_limit` significant decimal digits into n_mp
    // using Horner's method (n_mp = n_mp × 10 + digit), then scale by
    // 10^dec_exp where dec_exp = int_exp - frac_count + str_exp.
    //
    // Using guard precision p2 = prec + 2 throughout to suppress rounding
    // error in intermediate multiplications, truncating to prec at the end.
    //
    // log10(256^prec) = prec × log10(256) ≈ prec × 2.408 decimal digits;
    // add 4 guard digits to absorb Horner rounding.
    let p2: USize = prec + 2
    let sig_limit: USize = ((prec.f64() * 2.41).usize() + 4).max(1)

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
    let ten_mp: MPFloat = MPFloat.from_f64(10.0, p2, rnd)
    var n_mp: MPFloat = MPFloat.create(p2)
    var sig_count: USize = 0
    // int_exp: extra powers of 10 accumulated when skipping digits beyond
    //          sig_limit in the integer part.
    var int_exp: I64 = 0
    // frac_count: decimal places consumed (significant or leading-zero).
    var frac_count: I64 = 0
    var has_digit: Bool = false

    // Integer digits.
    while (pos < sz) and (st(pos)? >= '0') and (st(pos)? <= '9') do
      let d: U8 = st(pos)? - '0'
      has_digit = true
      if (d != 0) or (sig_count > 0) then
        if sig_count < sig_limit then
          n_mp = n_mp.mul(ten_mp)._trunc(p2)
                      .add(MPFloat.from_f64(d.f64(), p2, rnd))._trunc(p2)
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
      while (pos < sz) and (st(pos)? >= '0') and (st(pos)? <= '9') do
        let d: U8 = st(pos)? - '0'
        has_digit = true
        if (d != 0) or (sig_count > 0) then
          if sig_count < sig_limit then
            n_mp = n_mp.mul(ten_mp)._trunc(p2)
                        .add(MPFloat.from_f64(d.f64(), p2, rnd))._trunc(p2)
            sig_count = sig_count + 1
            frac_count = frac_count + 1
          end
          // Digits beyond sig_limit in the fractional part are ignored;
          // they fall below the requested precision.
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

    // Decimal exponent ('e', 'E', or '@' for base ≤ 10).
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
      while (pos < sz) and (st(pos)? >= '0') and (st(pos)? <= '9') do
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

    // True value = n_mp × 10^dec_exp.
    let dec_exp: I64 = (int_exp - frac_count) + str_exp

    // Scale by 10^|dec_exp| using fast (binary) exponentiation, keeping
    // only p2 significant base-256 digits at each step.
    if dec_exp != 0 then
      var scale: MPFloat = MPFloat.from_f64(1.0, p2, rnd)
      var sbase: MPFloat = ten_mp
      var sn: I64 = if dec_exp > 0 then dec_exp else -dec_exp end
      while sn > 0 do
        if (sn and 1) == 1 then
          scale = scale.mul(sbase)._trunc(p2)
        end
        sbase = sbase.mul(sbase)._trunc(p2)
        sn = sn / 2
      end
      n_mp =
        if dec_exp > 0 then
          n_mp.mul(scale)._trunc(p2)
        else
          n_mp.div(scale)._trunc(p2)
        end
    end

    let tmp: MPFloat = n_mp._trunc(prec)
    _sign = fsign
    _nan = tmp._nan
    _inf = tmp._inf
    _exponent = tmp._exponent
    _digits = tmp._digits
    _rounding = rnd


  new val from_mpint(n: MPInt, prec: USize = 8, rnd: RoundingMode = RoundingNearest) =>
    """
    Create a new `MPFloat` from the `MPInt` value `n` with `prec` bytes of
    precision (default 8, giving ~19 decimal digits) and rounding mode `rnd`.

    The conversion is exact up to the requested precision: the magnitude of
    `n` is represented without error as long as `n` fits within `prec` base-256
    digits (i.e. `|n| < 256^prec`); larger values are truncated to the `prec`
    most-significant bytes.

    Special cases:
    - Zero → `+0` (positive zero regardless of any sign on the MPInt zero).
    - Sign is preserved: a negative `MPInt` produces a negative `MPFloat`.

    Algorithm: call `MPInt.raw_digits()` to get the absolute value as
    a big-endian `Array[U8]` (each base-65536 word split into two bytes, MSW
    first, leading zeros stripped).  The result maps directly onto the MPFloat
    `_digits` layout:
    - `_exponent` = total byte count of the full magnitude (before truncation).
    - `_digits` = the first `prec` most-significant bytes (truncated).
    This is O(n) in the word count, vs O(n²) for the former string-based route.
    The rounding mode `rnd` is accepted for API compatibility; rounding support
    is a future TODO.
    """
    _rounding = rnd

    if n.is_zero() then
      _sign = false
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, prec)
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

    // Keep only the `prec` most-significant bytes.
    _digits = recover
      let keep: USize = total.min(prec)
      let d = Array[U8].create(keep)
      var i: USize = 0
      while i < keep do
        try d.push(mag(i)?) end
        i = i + 1
      end
      d
    end


  new min_normalized(prec: USize = 8, rnd: RoundingMode = RoundingNearest) =>
    """
    The smallest normalized floating point number.

    TODO
    """
    _sign = false
    _nan = false
    _inf = false
    _exponent = 0
    _digits = Array[U8].init(0, prec)
    _rounding = rnd


  new epsilon(prec: USize = 8, rnd: RoundingMode = RoundingNearest) =>
    """
    `epsilon` is the "floating-point precision", that is the smallest positive
    floating-point number such that `1.0 + epsilon != 1.0`. In infinite-precision
    or `MPFloat`, `epsilon` can be made as small as desired by increasing `prec`,
    and `0` becomes this limit value.

    TODO
    """
    _sign = false
    _nan = false
    _inf = false
    _exponent = 0
    _digits = Array[U8].init(0, prec)
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


  new val pi(prec: USize = 64, rnd: RoundingMode = RoundingNearest) =>
    """
    The pi constant, calculated with the specified `size` accuracy (number of
    base-256 digits, default 64 giving ~154 decimal digits). The rounding
    mode `rnd` is not used yet (TODO).

    Uses the Borwein quadratic convergence algorithm:

    ```
    x_0 = √2,  π_0 = 2 + √2,  y_0 = ⁴√2

    x_{n+1} = (√x_n + 1/√x_n) / 2
    π_{n+1} = π_n × (x_{n+1} + 1) / (y_n + 1)
    y_{n+1} = (y_n × √x_{n+1} + 1/√x_{n+1}) / (y_n + 1)
    ```
    """
    _sign = false
    _nan = false
    _inf = false
    _exponent = 1
    _rounding = rnd

    let prc: USize = prec + 2

    // Build two = 2.0  (first digit = 2, rest = 0, exponent = 1)
    let two: MPFloat = MPFloat._create(false, false, false, 1, recover
        let d = Array[U8].init(0, prc)
        try d.update(0, 2)? end
        d
      end, rnd)

    // x_0 = sqrt(2)
    var x_n: MPFloat = two.sqrt()

    // pi_0 = 2 + sqrt(2)  (digit_shl(1) removes the integer carry digit)
    var pi_n: MPFloat = (two + x_n).digit_shl(1)

    // y_0 = sqrt(sqrt(2))
    var y_n: MPFloat = x_n.sqrt()

    // Recurrence — bound all intermediates to prc digits to prevent growth.
    var iters: USize = 0
    let max_iters: USize = prc * 4
    try
      while iters < max_iters do
        iters = iters + 1

        // x_{n+1} = (√x_n + 1/√x_n) / 2
        let s: MPFloat = x_n.sqrt()
        (let t, _) = ((s + s.inv())._trunc(prc))._short_div(2)
        x_n = t.digit_shl(1)

        // y_{n+1} = (y_n × √x_{n+1} + 1/√x_{n+1}) / (y_n + 1)
        let u: MPFloat = x_n.sqrt()
        let v: MPFloat =
          ((((y_n * u)._trunc(prc)).digit_shl(1)) + u.inv())._trunc(prc).digit_shl(1)

        // y_n ← v / (y_n + 1)
        y_n = v / y_n._inc_first()
        y_n = y_n.digit_shl(1)._trunc(prc)

        // π_{n+1} = π_n × (x_{n+1} + 1) / (y_n + 1)
        let xp1: MPFloat = x_n._inc_first()
        let w: MPFloat = xp1 / y_n
        pi_n = ((pi_n * (w.digit_shl(1)))._trunc(prc)).digit_shl(1)

        // Convergence check: w = x_{n+1} / y_n → 1 when converged.
        // In base-256 that means w._digits(0) = 1 and w._digits(1..] = 0.
        let m: U8 = w._digits(0)? - 1
        var i: USize = 1
        while i < (prec + 1) do
          if w._digits(i)? != m then
            break
          end
          i = i + 1
        end
        if i == (prec + 1) then
          break
        end
      end
    else
      // Index out of bounds — should not happen for valid inputs.
      try
        Assert(false, "[M]PFloat.pi] unexpected error during Newton iteration", true)?
      end
    end

    _digits = pi_n._digits


  //- Internal helpers ---------------------------------------------------------

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
          let bi: U8 =
            if (col >= shift) and ((col - shift) < nb) then
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
          let bi: U16 =
            if (col >= shift) and ((col - shift) < nb) then
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
      return MPFloat.create(na)
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
      MPFloat.create(_size(), _rounding)
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

    Algorithm:
    1. Estimate the decimal exponent `dec_exp` via F64.
    2. Scale `|this|` to the interval `[1, 10)` by multiplying by
       `10^(1 − dec_exp)` (for dec_exp ≤ 1) or dividing by
       `10^(dec_exp − 1)` (for dec_exp > 1) using binary fast-exponentiation.
       When dec_exp = 1 no scaling is needed, which avoids lossy division
       for exactly-representable integer values (2.0, 3.0, …).
    3. Verify the result is in `[1, 10)`; adjust `dec_exp` by ±1 if the
       F64 estimate was off.
    4. Extract `n_dec` decimal digits: the first digit is `_digits(0)` of
       the scaled value (which has `_exponent = 1`). Subsequent digits use
       `digit_shl(1)` to remove the current integer digit followed by
       `_short_mul(10)` to bring the next digit into `_digits(0)`.
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

    let prec: USize = _size()
    let p2: USize = prec + 2
    // n_dec: significant decimal digits to produce.
    // log₁₀(256^prec) = prec × log₁₀(256) ≈ prec × 2.408; add 2 guard digits.
    let n_dec: USize = ((prec.f64() * 2.41).usize() + 2).max(1)

    // ---- Step 1: approximate decimal exponent via F64. ----------------------
    // fg_hi ≈ d₀ + d₁/256 + … ∈ [1, 256) (= F × 256 where F = 0.d₀d₁…).
    let ng: USize = prec.min(4)
    var fg_hi: F64 = 0.0
    let bas = _base.f64()
    try
      var k: USize = ng
      repeat
        k = k - 1
        fg_hi = (fg_hi / bas) + _digits(k)?.f64()
      until k == 0 end
    end
    // log₁₀(|this|) = log₁₀(fg_hi) + (_exponent − 1) × log₁₀(256).
    let log10_256: F64 = bas.log10() // 2.408240961 for _base = 256
    let log10_val: F64 = fg_hi.log10() + ((_exponent.f64() - 1.0) * log10_256)
    var dec_exp: I64 = log10_val.floor().i64() + 1

    // ---- Exact-integer fast path -------------------------------------------
    // If all bytes from position _exponent onward are zero, |this| is an exact
    // integer.  Use integer long-division by 10 to extract decimal digits;
    // this avoids the repeating base-256 fraction problem (e.g. 2.56 in
    // base-256 is periodic, so _short_div(10) gives a value slightly < 2.56
    // which then rounds the third digit down from 6 to 5).
    // The fast path is only valid when _exponent ≤ prec, meaning ALL integer
    // bytes fit within the stored digit array.  When _exponent > prec the
    // number is a large integer that overflows the stored precision; the
    // stored bytes represent value / 256^(_exponent − prec), not the full
    // integer, so the fast path would give a wrong result smaller by
    // 256^(_exponent − prec).  Fall through to the regular scaling path in
    // that case.
    let int_bytes: USize = if _exponent > 0 then _exponent.usize().min(prec) else 0 end
    var is_exact_int: Bool = (int_bytes > 0) and (_exponent.usize() <= prec)
    if is_exact_int then
      try
        var i: USize = int_bytes
        while i < prec do
          if _digits(i)? != 0 then
            is_exact_int = false
            break
          end
          i = i + 1
        end
      end
    end

    if is_exact_int then
      // Pre-pass: compute the exact decimal digit count via long-division on a
      // copy of the integer bytes.  This corrects the F64 log10 approximation
      // which can be off by 1 at boundary values close to 10^k (e.g. 10^20
      // gives log10_val ≈ 19.9999… → dec_exp = 20 instead of the correct 21).
      var pre_ibuf: Array[U8] = Array[U8].init(0, int_bytes.max(1))
      try
        var i: USize = 0
        while i < int_bytes do
          pre_ibuf(i)? = _digits(i)?
          i = i + 1
        end
      end
      var pre_has_nz: Bool = false
      try
        var i: USize = 0
        while i < int_bytes do
          if pre_ibuf(i)? != 0 then pre_has_nz = true; break end
          i = i + 1
        end
      end
      var n_dec_exact: USize = 0
      while pre_has_nz do
        var remain: U16 = 0
        var q_nz: Bool = false
        try
          var i: USize = 0
          while i < pre_ibuf.size() do
            let xv: U16 = remain.shl(8) + pre_ibuf(i)?.u16()
            let q: U16 = xv / 10
            remain = xv - (q * 10)
            pre_ibuf(i)? = q.u8()
            if q != 0 then q_nz = true end
            i = i + 1
          end
        end
        n_dec_exact = n_dec_exact + 1
        pre_has_nz = q_nz
      end
      dec_exp = n_dec_exact.max(1).i64()

      let sign_prefix: String = if _sign then "-" else "" end
      // All mutable computation is inside recover so Pony's capability rules
      // are satisfied.  _digits (Array[U8]) and _exponent (I64) are val
      // and are readable inside the recover block.
      let int_result: String iso = recover
        // Copy integer digit bytes into a mutable buffer for long-division.
        let ibuf: Array[U8] = Array[U8].init(0, int_bytes.max(1))
        try
          var i: USize = 0
          while i < int_bytes do
            ibuf(i)? = _digits(i)?
            i = i + 1
          end
        end

        // Long-divide ibuf (big-endian base-256) by 10, collecting remainders
        // (= decimal digits, least-significant first).
        let drev: Array[U8] = Array[U8].create()
        var has_nz: Bool = false
        try
          var i: USize = 0
          while i < ibuf.size() do
            if ibuf(i)? != 0 then
              has_nz = true
              break
            end
            i = i + 1
          end
        end

        while has_nz do
          var remain: U16 = 0
          var q_nonzero: Bool = false
          try
            var i: USize = 0
            while i < ibuf.size() do
              let xv: U16 = remain.shl(8) + ibuf(i)?.u16()
              let q: U16 = xv / 10
              remain = xv - (q * 10)
              ibuf(i)? = q.u8()
              if q != 0 then
                q_nonzero = true
              end
              i = i + 1
            end
          end
          drev.push(remain.u8())
          has_nz = q_nonzero
        end

        if drev.size() == 0
          then drev.push(0)
        end

        // Build string: sign + digits most-significant first + zero-pad to n_dec
        // total digits (so string() uses fixed rather than scientific notation).
        let total: USize = n_dec.max(drev.size())
        let buf = String.create(sign_prefix.size() + total)
        buf.append(sign_prefix)
        var i: USize = drev.size()
        while i > 0 do
          i = i - 1
          try buf.push('0' + drev(i)?) end
        end

        i = drev.size()
        while i < n_dec do
          buf.push('0')
          i = i + 1
        end
        buf
      end

      // Use drev.size() as dec_exp would be, but we need the actual digit count.
      // dec_exp (from F64) is correct for integers within F64 precision; use it.
      return (consume int_result, dec_exp, false)
    end

    // ---- Step 2: scale |this| to [1, 10). ----------------------------------
    // s_val = |this| × 10^(1 − dec_exp) so that s_val ∈ [1, 10).
    // For dec_exp = 1: multiply by 10^0 = 1, i.e. no scaling (avoids a
    //   lossy div for exactly-representable values like 2.0, 3.0, etc.).
    // For dec_exp < 1: multiply by 10^(1 − dec_exp).
    // For dec_exp > 1: divide by 10^(dec_exp − 1).
    // Binary exponentiation keeps the cost O(log |dec_exp|) multiplications.
    let abs_val = MPFloat._create(false, false, false, _exponent, _digits, _rounding)
    let ten_mp = MPFloat.from_f64(10.0, p2, _rounding)
    // TODO: Evaluate rounding
    let new_rnd = _rounding
    var s_val: MPFloat =
      if dec_exp == 1 then
        abs_val._trunc(p2)
      elseif dec_exp < 1 then
        // multiply by 10^(1 − dec_exp)
        var scale = MPFloat.from_f64(1.0, p2, _rounding)
        var sb: MPFloat = ten_mp
        var sn: I64 = 1 - dec_exp
        while sn > 0 do
          if (sn and 1) == 1 then
            scale = scale.mul(sb)._trunc(p2)
          end
          sb = sb.mul(sb)._trunc(p2)
          sn = sn / 2
        end
        abs_val.mul(scale)._trunc(p2)
      else
        // dec_exp > 1: exact short-division by 10 repeated (dec_exp − 1)
        // times.  Using _short_div avoids Newton's-method error.
        //
        // Precision note: 256/100 = 2.56 repeats in base-256 (period 5),
        // so each _short_div step may lose one byte of precision when a
        // leading zero is stripped.  We compensate by starting with
        // (dec_exp − 1) extra bytes so that after at most (dec_exp − 1)
        // normalisation steps the digit array still has ≥ p2 bytes.
        let extra: USize = (dec_exp - 1).usize()
        let work_size: USize = p2 + extra
        // Extend abs_val to work_size bytes, zero-padding on the right.
        // For exact integers this is mathematically correct; for others it
        // adds a small tail of zeros which is fine for the guard-digit role.
        var sv: MPFloat =
          if abs_val._size() >= work_size then
            abs_val._trunc(work_size)
          else
            let d_ext: Array[U8] val = recover
              let a = Array[U8].init(0, work_size)
              abs_val._digits.copy_to(a, 0, 0, abs_val._size())
              a
            end
            MPFloat._create(false, false, false, _exponent, d_ext, _rounding)
          end

        var sn: I64 = dec_exp - 1
        while sn > 0 do
          (let sv2, _) = sv._short_div(10)
          sv = sv2
          // Normalise: strip any leading zero digit(s) produced by the
          // division, adjust the exponent, and zero-pad on the RIGHT to
          // keep the array at the same size (preventing precision loss).
          var i: USize = 0
          try
            while (i < sv._size()) and (sv._digits(i)? == 0) do
              i = i + 1
            end
          end

          if i > 0 then
            let total_sz: USize = sv._size()
            let new_exp: I64 = sv._exponent - i.i64()
            let new_d: Array[U8] val = recover
              let a = Array[U8].init(0, total_sz) // zero-filled
              sv._digits.copy_to(a, i, 0, total_sz - i)
              a
            end
            sv = MPFloat._create(false, false, false, new_exp, new_d, new_rnd)
          end
          sn = sn - 1
        end
        sv
      end

    // ---- Step 3: verify s_val ∈ [1, 10); correct dec_exp if off by ±1. ----
    let one = MPFloat.from_f64(1.0, p2, _rounding)
    match s_val._cmp_mag(ten_mp)
    | Greater | Equal =>
      s_val = s_val.div(ten_mp)._trunc(p2)
      dec_exp = dec_exp + 1
    end
    match s_val._cmp_mag(one)
    | Less =>
      s_val = s_val.mul(ten_mp)._trunc(p2)
      dec_exp = dec_exp - 1
    end

    // ---- Step 4: extract n_dec decimal digits. ------------------------------
    // After step 3, s_val ∈ [1, 10): _exponent = 1, _digits(0) ∈ [1, 9].
    // First digit: _digits(0) directly.
    // Subsequent digits: digit_shl(1) removes the top digit (leaving the
    //   fractional remainder scaled ×256 in _exponent=1 form), then
    //   _short_mul(10) yields the next digit in _digits(0) as
    //   floor(frac_256 × 10 / 256) ∈ [0, 9].
    let sign_prefix: String = if _sign then "-" else "" end
    let result: String iso = recover
      let buf = String.create(sign_prefix.size() + n_dec)
      buf.append(sign_prefix)

      // First decimal digit from the integer part of s_val.
      buf.push('0' + (try s_val._digits(0)? else 0 end))
      // Remaining digits via repeated ×10 extraction.
      var frac: MPFloat = s_val.digit_shl(1)
      var i: USize = 1
      while i < n_dec do
        frac = frac._short_mul(10)
        buf.push('0' + (try frac._digits(0)? else 0 end))
        frac = frac.digit_shl(1)
        i = i + 1
      end
      buf
    end
    (consume result, dec_exp, false)


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
      MPFloat.create(_size(), _rounding)
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
    // pow2 must be a power of 2 and >= res_size to avoid circular aliasing.
    let pow2: USize = res_size.next_pow2()

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
    let nan: MPFloat = MPFloat.nan_val()
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
      return (MPFloat.create(_size(), _rounding), MPFloat._create(_sign, false, false, _exponent, _digits, _rounding))
    end
    if that.is_zero() then
      if is_zero() then
        return (nan, nan)
      end
      return (MPFloat._create(_sign != that._sign, false, true, 0, Array[U8].create(), _rounding), nan)
    end
    if is_zero() then
      let z: MPFloat = MPFloat.create(_size(), _rounding)
      return (z, z)
    end
    var q: MPFloat = div(that)._trunc_frac()
    var r: MPFloat = sub(q.mul(that))
    // Post-correction: Newton's inv() always undershoots, so q may be 1 too
    // small when this/that is an exact integer (e.g. 10/5 → q≈1.999 → trunc=1).
    // If |r| ≥ |that|, increment |q| by 1 and adjust r.  At most one step is
    // needed because the Newton error is bounded by one base-256 ULP.
    if (not r.is_zero()) and (not r.abs().lt(that.abs())) then
      let one: MPFloat = MPFloat.from_f64(1.0, _size(), _rounding)
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
        return q.sub(MPFloat.from_f64(1.0, _size(), _rounding))
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
    let g: MPFloat = MPFloat._create(false, false, false, 1, _digits, _rounding)

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
    var res: MPFloat = MPFloat.from_f64(1.0 / fg, size, _rounding)

    // Newton refinement: y_{n+1} = y_n × (2 − G × y_n).
    let two: MPFloat = MPFloat.from_f64(2.0, size, _rounding)
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
      q.sub_unsafe(MPFloat.from_f64(1.0, _size(), _rounding))
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
    let h: MPFloat = MPFloat._create(false, false, false, h_exp, _digits, _rounding)

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
    var res = MPFloat.from_f64(1.0 / fg_h.sqrt(), size, _rounding)

    // Newton refinement: y_{n+1} = y_n × (3 − H × y_n²) / 2.
    let three = MPFloat.from_f64(3.0, size, _rounding)
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
      return MPFloat.create(_size(), _rounding)
    end
    let e: USize = _exponent.usize()
    if e >= _size() then
      return MPFloat._create(_sign, false, false, _exponent, _digits, _rounding)
    end
    let new_digits: Array[U8] val = recover
      let d = Array[U8].create(e)
      var i: USize = 0
      while i < e do
        try d.push(_digits(i)?) end
        i = i + 1
      end
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
      return t.sub(MPFloat.from_f64(1.0, _size(), _rounding))
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
      return t.add(MPFloat.from_f64(1.0, _size(), _rounding))
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
    let half: MPFloat = MPFloat.from_f64(0.5, _size(), _rounding)
    if _sign then
      return sub(half).ceil()
    end
    add(half).floor()

