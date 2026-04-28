// Multi-precision floating point: pure representation layer.

use "../assertx"
use "../formatx"

use "collections"
use "debug"


class val MPFRep
  """
  The pure representation of an arbitrary-precision floating-point number.

  `MPFRep` stores the sign, special-value flags, base-256 exponent and mantissa
  digits of a number. The type is `class val` (globally immutable) so that instances
  can be freely shared across actors and stored in collections.

  Representation invariant (for a finite, non-zero value):

    value = (−1)^_sign × 0.d₀d₁d₂… × 256^_exponent

  where `d₀…d_{n−1}` are base-256 digits stored big-endian in `_digits`
  (most-significant digit at index 0), and `d₀ ≠ 0` (normalised form).
  Zero has an all-zero digit array. `_exponent` is a signed base-256 exponent.

  Special values: NaN (`_nan = true`) and ±∞ (`_inf = true`) are supported,
  mirroring IEEE 754 / GMP/MPFR conventions.

  The precision is encoded implicitly in `_digits.size()`: a representation
  with `p` bytes has approximately `p × log₁₀(256) ≈ p × 2.408` decimal
  significant digits. Construction routines accept a byte count `p_bytes`
  rather than a bit count; the caller is responsible for the conversion
  `p_bytes = ceil(prec_bits / 8)`.
  """

  let _sign: Bool
    """
    Sign bit: `true` when the value is strictly negative or is −0.
    Ignored when `_nan` is `true`.
    """

  let _nan: Bool
    """
    `true` when the number is Not-a-Number. When set, `_inf`, `_sign` and
    `_digits` are irrelevant.
    """

  let _inf: Bool
    """
    `true` when the number is ±∞. Sign is given by `_sign`.
    Ignored when `_nan` is `true`.
    """

  let _exponent: I64
    """
    Base-256 exponent `e` such that `|value| = 0.d₀d₁… × 256^e`.
    For the common arithmetic operations `e = 1` is used internally
    (the first digit `d₀` is the integer part).
    """

  let _digits: Array[U8] val
    """
    Big-endian base-256 mantissa digits (index 0 = most significant).
    For a finite non-zero number `_digits(0) ≠ 0` (normalised).
    For zero the array may be all-zeros or empty.
    """


  new val _create(
    sgn: Bool,
    nan: Bool,
    inf: Bool,
    expn: I64,
    digits: Array[U8] val)
  =>
    """
    Private canonical constructor. All constructors and internal operations
    produce their result through this entry point.
    """
    _sign = sgn
    _nan = nan
    _inf = inf
    _exponent = expn
    _digits = digits


  new val create(p_bytes: USize = 14) =>
    """
    Create a positive zero with `p_bytes` base-256 mantissa bytes
    (default 14, corresponding to the 112-bit F128 significand).

    `MPFRep(14)` gives a 14-byte positive zero; `MPFRep(0)` gives a
    precision-free zero (empty digit array).
    """
    _sign = false
    _nan = false
    _inf = false
    _exponent = 0
    _digits = Array[U8].init(0, p_bytes)


  new val nan_val() =>
    """
    Create a Not-a-Number representation.
    """
    _sign = false
    _nan = true
    _inf = false
    _exponent = 0
    _digits = Array[U8].create()


  new val inf_val(positive: Bool = true) =>
    """
    Create an infinite representation. Pass `positive = false` for −∞.
    """
    _sign = not positive
    _nan = false
    _inf = true
    _exponent = 0
    _digits = Array[U8].create()


  new val from_f64(f: F64, p_bytes: USize = 14) =>
    """
    Create a new `MPFRep` from the `F64` value `f` using `p_bytes` base-256
    mantissa bytes (default 14 ≈ 112 bits ≈ 33 decimal digits).

    Special values (NaN, ±∞, ±0) are preserved. The conversion normalises
    `f` so that `0.d₀d₁… × 256^_exponent` with `d₀ ≠ 0` for non-zero values.
    """
    if f.nan() then
      _sign = false
      _nan = true
      _inf = false
      _exponent = 0
      _digits = Array[U8].create()
    elseif f.infinite() then
      _sign = f < 0.0
      _nan = false
      _inf = true
      _exponent = 0
      _digits = Array[U8].create()
    elseif f == 0.0 then
      _sign = (f.bits() == 0x8000000000000000)
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, p_bytes)
    else
      _sign = f < 0.0
      _nan = false
      _inf = false

      // Normalise exactly using frexp.
      // value = m × 2^e = 0.d0d1… × 256^expn = frac × 2^(8×expn)
      // We want 1/256 ≤ frac < 1.
      // Since 0.5 ≤ m < 1, e − 8×expn must be in [−7, 0].
      (let m, let e_u32) = f.abs().frexp()
      let e = e_u32.i32().i64()
      let expn = (e.f64() / 8.0).ceil().i64()
      let shift = e - (expn * 8)
      var frac = m.f64() * F64(2).pow(shift.f64())
      _exponent = expn

      let base: F64 = 256.0
      _digits = recover
        let d = Array[U8].init(0, p_bytes)
        var i: USize = 0
        while i < p_bytes do
          frac = frac * base
          let di: U8 = frac.u8()
          try d.update(i, di)? end
          frac = frac - di.f64()
          i = i + 1
        end
        d
      end
    end


  new val from_f32(f: F32, p_bytes: USize = 14) =>
    """
    Create a new `MPFRep` from the `F32` value `f` using `p_bytes` base-256
    mantissa bytes (default 14 ≈ 112 bits ≈ 33 decimal digits).

    Special values (NaN, ±∞, ±0) are preserved. The conversion normalises
    `f` so that `0.d₀d₁… × 256^_exponent` with `d₀ ≠ 0` for non-zero values.
    """
    if f.nan() then
      _sign = false
      _nan = true
      _inf = false
      _exponent = 0
      _digits = Array[U8].create()
    elseif f.infinite() then
      _sign = f < 0.0
      _nan = false
      _inf = true
      _exponent = 0
      _digits = Array[U8].create()
    elseif f == 0.0 then
      _sign = (f.bits() == 0x80000000)
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, p_bytes)
    else
      _sign = f < 0.0
      _nan = false
      _inf = false

      (let m, let e_u32) = f.abs().frexp()
      let e = e_u32.i32().i64()
      let expn = (e.f32() / 8.0).ceil().i64()
      let shift = e - (expn * 8)
      var frac = m * F32(2).pow(shift.f32())
      _exponent = expn

      let base: F32 = 256.0
      _digits = recover
        let d = Array[U8].init(0, p_bytes)
        var i: USize = 0
        while i < p_bytes do
          frac = frac * base
          let di: U8 = frac.u8()
          try d.update(i, di)? end
          frac = frac - di.f32()
          i = i + 1
        end
        d
      end
    end


  new val from_mpint(n: MPInt, p_bytes: USize = 14) =>
    """
    Create a new `MPFRep` from the `MPInt` value `n` using `p_bytes` base-256
    mantissa bytes (default 14).

    The conversion is exact up to the requested precision: the magnitude of
    `n` is represented without error as long as `n` fits within `p_bytes × 8`
    bits. Larger values are truncated to `p_bytes` bytes.

    Special cases:
    - Zero → `+0` (positive zero regardless of any sign on the MPInt zero).
    - Sign is preserved: a negative `MPInt` produces a negative `MPFRep`.

    Algorithm: `MPInt.raw_digits()` provides the absolute value as a big-endian
    `Array[U8]` (each base-65536 word split into two bytes, MSW first, leading
    zeros stripped). The result maps directly onto the `MPFRep` layout:
    `_exponent` = total byte count of the full magnitude (before truncation).
    This is O(n) in the word count.
    """
    if n.is_zero() then
      _sign = false
      _nan = false
      _inf = false
      _exponent = 0
      _digits = Array[U8].init(0, p_bytes)
      return
    end

    _sign = n.is_negative()
    _nan  = false
    _inf  = false

    let mag = n.raw_digits()
    let total: USize = mag.size()
    _exponent = total.i64()

    _digits = recover
      let keep: USize = total.min(p_bytes)
      let d = Array[U8].create(keep)
      mag.copy_to(d, 0, 0, keep)
      d
    end


  new val from_mpfloat_rep(f: MPFRep, p_bytes: USize = 14) =>
    """
    Create a new `MPFRep` whose value equals `f` but with `p_bytes` base-256
    mantissa bytes (default 14). Useful for changing the working precision of
    an existing representation.

    The result is truncated to `p_bytes` (no rounding — rounding is the
    responsibility of `MPFContext`).
    """
    _sign = f._sign
    _nan = f._nan
    _inf = f._inf
    _exponent = f._exponent

    _digits = recover
      let size = p_bytes.min(f._size())
      let d = Array[U8].init(0, p_bytes)
      f._digits.copy_to(d, 0, 0, size)
      d
    end


  new val from_ulong(n: ULong, p_bytes: USize = 14) =>
    """
    Create a new `MPFRep` whose value equals the unsigned integer `n` with
    `p_bytes` base-256 mantissa bytes (default 14).
    """
    _sign = false
    _nan = false
    _inf = false

    let base: ULong = ULong(256)
    var q: ULong = n

    _digits = recover
      let d: Array[U8] = Array[U8].create(p_bytes)
      while q >= base do
        (q, let r) = q.divrem(base)
        d.push(r.u8())
      end
      d.push(q.u8())
      d.reverse_in_place()
      d
    end
    _exponent = _digits.size().i64()


  new \do_not_use\ val min_normalized(p_bytes: USize = 14) =>
    """
    The smallest normalised representation.

    As `MPFRep` has no notion of a minimum normalised value in the IEEE 754
    sense (the exponent is unbounded), the smallest normalised `MPFRep` is
    defined as `+0.0`.
    """
    _sign = false
    _nan = false
    _inf = false
    _exponent = 0
    _digits = Array[U8].init(0, p_bytes)


  new val epsilon(p_bytes: USize = 14) =>
    """
    The machine epsilon for the given precision: the smallest positive value ε
    such that `1 + ε ≠ 1` in base 256, i.e. `ε = 256^(1 − p_bytes)`.

    For the default precision of 14 bytes (112 bits, ≈33 significant decimal
    digits): `ε ≈ 4.93 × 10^{−32}`.

    This is the base-256 analogue of `F64.epsilon() ≈ 2.22 × 10^{−16}`.
    """
    _sign = false
    _nan = false
    _inf = false
    // ε = 0.1_{256} × 256^{2−p_bytes} = 256^{1−p_bytes}
    let pb: USize = p_bytes.max(1)
    _exponent = 2 - pb.i64()
    _digits = recover
      let a = Array[U8].init(0, pb)
      try a(0)? = 1 end
      a
    end


  new val min_value() =>
    """
    The minimum representable value is −∞ (`-inf`).
    """
    _sign = true
    _nan = false
    _inf = true
    _exponent = 0
    _digits = Array[U8].create()


  new val max_value() =>
    """
    The maximum representable value is +∞ (`+inf`).
    """
    _sign = false
    _nan = false
    _inf = true
    _exponent = 0
    _digits = Array[U8].create()


  //- Size and digit helpers ---------------------------------------------------

  fun _size(): USize =>
    """
    Number of base-256 digits (bytes) in the mantissa.
    """
    _digits.size()


  fun _lowb(a: U16): U8 =>
    """
    Extract the lower byte of the `U16` value `a`.
    """
    a.u8()


  fun _highb(a: U16): U8 =>
    """
    Extract the upper byte of the `U16` value `a` (bits 15–8).
    """
    a.shr(8).u8()
    
    
  fun _addc(a: U8, b: U8, c: U16): (U8, U16) =>
    """
    Single-digit addition `a + b + c`, returning `(digit, carry)`.
    The carry is 0 or 1 and propagates to the next more-significant column.
    """
    let r: U16 = a.u16() + b.u16() + c
    (r.u8(), r.shr(8))


  fun _subc(a: U8, b: U8, c: U16): (U8, U16) =>
    """
    Single-digit subtraction `a − b` with borrow `c`, returning
    `(digit, new_borrow)`. Uses unsigned-complement arithmetic so that the
    borrow is encoded in the high bit of a `U16`.
    """
    let r: U16 = ((a.max_value().u16() + a.u16()) - b.u16()) + c.shr(8)
    (r.u8(), r.shr(8))


  fun _inc_first(): MPFRep =>
    """
    Return a new `MPFRep` identical to `this` but with the most-significant
    digit (`_digits(0)`) incremented by 1. No carry propagation is performed;
    the caller must ensure `_digits(0) < 255`. Used internally by Newton
    iteration updates of the form `x ← x + 1`.
    """
    let size = _size()
    let d: Array[U8] val = recover
      let a = Array[U8].init(0, size)
      _digits.copy_to(a, 0, 0, size)
      try a.update(0, (try a(0)? else 0 end) + 1)? end
      a
    end
    MPFRep._create(_sign, _nan, _inf, _exponent, d)


  fun _trunc(n: USize): MPFRep =>
    """
    Return a new `MPFRep` containing only the `n` most-significant base-256
    digits of `this`. Used internally to bound intermediate results in Newton
    iterations so that successive multiplications do not grow the array without
    limit.

    This is a pure truncation (no rounding). Output rounding is the
    responsibility of `MPFContext._round_to`.
    """
    let s: USize = _size().min(n)
    let d: Array[U8] val = recover
      let a = Array[U8].init(0, s)
      _digits.copy_to(a, 0, 0, s)
      a
    end
    MPFRep._create(_sign, _nan, _inf, _exponent, d)


  fun _neg_comp(): MPFRep =>
    """
    Return the byte-level two's-complement negation of `_digits`. This is
    NOT a mathematical float negation — it is a fixed-point unsigned complement
    used internally by Newton iterations to compute `2 − x` in unsigned
    arithmetic. Use `neg()` (on `MPFloat`) for the mathematical sign flip.
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
    MPFRep._create(_sign, _nan, _inf, _exponent, d)


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

    - NaN → `false`. ±∞ → `false`.
    - ±0 → `true` (zero is an integer).
    - Any value whose base-256 representation has no non-zero fractional
      bytes (bytes at index ≥ `_exponent`) → `true`.
    - Any value with at least one non-zero fractional byte → `false`.

    Equivalent to `is_finite() and not _has_frac()`.
    """
    is_finite() and not _has_frac()


  fun _has_frac(): Bool =>
    """
    Return `true` if `this` has at least one non-zero fractional base-256 byte,
    i.e. `this ≠ trunc(this)`.

    - NaN or ±∞ → `false` (no fractional part defined for non-finite values).
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
    try
      for i in Range(e, _size()) do
        if _digits(i)? != 0 then
          return true
        end
      end
    end
    false


  fun _trunc_frac(): MPFRep =>
    """
    Return `this` with the fractional part removed (truncated toward zero).
    Internal alias equivalent to `trunc()` on `MPFloat`. Called by division
    and floored-division routines to obtain the integer part.
    """
    if not is_finite() then
      return MPFRep._create(_sign, _nan, _inf, _exponent, _digits)
    end
    if is_zero() then
      return MPFRep._create(_sign, false, false, 0, _digits)
    end
    if _exponent <= 0 then
      return MPFRep.create(_size())
    end
    let e: USize = _exponent.usize().min(_size())
    let d: Array[U8] val = recover
      let a = Array[U8].init(0, _size())
      _digits.copy_to(a, 0, 0, e)
      a
    end
    MPFRep._create(_sign, false, false, _exponent, d)


  //- Accessors ----------------------------------------------------------------

  fun exponent(): I64 =>
    """
    Return the base-256 exponent `e` such that `|value| = 0.d₀d₁… × 256^e`.

    The integer part of the value occupies the first `e` bytes of `raw_digits()`
    (when `e > 0`); bytes at position `e` and beyond represent the fractional
    part. A non-positive exponent means the value is strictly in (−1, 1).
    """
    _exponent


  fun raw_digits(): Array[U8] val =>
    """
    Return the internal base-256 mantissa as a big-endian `Array[U8]`.

    Together with `exponent()` this determines the magnitude exactly:
    `|value| = 0.d₀d₁… × 256^exponent()`. The first byte `d₀` is always
    non-zero for finite non-zero values.

    This accessor is intended for low-level conversions such as
    `MPInt.from_mpfloat`; prefer the high-level API for ordinary use.
    """
    _digits


  fun sign_bit(): Bool =>
    """
    Return the sign bit: `true` for negative values or −0, `false` otherwise.
    Always `false` for NaN.
    """
    _sign and not _nan


  //- Magnitude arithmetic (structural, no rounding) --------------------------

  fun _cmp_mag(that: MPFRep): Compare =>
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
    // the longer operand are all zero. If so the magnitudes are equal.
    if na > nb then
      var all_zero: Bool = true
      try
        while i < na do
          if _digits(i)? != 0 then
            all_zero = false
            break
          end
          i = i + 1
        end
      end
      if all_zero then Equal else Greater end
    elseif na < nb then
      var all_zero: Bool = true
      try
        while i < nb do
          if that._digits(i)? != 0 then
            all_zero = false
            break
          end
          i = i + 1
        end
      end
      if all_zero then Equal else Less end
    else
      Equal
    end


  fun _add_mag(that: MPFRep, sgn: Bool): MPFRep =>
    """
    Add the magnitudes `|this|` and `|that|` and return the sum with the
    given sign `sgn`. Both operands must be finite. The result exponent
    equals `max(this._exponent, that._exponent)`, incremented by one when a
    carry propagates out of the most-significant column.

    Operands are aligned on their most-significant digit before addition.
    The result precision equals `max(|this|_digits, shift + |that|_digits)`
    where `shift` is the exponent difference, so no precision is lost.

    The result is returned at full working precision (no output rounding).
    Rounding is applied by the caller (`MPFContext`).
    """
    // Align both operands: put the one with larger exponent first.
    (let ea, let eb, let ad, let bd, let na, let nb) =
      if _exponent >= that._exponent then
        (_exponent, that._exponent, _digits, that._digits, _size(), that._size())
      else
        (that._exponent, _exponent, that._digits, _digits, that._size(), _size())
      end
    let shift: USize = (ea - eb).usize()

    // When b falls entirely below a's precision, the sum equals a.
    if shift >= na then
      return MPFRep._create(sgn, false, false, ea, ad)
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
      MPFRep._create(sgn, false, false, ea, d)
    else
      // Carry out: one extra leading digit, exponent becomes ea + 1.
      MPFRep._create(sgn, false, false, ea + 1, raw)
    end


  fun _sub_mag(that: MPFRep, sgn: Bool): MPFRep =>
    """
    Subtract the magnitude `|that|` from `|this|`, where `|this| ≥ |that|`
    (enforced by the caller via `_cmp_mag`). Returns the difference with sign
    `sgn`. Both operands must be finite. The result is normalised: leading
    zero bytes are stripped and the exponent adjusted accordingly.

    The result is returned at full working precision (no output rounding).
    Rounding is applied by the caller (`MPFContext`).
    """
    let ea: I64 = _exponent
    let eb: I64 = that._exponent
    let shift: USize = (ea - eb).usize()
    let na: USize = _size()
    let nb: USize = that._size()
    let result_size: USize = na.max(shift + nb)
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
      return MPFRep.create(na)
    end

    let new_size: USize = raw_size - leading
    let new_exp: I64 = ea - leading.i64()
    let d: Array[U8] val = recover
      let a2 = Array[U8].init(0, new_size)
      raw.copy_to(a2, leading, 0, new_size)
      a2
    end
    MPFRep._create(sgn, false, false, new_exp, d)


  //- Short arithmetic (structural, scalar operand) ---------------------------

  fun _short_add(b: U8): MPFRep =>
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
        Fail(Format("[MPFRep._short_add] index out of bounds at i={}", err_i))
      end
      res
    end
    MPFRep._create(_sign, _nan, _inf, _exponent, d)


  fun _short_mul(b: U8): MPFRep =>
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
        Fail(Format("[MPFRep._short_mul] index out of bounds at i={}", err_i))
      end
      res
    end
    MPFRep._create(_sign, _nan, _inf, _exponent, d)


  fun _short_div(b: U8): (MPFRep, U8) =>
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
          (let quotient, let r2) =
            (remain.shl(8) + _digits(i)?.u16()).divrem(b.u16())
          remain = r2
          res.update(i, _lowb(quotient))?
          i = i + 1
        end
      else
        Fail(Format("[MPFRep._short_div] index out of bounds at i={}", err_i))
      end
      res
    end
    (MPFRep._create(_sign, _nan, _inf, _exponent, d), remain.u8())


  fun digit_shl(n: USize = 1): MPFRep =>
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
        Fail(Format("[MPFRep.digit_shl] index out of bounds at i={}", i))
      end
      a
    end
    MPFRep._create(_sign, _nan, _inf, _exponent, d)


  //- Hashing -----------------------------------------------------------------

  fun hash(): USize =>
    """
    Calculate a hash of this `MPFRep` based on the full mantissa, exponent and
    sign. Two representations that compare equal will have the same hash.
    """
    hash64().usize()


  fun hash64(): U64 =>
    """
    Calculate a 64-bit hash of this `MPFRep` based on the full mantissa,
    exponent and sign. Two representations that compare equal will have the
    same hash.
    """
    // FNV-1a hash: hash = offset_basis; for each byte b: hash = (hash XOR b) * prime.
    // The offset basis seeds the hash so an empty input doesn't hash to 0 and short
    // inputs don't cluster near 0. The constant 14695981039346656037 is the standard
    // 64-bit FNV offset basis, chosen to minimise collisions across the full U64 space.
    var h: U64 = 14695981039346656037
    let prime: U64 = 1099511628211

    // Mix sign, nan, inf, exponent into the hash.
    h = (h xor (if _sign then 1 else 0 end)) * prime
    h = (h xor (if _nan  then 2 else 0 end)) * prime
    h = (h xor (if _inf  then 4 else 0 end)) * prime
    h = (h xor _exponent.u64()) * prime

    // Mix all mantissa bytes.
    for d in _digits.values() do
      h = (h xor d.u64()) * prime
    end
    h


  //- Conversions: to floating-point ------------------------------------------

  fun f64(): F64 =>
    """
    Convert the current `MPFRep` to `F64`. Overflows are converted to ±∞;
    underflows to 0.

    This always holds: `MPFRep.from_f64(f).f64() == f` for any finite `F64`.

    The algorithm limits mantissa accumulation to 16 digits (128 bits, covering
    F64's 53-bit mantissa) to prevent intermediate overflow. The exponent
    scaling is split into two halves when outside the safe range [−125, 125]
    to avoid intermediate underflow for near-subnormal values.
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

    var result: F64 = 0.0
    let size = _size()
    let p = size.min(16)
    let base_f64: F64 = 256.0
    let inv_base = 1.0 / base_f64

    try
      var i: USize = p
      while i > 0 do
        i = i - 1
        result = (result + _digits(i)?.f64()) * inv_base
      end
    end

    let final_value =
      if (_exponent >= -125) and (_exponent <= 125) then
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
    Convert the current `MPFRep` to `F32`. Overflows are converted to ±∞;
    underflows to 0.

    This always holds: `MPFRep.from_f32(f).f32() == f` for any finite `F32`.

    The algorithm limits mantissa accumulation to 16 digits to prevent
    intermediate overflow. The exponent scaling is split when outside
    [−125, 125] to avoid intermediate underflow.
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

    var result: F64 = 0.0
    let size = _size()
    let p = size.min(16)
    let base_f64: F64 = 256.0
    let inv_base = 1.0 / base_f64

    try
      var i: USize = p
      while i > 0 do
        i = i - 1
        result = (result + _digits(i)?.f64()) * inv_base
      end
    end

    let final_value_f64 =
      if (_exponent >= -125) and (_exponent <= 125) then
        result * base_f64.pow(_exponent.f64())
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


  //- Conversions: to integer (safe, saturating) ------------------------------

  fun i128(): I128 =>
    """
    Convert `this` to an `I128` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `I128.max_value()` or
    `I128.min_value()`.

    - NaN → 0
    - −∞ → `I128.min_value()`
    - +∞ → `I128.max_value()`
    - `]−1, +1[` → 0
    - Other values are saturated to the `I128` range.
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

    // I128 range is [−2^127, 2^127 − 1]. 2^127 = 128 × 256^15.
    if (_exponent > 16) or
       ((_exponent == 16) and (try _digits(0)? >= 128 else false end))
    then
      if _sign then
        return I128.min_value()
      else
        return I128.max_value()
      end
    end

    var res: U128 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << U128(8)) or _digits(i)?.u128()
      end
    end

    if _exponent.usize() > n then
      res = res << (8 * (_exponent.usize() - n)).u128()
    end

    if _sign then
      (-res).i128()
    else
      res.i128()
    end


  fun i64(): I64 =>
    """
    Convert `this` to an `I64` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `I64.max_value()` or
    `I64.min_value()`.

    - NaN → 0
    - −∞ → `I64.min_value()`
    - +∞ → `I64.max_value()`
    - `]−1, +1[` → 0
    - Other values are saturated to the `I64` range.
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

    // I64 range is [−2^63, 2^63 − 1]. 2^63 = 128 × 256^7.
    if (_exponent > 8) or
       ((_exponent == 8) and (try _digits(0)? >= 128 else false end))
    then
      if _sign then
        return I64.min_value()
      else
        return I64.max_value()
      end
    end

    var res: U64 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << U64(8)) or _digits(i)?.u64()
      end
    end

    if _exponent.usize() > n then
      res = res << (8 * (_exponent.usize() - n)).u64()
    end

    if _sign then
      (-res).i64()
    else
      res.i64()
    end


  fun i32(): I32 =>
    """
    Convert `this` to an `I32` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `I32.max_value()` or
    `I32.min_value()`.

    - NaN → 0
    - −∞ → `I32.min_value()`
    - +∞ → `I32.max_value()`
    - `]−1, +1[` → 0
    - Other values are saturated to the `I32` range.
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

    // I32 range is [−2^31, 2^31 − 1]. 2^31 = 128 × 256^3.
    if (_exponent > 4) or
       ((_exponent == 4) and (try _digits(0)? >= 128 else false end))
    then
      if _sign then
        return I32.min_value()
      else
        return I32.max_value()
      end
    end

    var res: U32 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << U32(8)) or _digits(i)?.u32()
      end
    end

    if _exponent.usize() > n then
      res = res << (8 * (_exponent.usize() - n)).u32()
    end

    if _sign then
      (-res).i32()
    else
      res.i32()
    end


  fun i16(): I16 =>
    """
    Convert `this` to an `I16` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `I16.max_value()` or
    `I16.min_value()`.

    - NaN → 0
    - −∞ → `I16.min_value()`
    - +∞ → `I16.max_value()`
    - `]−1, +1[` → 0
    - Other values are saturated to the `I16` range.
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

    // I16 range is [−2^15, 2^15 − 1]. 2^15 = 128 × 256^1.
    if (_exponent > 2) or
       ((_exponent == 2) and (try _digits(0)? >= 128 else false end))
    then
      if _sign then
        return I16.min_value()
      else
        return I16.max_value()
      end
    end

    var res: U16 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << U16(8)) or _digits(i)?.u16()
      end
    end

    if _exponent.usize() > n then
      res = res << (8 * (_exponent.usize() - n)).u16()
    end

    if _sign then
      (-res).i16()
    else
      res.i16()
    end


  fun i8(): I8 =>
    """
    Convert `this` to an `I8` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `I8.max_value()` or
    `I8.min_value()`.

    - NaN → 0
    - −∞ → `I8.min_value()`
    - +∞ → `I8.max_value()`
    - `]−1, +1[` → 0
    - Other values are saturated to the `I8` range.
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

    // I8 range is [−2^7, 2^7 − 1]. 2^7 = 128.
    if (_exponent > 1) or
       ((_exponent == 1) and (try _digits(0)? >= 128 else false end))
    then
      if _sign then
        return I8.min_value()
      else
        return I8.max_value()
      end
    end

    var res = try _digits(0)?.u8() else 0 end

    if _exponent > 1 then
      res = res << (U8(8) * (_exponent.u8() - 1))
    end

    if _sign then
      (-res).i8()
    else
      res.i8()
    end


  fun ilong(): ILong =>
    """
    Convert `this` to an `ILong` value. On LLP64 and ILP32 platforms that is
    a conversion to `I32`; on all other platforms it is a conversion to `I64`.

    - NaN → 0
    - −∞ → `ILong.min_value()`
    - +∞ → `ILong.max_value()`
    - `]−1, +1[` → 0
    - Other values are saturated to the `ILong` range.
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
    Convert `this` to an `ISize` value. On ILP32 platforms that is a
    conversion to `I32`; on all other platforms it is a conversion to `I64`.

    - NaN → 0
    - −∞ → `ISize.min_value()`
    - +∞ → `ISize.max_value()`
    - `]−1, +1[` → 0
    - Other values are saturated to the `ISize` range.
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


  //- Conversions: to integer (unsafe) ----------------------------------------

  fun i128_unsafe(): I128 =>
    """
    Convert `this` to an `I128` value. Fractional digits are truncated.
    In case of overflow or special value the result is undefined.
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
    Convert `this` to an `I64` value. Fractional digits are truncated.
    In case of overflow or special value the result is undefined.
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
    Convert `this` to an `I32` value. Fractional digits are truncated.
    In case of overflow or special value the result is undefined.
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
    Convert `this` to an `I16` value. Fractional digits are truncated.
    In case of overflow or special value the result is undefined.
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
    Convert `this` to an `I8` value. Fractional digits are truncated.
    In case of overflow or special value the result is undefined.
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
    Convert `this` to an `ILong` value without checking for overflow or
    special values. On LLP64 and ILP32 platforms that is a conversion to
    `I32`; on all other platforms it is a conversion to `I64`.
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
    Convert `this` to an `ISize` value without checking for overflow or
    special values. On ILP32 platforms that is a conversion to `I32`;
    on all other platforms it is a conversion to `I64`.
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


  //- Conversions: to unsigned integer ----------------------------------------

  fun u128(): U128 =>
    """
    Convert `this` to a `U128` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `U128.max_value()` or
    `U128.min_value()` (0).

    - NaN → 0
    - −∞ → `U128.min_value()` (0)
    - +∞ → `U128.max_value()`
    - `]−1, +1[` → 0
    - negative values → 0
    - Other values are saturated to the `U128` range.
    """
    if _nan then
      return 0
    end
    if _inf then
      if _sign then
        return U128.min_value()
      else
        return U128.max_value()
      end
    end
    if is_zero() or (_exponent <= 0) then
      return 0
    end
    if _sign then
      return U128.min_value()
    end

    // U128 range is [0, 2^128 − 1]. 2^128 = 256^16.
    if _exponent > 16 then
      return U128.max_value()
    end

    var res: U128 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << U128(8)) or _digits(i)?.u128()
      end
    end

    if _exponent.usize() > n then
      res = res << (8 * (_exponent.usize() - n)).u128()
    end

    res


  fun u64(): U64 =>
    """
    Convert `this` to a `U64` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `U64.max_value()` or
    `U64.min_value()` (0).

    - NaN → 0
    - −∞ → `U64.min_value()` (0)
    - +∞ → `U64.max_value()`
    - `]−1, +1[` → 0
    - negative values → 0
    - Other values are saturated to the `U64` range.
    """
    if _nan then
      return 0
    end
    if _inf then
      if _sign then
        return U64.min_value()
      else
        return U64.max_value()
      end
    end
    if is_zero() or (_exponent <= 0) then
      return 0
    end
    if _sign then
      return U64.min_value()
    end

    // U64 range is [0, 2^64 − 1]. 2^64 = 256^8.
    if _exponent > 8 then
      return U64.max_value()
    end

    var res: U64 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << U64(8)) or _digits(i)?.u64()
      end
    end

    if _exponent.usize() > n then
      res = res << (8 * (_exponent.usize() - n)).u64()
    end

    res


  fun u32(): U32 =>
    """
    Convert `this` to a `U32` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `U32.max_value()` or
    `U32.min_value()` (0).

    - NaN → 0
    - −∞ → `U32.min_value()` (0)
    - +∞ → `U32.max_value()`
    - `]−1, +1[` → 0
    - negative values → 0
    - Other values are saturated to the `U32` range.
    """
    if _nan then
      return 0
    end
    if _inf then
      if _sign then
        return U32.min_value()
      else
        return U32.max_value()
      end
    end
    if is_zero() or (_exponent <= 0) then
      return 0
    end
    if _sign then
      return U32.min_value()
    end

    // U32 range is [0, 2^32 − 1]. 2^32 = 256^4.
    if _exponent > 4 then
      return U32.max_value()
    end

    var res: U32 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << U32(8)) or _digits(i)?.u32()
      end
    end

    if _exponent.usize() > n then
      res = res << (8 * (_exponent.usize() - n)).u32()
    end

    res


  fun u16(): U16 =>
    """
    Convert `this` to a `U16` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `U16.max_value()` or
    `U16.min_value()` (0).

    - NaN → 0
    - −∞ → `U16.min_value()` (0)
    - +∞ → `U16.max_value()`
    - `]−1, +1[` → 0
    - negative values → 0
    - Other values are saturated to the `U16` range.
    """
    if _nan then
      return 0
    end
    if _inf then
      if _sign then
        return U16.min_value()
      else
        return U16.max_value()
      end
    end
    if is_zero() or (_exponent <= 0) then
      return 0
    end
    if _sign then
      return U16.min_value()
    end

    // U16 range is [0, 2^16 − 1]. 2^16 = 256^2.
    if _exponent > 2 then
      return U16.max_value()
    end

    var res: U16 = 0
    let n: USize = _exponent.usize().min(_digits.size())
    try
      for i in Range(0, n) do
        res = (res << U16(8)) or _digits(i)?.u16()
      end
    end

    if _exponent.usize() > n then
      res = res << (8 * (_exponent.usize() - n)).u16()
    end

    res


  fun u8(): U8 =>
    """
    Convert `this` to a `U8` value. Fractional digits are truncated toward
    zero. In case of overflow the value is saturated to `U8.max_value()` or
    `U8.min_value()` (0).

    - NaN → 0
    - −∞ → `U8.min_value()` (0)
    - +∞ → `U8.max_value()`
    - `]−1, +1[` → 0
    - negative values → 0
    - Other values are saturated to the `U8` range.
    """
    if _nan then
      return 0
    end
    if _inf then
      if _sign then
        return U8.min_value()
      else
        return U8.max_value()
      end
    end
    if is_zero() or (_exponent <= 0) then
      return 0
    end
    if _sign then
      return U8.min_value()
    end

    // U8 range is [0, 255]. 2^8 = 256.
    if _exponent > 1 then
      return U8.max_value()
    end

    try _digits(0)? else 0 end


  fun ulong(): ULong =>
    """
    Convert `this` to a `ULong` value. On LLP64 and ILP32 platforms that is
    a conversion to `U32`; on all other platforms it is a conversion to `U64`.

    - NaN → 0
    - −∞ → `ULong.min_value()` (0)
    - +∞ → `ULong.max_value()`
    - negative values → 0
    - Other values are saturated to the `ULong` range.
    """
    ifdef lp64 then
      u64().ulong()
    elseif llp64 then
      u32().ulong()
    elseif ilp32 then
      u32().ulong()
    else
      u64().ulong()
    end


  fun usize(): USize =>
    """
    Convert `this` to a `USize` value. On ILP32 platforms that is a
    conversion to `U32`; on all other platforms it is a conversion to `U64`.

    - NaN → 0
    - −∞ → `USize.min_value()` (0)
    - +∞ → `USize.max_value()`
    - negative values → 0
    - Other values are saturated to the `USize` range.
    """
    ifdef lp64 then
      u64().usize()
    elseif llp64 then
      u64().usize()
    elseif ilp32 then
      u32().usize()
    else
      u64().usize()
    end


  //- String output -----------------------------------------------------------

  fun _digits_as_mpint(): MPInt =>
    """
    Convert `_digits` (the big-endian base-256 significand of this `MPFRep`,
    treated as a non-negative integer) into an `MPInt`.

    The result is always non-negative regardless of `_sign`. Uses
    `MPInt.from_bytes_be` so no internal `MPInt` detail is exposed.
    """
    MPInt.from_bytes_be(false, _digits)


  fun exact_string(base: U8 = 10): (String iso^, I64, Bool) =>
    """
    Return `(mantissa, exponent, inexact)` where:

      `|this| = 0.mantissa_digits × 10^exponent`

    The mantissa string begins with `"-"` when `this` is negative, followed
    by decimal digit characters (no decimal point). The number of digits is
    approximately `⌈_size() × log₁₀(256)⌉ + 2`; trailing zeros are NOT
    stripped. The `exponent` is the number of decimal digits to the left of
    the decimal point (so `exponent = 1` means the value is in `[0.1, 1)`
    after dividing by `10^(exponent − 1)`).

    Special values:
    - NaN → `("nan", 0, false)`, ±∞ → `("±inf", 0, false)`, ±0 → `("0", 0, false)`

    Only base 10 is currently implemented; other bases return `("", 0, true)`.
    The `inexact` flag is always `false` (rounding is a future TODO).

    Algorithm (exact, no Newton approximation):
    Let `N` be `_digits` interpreted as a base-256 integer and
    `k = _exponent − _size()`. Then `|this| = N × 256^k`.

    - Integer case (`k ≥ 0`): `|this| = N × 2^{8k}`.
      Compute `N << (8k)` as an `MPInt` and call `.string()`.
    - Fractional case (`k < 0`): `|this| = N × 5^b / 10^b` where `b = 8|k|`.
      Compute `(N × 5^b).string()` via `MPInt` and set
      `dec_exp = len(str) − b`.

    Both paths are always exact — no scaling error, no Newton undershooting.
    """
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
    let n_dec: USize = ((prec.f64() * 2.41).usize() + 2).max(1)
    let sign_prefix: String = if _sign then "-" else "" end

    let n_int: MPInt = _digits_as_mpint()

    let k_bytes: I64 = _exponent - prec.i64()

    if (k_bytes >= 0) and (k_bytes < 150) then
      let b: USize = k_bytes.usize() * 8
      let shifted: MPInt =
        if b > 0 then
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
      let b: USize = ((-k_bytes) * 8).usize()
      let five_pow: MPInt =
        MPInt.from[ILong](5).pow(MPInt.from[ILong](b.ilong()))
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
      // Extreme values: use approximate path via F64 to avoid OOM in MPInt/String.
      ifdef debug then
        Debug(
          "[MPFRep.exact_string] Extreme value (k_bytes=" + k_bytes.string() +
          "). Using approximate F64 path.")
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

        let exp_val = e_str.i64()? +
          (if dot_seen then dot_pos.i64() else m_str.size().i64() end)
        (consume m_raw, exp_val, true)
      else
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
        let exp_val =
          if dot_seen then dot_pos.i64() else s_f_val.size().i64() end
        (consume m_raw, exp_val, true)
      end
    end


  fun string(): String iso^ =>
    """
    Return a decimal string representation of this `MPFRep`, formatted in the
    style of C's `printf` `%g` with `⌈_size() × log₁₀(256)⌉ + 2` significant
    digits. Trailing zeros are stripped; at least one digit after the decimal
    point is always shown.

    Special values: `"nan"`, `"+inf"`, `"-inf"`, `"0.0"`, `"-0.0"`.

    Scientific notation (`d.dddde±N`) is used when the decimal exponent is
    less than `−4` or greater than or equal to the number of significant
    digits; fixed notation (`ddd.ddd`) is used otherwise.
    """
    if _nan then
      return "nan".clone()
    end
    if _inf then
      return (if _sign then "-inf" else "+inf" end).clone()
    end
    if is_zero() then
      return (if _sign then "-0.0" else "0.0" end).clone()
    end

    (let raw_mant, let dec_exp, _) = exact_string(10)
    let raw: String = consume raw_mant

    var dig_start: USize = 0
    let sgn: String =
      if (raw.size() > 0) and (try raw(0)? == '-' else false end) then
        dig_start = 1
        "-"
      else
        ""
      end

    let all_digits: String = raw.trim(dig_start)
    var sig_end: USize = all_digits.size()
    while (sig_end > 1) and
          (try all_digits(sig_end - 1)? == '0' else false end)
    do
      sig_end = sig_end - 1
    end
    let digits: String = all_digits.trim(0, sig_end)
    let nd: USize = digits.size()

    let use_sci: Bool =
      (dec_exp < -4) or (dec_exp >= all_digits.size().i64())

    recover
      let s = String.create(sgn.size() + nd + 8)
      s.append(sgn)
      if use_sci then
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
          s.append("0.")
          var zi: ISize = 0
          while zi < -e do
            s.push('0')
            zi = zi + 1
          end
          s.append(digits)
        elseif e >= ni then
          s.append(digits)
          var zi: ISize = e - ni
          while zi > 0 do
            s.push('0')
            zi = zi - 1
          end
          s.append(".0")
        else
          s.append(digits.trim(0, e.usize()))
          s.append(".")
          s.append(digits.trim(e.usize()))
        end
      end
      s
    end
