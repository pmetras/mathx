// Multi-precision floating point: pure representation layer.

use "../assertx"
use "../formatx"

use "collections"
use "debug"


class ref MPFRep is (Equatable[MPFRep] & Formattable & Stringable)
  """
  The pure representation of an arbitrary-precision floating-point number.

  `MPFRep` stores the sign, special-value flags, base-256 exponent and mantissa
  digits of a number. The default capability is `ref` (mutable), which allows
  internal arithmetic loops in `_MPFAlgo` to update intermediate results in
  place without allocating a new instance at every step. Callers that need a
  sendable, globally-immutable value consume the `iso^` returned by the
  constructors into a `val` viewpoint before storing it (e.g. in `MPFloat._rep`).

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

  ## Rounding guarantees (via `_MPFAlgo`)

  | Operation                        | Guarantee          | Method                    |
  | -------------------------------- | ------------------ | ------------------------- |
  | `add`, `sub`, `mul`              | Correct (≤ 0.5 ULP) | IEEE 754 §5.4: exact intermediate + one round |
  | `div`, `inv`, `sqrt`             | Correct (≤ 0.5 ULP) | Newton + one refinement step                  |
  | `ln`, `log2`, `log10`            | Correct (≤ 0.5 ULP) | Ziv's iteration                               |
  | `exp`, `exp2`                    | Correct (≤ 0.5 ULP) | Ziv's iteration                               |
  | `powi`, `pow`                    | Correct (≤ 0.5 ULP) | Ziv's iteration                               |
  | `logb`                           | Faithful (≤ 1 ULP)  | `ln/ln` — two roundings                       |
  | `sin`, `cos`, `tan`              | Faithful (≤ 1 ULP)  | Taylor + arg. reduction                       |
  | `sinh`, `cosh`, `tanh`           | Faithful (≤ 1 ULP)  | Taylor series                                 |
  | `sech`, `csch`, `coth`           | Faithful (≤ 1 ULP)  | composed from above                           |
  | `pi`, `pi_bbp`                   | Faithful (≤ 1 ULP)  | iterative formula                             |

  **Correctly rounded** means the result equals the infinite-precision value
  rounded under the requested `RoundingMode` — identical to what MPFR
  guarantees. **Faithfully rounded** means the error is at most 1 ULP in
  either direction.
  """

  var _sign: Bool
    """
    Sign bit: `true` when the value is strictly negative or is −0.
    Ignored when `_nan` is `true`.
    """

  var _nan: Bool
    """
    `true` when the number is Not-a-Number. When set, `_inf`, `_sign` and
    `_digits` are irrelevant.
    """

  var _inf: Bool
    """
    `true` when the number is ±∞. Sign is given by `_sign`.
    Ignored when `_nan` is `true`.
    """

  var _exponent: I64
    """
    Base-256 exponent `e` such that `|value| = 0.d₀d₁… × 256^e`.
    For the common arithmetic operations `e = 1` is used internally
    (the first digit `d₀` is the integer part).
    """

  var _digits: Array[U8] val
    """
    Big-endian base-256 mantissa digits (index 0 = most significant).
    For a finite non-zero number `_digits(0) ≠ 0` (normalised).
    For zero the array may be all-zeros or empty.
    Immutable `val` so the array can safely be captured inside `recover`
    blocks and passed to sendable-parameter constructors. `_update` swaps
    the pointer (no copy) to reuse the `MPFRep` struct across Newton/Taylor
    iterations without allocating a new struct on every step.
    """


  new iso _create(
    sgn: Bool,
    nan: Bool,
    inf: Bool,
    expn: I64,
    digits: Array[U8] val)
  =>
    """
    Private canonical constructor. All constructors and internal operations
    produce their result through this entry point.

    Returns an `iso^` so that callers may consume it to `val` (for storage in
    `MPFloat._rep`) or to `ref` (for mutable accumulator reuse in `_MPFAlgo`).
    """
    _sign = sgn
    _nan = nan
    _inf = inf
    _exponent = expn
    _digits = digits


  new iso create(p_bytes: USize = 14) =>
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


  new iso nan_val() =>
    """
    Create a Not-a-Number representation.
    """
    _sign = false
    _nan = true
    _inf = false
    _exponent = 0
    _digits = Array[U8].create()


  new iso inf_val(positive: Bool = true) =>
    """
    Create an infinite representation. Pass `positive = false` for −∞.
    """
    _sign = not positive
    _nan = false
    _inf = true
    _exponent = 0
    _digits = Array[U8].create()


  new iso from[A: ((Number | MPInt | MPFRep) & Real[A] val)](n: A,
                                                             p_bytes: USize = 14) =>
    """
    Create a new `MPFRep` from any numeric value `n` with `p_bytes` base-256
    mantissa bytes (default 14 ≈ 112 bits ≈ 33 decimal digits).

    Accepts:
    - `MPFRep` — precision-change copy (truncates, no rounding)
    - `MPInt`  — exact integer conversion (truncated to `p_bytes`)
    - `F64`    — preserves special values (NaN, ±∞, ±0)
    - `F32`    — preserves special values (NaN, ±∞, ±0)
    - `ULong`, `USize`, `U8`…`U128` — exact unsigned integer
    - `ILong`, `ISize`, `I8`…`I128` — exact signed integer (via MPInt)

    The conversion bodies live in `_MPFAlgo._from_*`; this constructor
    delegates to them and unpacks the returned `val` into `this`.
    """
    let rep: MPFRep val =
      iftype A <: MPFRep then
        let r: MPFRep val = n
        r._trunc(p_bytes)
      elseif A <: MPInt then
        let m: MPInt = n
        _MPFAlgo._from_mpint(m, p_bytes)
      elseif A <: F64 then
        _MPFAlgo._from_f64(n.f64(), p_bytes)
      elseif A <: F32 then
        _MPFAlgo._from_f32(n.f32(), p_bytes)
      elseif A <: (I8 | I16 | I32 | I64 | I128 | ILong | ISize) then
        _MPFAlgo._from_mpint(MPInt.from[I128](n.i128()), p_bytes)
      elseif A <: U128 then
        _MPFAlgo._from_mpint(MPInt.from[U128](n.u128()), p_bytes)
      elseif A <: (U8 | U16 | U32 | U64 | ULong | USize) then
        _MPFAlgo._from_ulong(n.ulong(), p_bytes)
      else
        Debug("[MPFRep.from] Unknown type. Constructor failed! Adapt MPFRep code.")
        MPFRep._create(false, false, false, 0, Array[U8].init(0, p_bytes))
      end
    _sign = rep._sign
    _nan = rep._nan
    _inf = rep._inf
    _exponent = rep._exponent
    _digits = rep._digits


  new \do_not_use\ iso min_normalized(p_bytes: USize = 14) =>
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


  new iso epsilon(p_bytes: USize = 14) =>
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


  new iso min_value() =>
    """
    The minimum representable value is −∞ (`-inf`).
    """
    _sign = true
    _nan = false
    _inf = true
    _exponent = 0
    _digits = Array[U8].create()


  new iso max_value() =>
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


  fun _inc_first(): MPFRep val =>
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


  fun _trunc(n: USize): MPFRep val =>
    """
    Return a new `MPFRep` containing only the `n` most-significant base-256
    digits of `this`. Used internally to bound intermediate results in Newton
    iterations so that successive multiplications do not grow the array without
    limit.

    This is a pure truncation (no rounding). Output rounding is the
    responsibility of `MPFContext._round_to`.

    Uses `Array.trim` on the `val` field — shares the backing store, no copy.
    """
    MPFRep._create(_sign, _nan, _inf, _exponent, _digits.trim(0, _size().min(n)))


  fun _neg_comp(): MPFRep val =>
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


  fun _trunc_frac(): MPFRep val =>
    """
    Return `this` with the fractional part removed (truncated toward zero).
    Internal alias equivalent to `trunc()` on `MPFloat`. Called by division
    and floored-division routines to obtain the integer part.
    """
    if not is_finite() then
      // Reuse the val reference — no copy needed.
      return MPFRep._create(_sign, _nan, _inf, _exponent, _digits)
    end
    if is_zero() then
      return MPFRep._create(_sign, false, false, 0, _digits)
    end
    if _exponent <= 0 then
      return MPFRep.create(_size())
    end
    let e: USize = _exponent.usize().min(_size())
    // Zero the fractional bytes in a new array (trim can't zero-fill the tail).
    let d: Array[U8] val = recover
      let a = Array[U8].init(0, _size())
      _digits.copy_to(a, 0, 0, e)
      a
    end
    MPFRep._create(_sign, false, false, _exponent, d)


  fun _clone(): MPFRep iso^ =>
    """
    Return a fresh copy of this `MPFRep` as an owned `iso^`.
    Useful when a loop needs a mutable `ref` accumulator seeded from a `box`
    or `val` source: `var acc: MPFRep ref = some_box._clone()`.
    """
    MPFRep._create(_sign, _nan, _inf, _exponent, _digits)


  fun ref _update(that: MPFRep box) =>
    """
    Replace every field of this `MPFRep` with the values from `that` in place,
    reusing the existing struct instead of allocating a new one.
    Called in Newton and Taylor loops to avoid per-iteration struct allocation.

    The digit array pointer is swapped (no byte copy) since `_digits` is `val`:
    `that._digits` (a fresh `val` produced by each arithmetic step) is installed
    directly, saving the cost of allocating a new `MPFRep` struct on every step.
    """
    _sign     = that.sign_bit()
    _nan      = that.is_nan()
    _inf      = that.is_infinite()
    _exponent = that.exponent()
    _digits   = that.raw_digits()


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
    Return the internal base-256 mantissa as a big-endian `Array[U8] val`.

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


  fun eq(that: MPFRep box): Bool =>
    """
    Value equality: `true` iff both values represent the same number.

    IEEE 754 rules apply:
    - NaN ≠ NaN (always `false` when either operand is NaN).
    - +0 = −0.
    - Two infinities are equal iff they have the same sign.
    - Finite values are equal iff sign, exponent, and all digits agree
      (after stripping trailing zero bytes).
    """
    if _nan or that.is_nan() then return false end
    if _inf or that.is_infinite() then
      return _inf and that.is_infinite() and (_sign == that.sign_bit())
    end
    // +0 = -0: ignore sign for zero values.
    let a_zero = is_zero()
    let b_zero = that.is_zero()
    if a_zero or b_zero then return a_zero and b_zero end
    if _sign != that.sign_bit() then return false end
    if _exponent != that.exponent() then return false end
    // Compare digits, ignoring trailing zeros on either side.
    let ad = _digits
    let bd = that.raw_digits()
    let n: USize = ad.size().max(bd.size())
    try
      for i in Range(0, n) do
        let a_byte: U8 = if i < ad.size() then ad(i)? else 0 end
        let b_byte: U8 = if i < bd.size() then bd(i)? else 0 end
        if a_byte != b_byte then return false end
      end
    end
    true

  fun ne(that: MPFRep box): Bool => not eq(that)


  //- Magnitude arithmetic (structural, no rounding) --------------------------

  fun _cmp_mag(that: MPFRep box): Compare =>
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


  fun _add_mag(that: MPFRep box, sgn: Bool, w: USize = USize.max_value()): MPFRep val =>
    """
    Add the magnitudes `|this|` and `|that|` and return the sum with sign
    `sgn`. Both operands must be finite.

    `w` caps the output to at most `w` most-significant bytes (default:
    unlimited). Pass the working precision to skip a separate `_trunc` call.

    The result exponent equals `max(this._exponent, that._exponent)`,
    incremented by one when a carry propagates out of the MSB column.
    """
    // Align both operands: put the one with larger exponent first.
    (let ea, let eb, let ad, let bd, let na, let nb) =
      if _exponent >= that._exponent then
        (_exponent, that._exponent, _digits, that._digits, _size(), that._size())
      else
        (that._exponent, _exponent, that._digits, _digits, that._size(), _size())
      end
    let shift: USize = (ea - eb).usize()

    // When b falls entirely below a's precision, the sum equals a (capped at w).
    if shift >= na then
      return MPFRep._create(sgn, false, false, ea, ad.trim(0, na.min(w)))
    end

    // Run the addition loop over all prec columns (carry propagates LSB→MSB).
    // Write the carry out of the MSB to an outer var (value type — Pony allows
    // writing outer var bindings of machine-word types from inside recover).
    // Trim to w bytes in-place before leaving the recover block (zero extra
    // allocation for the common no-carry case).
    let prec: USize = na.max(shift + nb)
    var carry_out: U8 = 0
    let d: Array[U8] val = recover
      let res = Array[U8].init(0, prec)
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
          res.update(col, sum)?
        until col == 0 end
        carry_out = _lowb(carry)
      end
      // Trim in-place to w bytes: adjusts internal length, no copy.
      res.trim_in_place(0, prec.min(w))
      res
    end

    if carry_out == 0 then
      // No carry: result is exactly d (already trimmed), exponent stays ea.
      MPFRep._create(sgn, false, false, ea, d)
    else
      // Carry out (rare): prepend the carry byte, exponent becomes ea + 1.
      let d2: Array[U8] val = recover
        let out = Array[U8].init(0, d.size() + 1)
        try out.update(0, carry_out)? end
        d.copy_to(out, 0, 1, d.size())
        out
      end
      MPFRep._create(sgn, false, false, ea + 1, d2)
    end


  fun _sub_mag(that: MPFRep box, sgn: Bool, w: USize = USize.max_value()): MPFRep val =>
    """
    Subtract the magnitude `|that|` from `|this|`, where `|this| ≥ |that|`
    (enforced by the caller via `_cmp_mag`). Returns the difference with sign
    `sgn`. Both operands must be finite. The result is normalised: leading
    zero bytes are stripped and the exponent adjusted accordingly.

    `w` caps the output to at most `w` most-significant bytes (default:
    unlimited). Pass the working precision to skip a separate `_trunc` call.
    """
    let ea: I64 = _exponent
    let eb: I64 = that._exponent
    let shift: USize = (ea - eb).usize()
    let na: USize = _size()
    let nb: USize = that._size()
    let result_size: USize = na.max(shift + nb)
    let max_b: U16 = U8.max_value().u16()
    // Compute the subtraction, find leading zeros, strip them, and trim to w
    // bytes in a single recover block — no separate allocation.
    var leading: USize = 0
    let d: Array[U8] val = recover
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
        // Find first non-zero byte inline.
        var lead: USize = 0
        while (lead < result_size) and (res(lead)? == 0) do
          lead = lead + 1
        end
        leading = lead
        if lead == result_size then
          // All zero — return empty array as sentinel; handled below.
          Array[U8].create()
        else
          let keep = (result_size - lead).min(w)
          // trim_in_place: adjust length to [lead, lead+keep), no copy.
          res.trim_in_place(lead, lead + keep)
          res
        end
      else
        Array[U8].create()
      end
    end

    if d.size() == 0 then
      return MPFRep.create(na)
    end
    let new_exp: I64 = ea - leading.i64()
    MPFRep._create(sgn, false, false, new_exp, d)


  //- Short arithmetic (structural, scalar operand) ---------------------------

  fun _short_add(b: U8): MPFRep val =>
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


  fun _short_mul(b: U8): MPFRep val =>
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


  fun _short_div(b: U8): (MPFRep val, U8) =>
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


  fun digit_shl(n: USize = 1): MPFRep val =>
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

    This always holds: `MPFRep.from[F64](f).f64() == f` for any finite `F64`.

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

    This always holds: `MPFRep.from[F32](f).f32() == f` for any finite `F32`.

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
    - NaN → `("nan", 0, false)`
    - ±∞ → `("±inf", 0, false)`
    - ±0 → `("0", 0, false)`

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

    **DoS counter-measure — exact vs approximate path:**

    Let `k = _exponent − _size()` and `gate = _size() + 150`. The exact MPInt
    paths above require computing either `N × 2^{8k}` (integer shift) or
    `N × 5^{8|k|}` (fractional scaling). The `5^{8|k|}` term grows rapidly:
    when `|k|` is large relative to `_size()` the intermediate `MPInt` can
    reach millions of digits, making the call both slow and memory-hungry for
    adversarially constructed values.

    When `|k| ≥ gate` (`k ≤ −gate` or `k ≥ gate`) the method falls back to an
    **approximate F64 path**:
    - `inexact` is always `true`.
    - The mantissa carries only ≈ 15–17 significant decimal digits (F64 precision),
      regardless of the stored `MPFRep` precision.
    - Values that overflow F64 (roughly `|value| > 1.8 × 10^308`) return the
      sentinel tuple `("1", 1000, true)` — do **not** interpret the dec_exp
      literally; it signals "too large for exact representation".
    - Values that underflow F64 (roughly `|value| < 5 × 10^-324`) return the
      sentinel tuple `("0", −1000, true)` — similarly a non-literal sentinel for
      "too small for exact representation".
    - A `Debug` line is emitted on every fallback invocation.

    The proportional gate keeps the exact-path output bounded to at most
    `(2×_size() + 150) × 2.41` decimal digits — feasible for any sane precision.
    A flat gate of 150 was too tight for high-precision values (e.g. π stored at
    437 bytes has `k = −436` but only produces a ~3500-digit string, well within
    the capabilities of `MPInt`).

    For a 14-byte (112-bit, ≈ F128) `MPFRep`, `gate = 164`, so the exact path
    covers values roughly in `(10^{−395}, 10^{+395})`.
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
    // Gate: allow exact MPInt path when |k_bytes| < prec + 150.
    // Flat 150 was too tight for high-precision values (e.g. π at 437 bytes has
    // k_bytes = −436).  The proportional gate keeps the output string bounded to
    // at most ~(2×prec + 150) × 2.41 decimal digits — feasible for any sane prec.
    let gate: I64 = prec.i64() + 150

    if (k_bytes >= 0) and (k_bytes < gate) then
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
        for i in Range(dec_exp.usize(), n_dec) do
          s.push('0')
        end
        s
      end
      (consume result, dec_exp, false)
    elseif (k_bytes < 0) and (k_bytes > -gate) then
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
        for i in Range(num_len, n_dec) do
          s.push('0')
        end
        s
      end
      (consume result, dec_exp, false)
    else
      // Extreme values: use approximate path via F64 to avoid OOM in MPInt/String.
      // This is a DoS counter-measure.
        Debug("[MPFRep.exact_string] Extreme value (k_bytes=" + k_bytes.string() +
              "). Using approximate F64 path.")

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


  fun format(spec: String = ""): String =>
    """
    Format this value according to a `FormatSpec`.

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
    let fspec = FormatSpec(spec)

    // ---- special values ----
    if _nan then
      let upper = (fspec.type_char == 'E') or (fspec.type_char == 'F')
        or (fspec.type_char == 'G')
      let s: String val = if upper then "NAN" else "nan" end
      return _fmt_pad(s, fspec)
    end
    if _inf then
      let upper = (fspec.type_char == 'E') or (fspec.type_char == 'F')
        or (fspec.type_char == 'G')
      let inf_s: String val = if upper then "INF" else "inf" end
      let sgn: String val =
        if _sign then "-"
        elseif fspec.sign is SignPlus then "+"
        elseif fspec.sign is SignSpace then " "
        else "" end
      return _fmt_pad(sgn + inf_s, fspec)
    end

    // ---- integer type codes: truncate then delegate ----
    let tc = fspec.type_char
    if (tc == 'd') or (tc == 'b') or (tc == 'o') or (tc == 'x') or (tc == 'X') then
      // Truncate toward zero: keep only bytes at positive exponent positions.
      let int_digits: Array[U8] val =
        if _exponent <= 0 then
          recover [U8(0)] end
        elseif _exponent.usize() >= _size() then
          _digits
        else
          _digits.trim(0, _exponent.usize())
        end
      // Build a temporary MPFRep with no fractional part and format via MPInt path.
      let int_rep = MPFRep._create(_sign, false, false,
        int_digits.size().i64(), int_digits)
      return int_rep._fmt_int(fspec)
    end

    // ---- floating-point formatting ----
    if is_zero() then
      let zero_s: String val = if _sign then "-0" else "0" end
      let zero_rep: String val = recover
        let s = String
        let neg: Bool = _sign and not (fspec.z)
        if neg then s.push('-')
        elseif fspec.sign is SignPlus then s.push('+')
        elseif fspec.sign is SignSpace then s.push(' ')
        end
        s.append("0")
        if (tc == 'f') or (tc == 'F') then
          let p: USize = match fspec.precision | let p2: USize => p2 else 1 end
          if p > 0 then
            s.push('.')
            for i in Range(0, p) do
              s.push('0')
            end
          elseif fspec.hash then
            s.push('.')
          end
        elseif (tc == 'e') or (tc == 'E') then
          let p: USize = match fspec.precision | let p2: USize => p2 else 1 end
          s.push('.')
          for i in Range(0, p) do
            s.push('0')
          end
          s.append(if (tc == 'E') then "E+00" else "e+00" end)
        elseif (tc == '%') then
          let p: USize = match fspec.precision | let p2: USize => p2 else 6 end
          s.push('.')
          for i in Range(0, p) do
            s.push('0')
          end
          s.push('%')
        else
          // Default (no type code) and 'g'/'G': show "0.0".
          s.append(".0")
        end
        s
      end
      return _fmt_num_pad("", zero_rep, fspec)
    end

    // Get the full exact decimal mantissa and exponent.
    (let mant_iso, let dec_exp, _) = exact_string(10)
    let mant_val: String val = consume mant_iso

    // Strip sign from mantissa; record separately.
    var dig_start: USize = 0
    let neg: Bool = (mant_val.size() > 0) and
      (try mant_val(0)? == '-' else false end)
    if neg then
      dig_start = 1
    end
    let all_digits: String val = mant_val.trim(dig_start)

    let sign_str: String val = if neg then
        "-"
      elseif fspec.sign is SignPlus then
        "+"
      elseif fspec.sign is SignSpace then
        " "
      else
        ""
      end

    let upper = (tc == 'E') or (tc == 'F') or (tc == 'G')

    // Determine output precision: explicit spec, or full stored precision.
    let full_prec: USize = all_digits.size()

    match tc
    | 'f' | 'F' =>
      let p: USize = match fspec.precision
        | let p2: USize => p2
        else
          full_prec
        end
      let digits: String val = _fmt_fixed(all_digits, dec_exp, p, fspec.hash)
      let grouped: String val = _fmt_group_fixed(digits, fspec.grouping)
      _fmt_num_pad(sign_str, grouped, fspec)

    | 'e' | 'E' =>
      let p: USize = match fspec.precision
        | let p2: USize => p2
        else
          full_prec - 1
        end
      let digits: String val = _fmt_sci(all_digits, dec_exp, p, upper, fspec.hash)
      _fmt_num_pad(sign_str, digits, fspec)

    | 'g' | 'G' =>
      let p: USize = match fspec.precision
        | let p2: USize => if p2 == 0 then 1 else p2 end
        else
          full_prec
        end
      // Use scientific if exp < -4 or exp >= p.
      let use_sci: Bool = (dec_exp < -4) or (dec_exp >= p.i64())
      let digits: String val = if use_sci then
          let sci_p: USize = if p > 0 then p - 1 else 0 end
          let s = _fmt_sci(all_digits, dec_exp, sci_p, upper, fspec.hash)
          if fspec.hash then
            s
          else
            _fmt_trim_zeros_sci(s, upper)
          end
        else
          let fix_p_i: ISize = (p.isize() - dec_exp.isize()).max(0)
          let fix_p: USize = fix_p_i.usize()
          let s = _fmt_fixed(all_digits, dec_exp, fix_p, fspec.hash)
          if fspec.hash then
            _fmt_ensure_dot(s)
          else
            _fmt_trim_zeros_fixed(s)
          end
        end
      _fmt_num_pad(sign_str, digits, fspec)

    | '%' =>
      // Shift dec_exp by 2 (multiply by 100) and format as 'f'.
      let p: USize = match fspec.precision
        | let p2: USize => p2
        else
          6
        end
      let digits: String val = _fmt_fixed(all_digits, dec_exp + 2, p, fspec.hash)
      let grouped: String val = _fmt_group_fixed(digits, fspec.grouping)
      let with_pct: String val = recover
        let s = String(grouped.size() + 1)
        s.append(grouped)
        s.push('%')
        s
      end
      _fmt_num_pad(sign_str, with_pct, fspec)

    else
      // Default: general with full precision, no trailing zeros.
      let p: USize = full_prec
      let use_sci: Bool = (dec_exp < -4) or (dec_exp >= p.i64())
      let digits: String val = if use_sci then
          _fmt_trim_zeros_sci(_fmt_sci(all_digits, dec_exp, p - 1, false, false), false)
        else
          let fix_p_i: ISize = (p.isize() - dec_exp.isize()).max(0)
          let fix_p: USize = fix_p_i.usize()
          _fmt_trim_zeros_fixed(_fmt_fixed(all_digits, dec_exp, fix_p, false))
        end
      _fmt_num_pad(sign_str, digits, fspec)
    end


  // ---- private formatting helpers -------------------------------------------

  fun _fmt_round(digits: String val, keep: USize): (String val, Bool) =>
    """
    Round-to-nearest on `digits` (decimal digit chars), keeping `keep` digits.
    Returns `(rounded_digits, carry_out)` where `carry_out` is true when
    rounding propagates past the first digit (e.g. "999" rounded to 1 digit
    gives ("000", true), meaning the caller should add 1 to the integer part).
    When `keep >= digits.size()` no rounding is needed; returns `(digits, false)`
    with the original string (caller pads with zeros as needed).
    """
    if keep >= digits.size() then
      return (digits, false)
    end
    // Check the digit just after the last kept position.
    let next: U8 = try digits(keep)? else '0' end
    if next < '5' then
      // Truncate.
      return (digits.trim(0, keep), false)
    end
    // Round up: copy kept digits into a mutable string, then propagate carry.
    let rounded: String ref = String(keep)
    for i in Range(0, keep) do
      rounded.push(try digits(i)? else '0' end)
    end
    var carry: Bool = true
    var pos: ISize = keep.isize() - 1
    while carry and (pos >= 0) do
      let d = try rounded(pos.usize())? else '0' end
      if d < '9' then
        try rounded(pos.usize())? = d + 1 end
        carry = false
      else
        try rounded(pos.usize())? = '0' end
        pos = pos - 1
      end
    end
    (rounded.clone(), carry)


  fun _fmt_fixed(digits: String val, dec_exp: I64, prec: USize,
    alt: Bool): String val =>
    """
    Build a fixed-point string from `digits` (unsigned decimal mantissa),
    `dec_exp` (digits left of decimal point), and `prec` (digits after point).
    `alt` is the `#` flag — forces a decimal point even with zero fractional digits.
    `digits` may be shorter or longer than needed; missing digits are zeros.
    The last fractional digit is rounded to nearest.
    """
    // Total digits we need: e integer digits + prec fractional digits.
    // Round the source digits to (max(e,0) + prec) significant positions.
    let e = dec_exp.isize()
    let int_digits: USize = if e > 0 then e.usize() else 0 end
    let total_needed: USize = int_digits + prec
    (let rd, let carry) = _fmt_round(digits, total_needed)

    // If rounding carried out of all digits, adjust dec_exp effectively by +1.
    let actual_exp: I64 = if carry then dec_exp + 1 else dec_exp end
    let ae = actual_exp.isize()
    let nd = rd.size()

    recover
      let s = String

      // Integer part.
      if ae <= 0 then
        s.append("0")
      else
        for i in Range[ISize](0, ae) do
          s.push(if (not carry) and (i < nd.isize()) then
              try rd(i.usize())? else '0' end
            elseif carry then
              // carry means the rounded string is all zeros with a leading 1.
              if i == 0 then '1' else '0' end
            else '0' end)
        end
      end

      // Fractional part.
      if prec > 0 then
        s.push('.')
        for i in Range(0, prec) do
          let src_pos: ISize = ae + i.isize()
          let c: U8 =
            if carry then
              '0'
            elseif (src_pos >= 0) and (src_pos < nd.isize()) then
              try rd(src_pos.usize())? else '0' end
            else
              '0'
            end
          s.push(c)
        end
      elseif alt then
        s.push('.')
      end
      s
    end


  fun _fmt_sci(digits: String val, dec_exp: I64, prec: USize,
    upper: Bool, alt: Bool = false): String val =>
    """
    Build a scientific-notation string: `d.dddde±ee`.
    `prec` is the number of digits after the decimal point in the mantissa.
    `alt` is the `#` flag — forces a decimal point even with zero fractional digits.
    The last mantissa digit is rounded to nearest.
    """
    // We need 1 + prec significant digits total.
    let keep = prec + 1
    (let rd, let carry) = _fmt_round(digits, keep)
    // carry means mantissa rounded up to 10.0... → adjust exponent by 1.
    let actual_exp: I64 = if carry then dec_exp + 1 else dec_exp end
    let nd = rd.size()

    recover
      let s = String

      // Leading digit.
      s.push(if carry then '1' elseif nd > 0 then try rd(0)? else '0' end else '0' end)

      // Fractional digits of mantissa.
      if prec > 0 then
        s.push('.')
        for i in Range(0, prec) do
          s.push(
            if carry then '0'
            elseif (i + 1) < nd then try rd(i + 1)? else '0' end
            else '0' end)
        end
      elseif alt then
        s.push('.')
      end

      // Exponent: actual_exp - 1 (mantissa is d.xxx, not 0.dxxx).
      let exp_val: I64 = actual_exp - 1
      s.append(if upper then "E" else "e" end)
      s.append(if exp_val >= 0 then "+" else "-" end)
      let abs_exp: U64 = if exp_val >= 0 then exp_val.u64() else (-exp_val).u64() end
      let exp_str: String = abs_exp.string()
      if exp_str.size() < 2 then
        s.push('0')
      end
      s.append(exp_str)
      s
    end


  fun _fmt_trim_zeros_fixed(s: String val): String val =>
    """
    Remove trailing zeros and a lone decimal point from a fixed string.
    """
    var end_pos = s.size()
    while end_pos > 0 do
      let c = try s(end_pos - 1)? else break end
      if c == '0' then
        end_pos = end_pos - 1
      elseif c == '.' then
        end_pos = end_pos - 1
        break
      else
        break
      end
    end
    if end_pos == s.size() then
      s
    else
      s.trim(0, end_pos)
    end


  fun _fmt_trim_zeros_sci(s: String val, upper: Bool): String val =>
    """
    Remove trailing zeros from the mantissa part of a scientific string,
    stopping at the `e`/`E` separator.
    """
    let e_char: U8 = if upper then 'E' else 'e' end
    var e_pos: USize = s.size()
    for i in Range(0, s.size()) do
      if try s(i)? == e_char else false end then
        e_pos = i
        break
      end
    end
    if e_pos == s.size() then
      return s
    end

    let mantissa: String val = s.trim(0, e_pos)
    let exp_part: String val = s.trim(e_pos)

    let trimmed_m: String val = _fmt_trim_zeros_fixed(mantissa)
    if trimmed_m.size() == mantissa.size() then
      s
    else
      recover
        let r = String(trimmed_m.size() + exp_part.size())
        r.append(trimmed_m)
        r.append(exp_part)
        r
      end
    end


  fun _fmt_ensure_dot(s: String val): String val =>
    """
    Append a trailing '.' if the string contains no decimal point.
    """
    for i in Range(0, s.size()) do
      if try s(i)? == '.' else false end then
        return s
      end
    end
    recover
      let r = String(s.size() + 1)
      r.append(s)
      r.push('.')
      r
    end


  fun _fmt_group_fixed(fixed: String val, grouping: (U8 | None)): String val =>
    """
    Apply thousands grouping to the integer part of a fixed-point string.
    """
    let sep: U8 = match grouping
      | let g: U8 => g
      else
        return fixed
      end
    var dot_pos: USize = fixed.size()
    for i in Range(0, fixed.size()) do
      if try fixed(i)? == '.' else false end then
        dot_pos = i
        break
      end
    end
    let int_part: String val = fixed.trim(0, dot_pos)
    let frac_part: String val = fixed.trim(dot_pos)
    let n = int_part.size()
    if n <= 3 then
      return fixed
    end
    let num_seps = (n - 1) / 3
    recover
      let r = String(n + num_seps + frac_part.size())
      for i in Range(0, n) do
        if (i > 0) and (((n - i) % 3) == 0) then
          r.push(sep)
        end
        try r.push(int_part(i)?) end
      end
      r.append(frac_part)
      r
    end


  fun _fmt_int(fspec: FormatSpec val): String =>
    """
    Format this (assumed integer) MPFRep using integer type codes.
    """
    let base: USize = match fspec.type_char
      | 'b' => 2
      | 'o' => 8
      | 'x' | 'X' => 16
      else
        10
      end
    let upper = fspec.type_char == 'X'

    // Extract absolute integer value from digits.
    // Convert base-256 integer to target base string using the magnitude digits.
    // _digits holds the absolute value; sign does not affect the digit array.
    let raw: String val = _b256_to_base(_digits, base, upper)

    // Precision zero-padding for integers.
    let padded: String val = match fspec.precision
      | let p: USize if raw.size() < p =>
        recover
          let s = String(p)
          for _ in Range(raw.size(), p) do s.push('0') end
          s.append(raw)
          s
        end
      else
        raw
      end

    let prefix: String val =
      if fspec.hash then
        match fspec.type_char
        | 'b' => "0b"
        | 'o' => "0o"
        | 'x' => "0x"
        | 'X' => "0X"
        else "" end
      else "" end

    let sign_str: String val = if _sign then
        "-"
      elseif fspec.sign is SignPlus then
        "+"
      elseif fspec.sign is SignSpace then
        " "
      else
        ""
      end

    let fill_char: U32 = if fspec.zero then U32('0') else fspec.fill end
    let content = sign_str.size() + prefix.size() + padded.size()
    let width = fspec.width

    if content >= width then
      return recover
        let s = String(content)
        s.append(sign_str)
        s.append(prefix)
        s.append(padded)
        s
      end
    end

    let pad = width - content
    recover
      let s = String(width)
      match fspec.align
      | AlignLeft =>
        s.append(sign_str)
        s.append(prefix)
        s.append(padded)
        s.append(String.from_utf32(fill_char) * pad)
      | AlignNumeric =>
        s.append(sign_str)
        s.append(prefix)
        s.append(String.from_utf32(fill_char) * pad)
        s.append(padded)
      | AlignCenter =>
        let before = pad / 2
        s.append(String.from_utf32(fill_char) * before)
        s.append(sign_str)
        s.append(prefix)
        s.append(padded)
        s.append(String.from_utf32(fill_char) * (pad - before))
      else
        s.append(String.from_utf32(fill_char) * pad)
        s.append(sign_str)
        s.append(prefix)
        s.append(padded)
      end
      s
    end


  fun _b256_to_base(digits: Array[U8] val, base: USize, upper: Bool): String val =>
    """
    Convert a big-endian base-256 integer to the target base string.
    """
    ifdef debug then
      (base <= 36) or Fail(Format("[_b256_to_base] Base for conversion ({}) must be in range [2, 36]", base))
    end
    if digits.size() == 0 then
      return "0"
    end
    let charset: String val = if upper then
        "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      else
        "0123456789abcdefghijklmnopqrstuvwxyz"
      end
    // Work with a mutable copy.
    let d: Array[U8] ref = recover
      let a = Array[U8](digits.size())
      for b in digits.values() do
        a.push(b)
      end
      a
    end
    var result = recover String end
    while not _b256_is_zero(d) do
      let rem = _b256_short_div(d, base.u32())
      try result.unshift(charset(rem.usize())?) end
    end
    if result.size() == 0 then
      result.push('0')
    end
    consume result


  fun tag _b256_is_zero(d: Array[U8] ref): Bool =>
    """
    `true` when all digits in `d` are zero.
    """"
    for b in d.values() do
      if b != 0 then
        return false
      end
    end
    true


  fun tag _b256_short_div(d: Array[U8] ref, divisor: U32): U32 =>
    """
    Divide big-endian base-256 array in place by divisor; return remainder.
    """
    var rem: U32 = 0
    for i in Range(0, d.size()) do
      let cur: U32 = (rem * 256) + (try d(i)?.u32() else 0 end)
      try d(i)? = (cur / divisor).u8() end
      rem = cur % divisor
    end
    rem


  fun _fmt_pad(s: String val, fspec: FormatSpec val): String =>
    """
    Apply width/fill/alignment to an already-complete string (NaN, Inf).
    """
    let width = fspec.width
    if s.size() >= width then
      return s
    end
    let pad = width - s.size()
    recover
      let r = String(width)
      match fspec.align
      | AlignLeft =>
        r.append(s); r.append(String.from_utf32(fspec.fill) * pad)
      | AlignCenter =>
        let before = pad / 2
        r.append(String.from_utf32(fspec.fill) * before)
        r.append(s)
        r.append(String.from_utf32(fspec.fill) * (pad - before))
      else
        r.append(String.from_utf32(fspec.fill) * pad); r.append(s)
      end
      r
    end


  fun _fmt_num_pad(sign_str: String val, digits: String val,
    fspec: FormatSpec val): String =>
    """
    Apply width/fill/alignment keeping sign and digits separate.
    """
    let fill_char: U32 = if fspec.zero then U32('0') else fspec.fill end
    let content = sign_str.size() + digits.size()
    let width = fspec.width
    if content >= width then
      return recover
        let s = String(content)
        s.append(sign_str)
        s.append(digits)
        s
      end
    end
    let pad = width - content
    recover
      let s = String(width)
      match fspec.align
      | AlignLeft =>
        s.append(sign_str)
        s.append(digits)
        s.append(String.from_utf32(fill_char) * pad)
      | AlignNumeric =>
        s.append(sign_str)
        s.append(String.from_utf32(fill_char) * pad)
        s.append(digits)
      | AlignCenter =>
        let before = pad / 2
        s.append(String.from_utf32(fill_char) * before)
        s.append(sign_str)
        s.append(digits)
        s.append(String.from_utf32(fill_char) * (pad - before))
      else
        s.append(String.from_utf32(fill_char) * pad)
        s.append(sign_str)
        s.append(digits)
      end
      s
    end

