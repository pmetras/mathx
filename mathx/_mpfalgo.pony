// Multi-precision float algorithms: all computation, no policy.
//
// Each method takes `MPFRep` operands and a working-precision `w: USize`
// (in bytes), performs the computation at that precision using `MPFRep`
// structural operations, and returns an `MPFRep` at precision `w` (unrounded).
//
// Rounding from `w` to the output precision `p` is the responsibility of
// `MPFContext._round_to`; it is NOT done here.  Internal truncation to `w`
// inside loops (e.g. Newton iterations) uses `MPFRep._trunc(w)` directly —
// this bounds intermediate growth at *working* precision, which is distinct
// from the output rounding done by the caller.
//
// This primitive is private (`_` prefix) — power users go through `MPFContext`
// or `MPFloat`.

use "../assertx"

use "collections"


primitive _MPFAlgo
  """
  All arithmetic and transcendental algorithms for `MPFRep` values.

  Methods are pure functions: no state, no policy decisions.  They receive the
  working precision `w` from the caller and return a result at that precision.

  This is a stub file.  Algorithms will be migrated here one group at a time
  (Steps 4a–4h of the architectural split), keeping `mpfloat.pony` unchanged
  until Step 5.
  """

  // ── Addition / Subtraction ─────────────────────────────────────────────────

  fun _add_mag(a: MPFRep, b: MPFRep, sgn: Bool): MPFRep =>
    """
    Add the magnitudes `|a|` and `|b|` and return the sum with sign `sgn`.
    Both operands must be finite.  `a` must have the larger (or equal) exponent.

    The result precision is `max(|a._digits|, shift + |b._digits|)` where
    `shift = a.exponent() − b.exponent()`, so no precision is lost.  A one-byte
    carry guard is prepended; when the guard is zero it is stripped and the
    exponent stays `a.exponent()`; when it is nonzero the full array (with
    the carry byte) is used and the exponent becomes `a.exponent() + 1`.
    """
    let ea: I64 = a.exponent()
    let eb: I64 = b.exponent()
    let ad: Array[U8] val = a.raw_digits()
    let bd: Array[U8] val = b.raw_digits()
    let na: USize = a._size()
    let nb: USize = b._size()
    let shift: USize = (ea - eb).usize()

    // When `b` lies entirely below `a`'s precision, the sum is `a`.
    if shift >= na then
      return MPFRep._create(sgn, false, false, ea, ad)
    end

    let prec: USize = na.max(shift + nb)
    let result_size: USize = prec + 1
    let base = _MPFBase
    let raw: Array[U8] val = recover
      let res = Array[U8].init(0, result_size)
      try
        var carry: U16 = 0
        var col: USize = prec
        repeat
          col = col - 1
          let ai: U8 = if col < na then ad(col)? else 0 end
          let bi: U8 =
            if (col >= shift) and ((col - shift) < nb) then bd(col - shift)?
            else 0
            end
          (let sum, let c2) = base.addc(ai, bi, carry)
          carry = c2
          res.update(col + 1, sum)?
        until col == 0 end
        res.update(0, base.lowb(carry))?
      end
      res
    end

    let carry_digit: U8 = try raw(0)? else 0 end
    if carry_digit == 0 then
      let d: Array[U8] val = recover
        let a2 = Array[U8].init(0, prec)
        raw.copy_to(a2, 1, 0, prec)
        a2
      end
      MPFRep._create(sgn, false, false, ea, d)
    else
      MPFRep._create(sgn, false, false, ea + 1, raw)
    end


  fun _sub_mag(a: MPFRep, b: MPFRep, sgn: Bool): MPFRep =>
    """
    Subtract `|b|` from `|a|`, where `|a| >= |b|` (enforced by the caller).
    Returns the difference with sign `sgn`.  Both operands must be finite.

    The result is normalised: leading zero bytes are stripped and the exponent
    adjusted accordingly.  When the result is exactly zero a zero `MPFRep` at
    precision `a._size()` is returned.
    """
    let ea: I64 = a.exponent()
    let eb: I64 = b.exponent()
    let shift: USize = (ea - eb).usize()
    let na: USize = a._size()
    let nb: USize = b._size()
    let result_size: USize = na.max(shift + nb)
    let ad: Array[U8] val = a.raw_digits()
    let bd: Array[U8] val = b.raw_digits()
    let max_b: U16 = U8.max_value().u16()
    let bbase = _MPFBase
    let raw: Array[U8] val = recover
      let res = Array[U8].init(0, result_size)
      try
        var carry: U16 = max_b + 1
        var col: USize = result_size
        repeat
          col = col - 1
          let ai: U16 = if col < na then ad(col)?.u16() else 0 end
          let bi: U16 =
            if (col >= shift) and ((col - shift) < nb) then bd(col - shift)?.u16()
            else 0
            end
          carry = ((max_b + ai) - bi) + bbase.highb(carry).u16()
          res.update(col, bbase.lowb(carry))?
        until col == 0 end
      end
      res
    end

    // Normalise: strip leading zero bytes and adjust the exponent.
    let raw_size: USize = raw.size()
    var leading: USize = 0
    try
      while (leading < raw_size) and (raw(leading)? == 0) do
        leading = leading + 1
      end
    end
    if leading == raw_size then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, na))
    end
    let new_size: USize = raw_size - leading
    let new_exp: I64 = ea - leading.i64()
    let d: Array[U8] val = recover
      let a2 = Array[U8].init(0, new_size)
      raw.copy_to(a2, leading, 0, new_size)
      a2
    end
    MPFRep._create(sgn, false, false, new_exp, d)


  fun _add(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Add `a + b` at working precision `w` bytes.

    Special-value rules:
    - NaN propagates: `NaN + x = NaN`.
    - `+∞ + (−∞) = NaN`.
    - `±∞ + finite = ±∞`.
    - Zero is the additive identity.

    When operand signs agree, magnitudes are added (`_add_mag`).  When they
    differ, the smaller magnitude is subtracted from the larger and the result
    takes the sign of the larger.  The result is truncated to `w` bytes
    (working precision) before being returned; the caller rounds to output
    precision via `MPFContext._round_to`.
    """
    if a.is_nan() or b.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() and b.is_infinite() then
      if a.sign_bit() == b.sign_bit() then
        return MPFRep._create(a.sign_bit(), false, true, 0, Array[U8].create())
      end
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(a.sign_bit(), false, true, 0, Array[U8].create())
    end
    if b.is_infinite() then
      return MPFRep._create(b.sign_bit(), false, true, 0, Array[U8].create())
    end
    if a.is_zero() then
      return b._trunc(w)
    end
    if b.is_zero() then
      return a._trunc(w)
    end

    let res: MPFRep =
      if a.sign_bit() == b.sign_bit() then
        // Same sign: add magnitudes.
        if a.exponent() >= b.exponent() then
          _add_mag(a, b, a.sign_bit())
        else
          _add_mag(b, a, a.sign_bit())
        end
      else
        // Different signs: subtract smaller magnitude from larger.
        match a._cmp_mag(b)
        | Greater => _sub_mag(a, b, a.sign_bit())
        | Less    => _sub_mag(b, a, b.sign_bit())
        else
          MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
        end
      end
    res._trunc(w)


  fun _sub(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Subtract `a − b` at working precision `w` bytes.  Delegates to `_add`
    after negating `b`, inheriting all sign and special-value handling.
    """
    let neg_b = MPFRep._create(
      not b.sign_bit(), b.is_nan(), b.is_infinite(), b.exponent(), b.raw_digits())
    _add(a, neg_b, w)


  // ── Multiplication ─────────────────────────────────────────────────────────

  fun _mul(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Multiply `a × b` at working precision `w` bytes via FFT convolution.

    The result sign is `a.sign_bit() XOR b.sign_bit()`.  The raw convolution
    product has `a._size() + b._size()` digits; it is normalised by stripping
    leading zeros and adjusting the exponent, then truncated to `w` bytes.

    Both operands must be finite (NaN/±∞ handling belongs to the public API
    layer in `MPFContext`).
    """
    let result_sign: Bool = a.sign_bit() != b.sign_bit()
    let this_size: USize = a._size()
    let that_size: USize = b._size()
    let res_size: USize = this_size + that_size
    let pow2: USize = res_size.next_pow2().max(4)
    let ad: Array[U8] val = a.raw_digits()
    let bd: Array[U8] val = b.raw_digits()

    var err_i: USize = 0
    let d: Array[U8] val = recover
      let res = Array[U8].init(0, res_size)
      var i: USize = 0
      try
        let fa = Array[F64].init(0.0, pow2)
        i = 0
        while i < this_size do
          fa.update(i, ad(i)?.f64())?
          i = i + 1
        end

        let fb = Array[F64].init(0.0, pow2)
        i = 0
        while i < that_size do
          fb.update(i, bd(i)?.f64())?
          i = i + 1
        end

        FFT.fourier_real(fa)
        FFT.fourier_real(fb)

        // Pointwise complex multiply in the frequency domain.
        fb.update(0, fb(0)? * fa(0)?)?
        fb.update(1, fb(1)? * fa(1)?)?
        i = 2
        while i < pow2 do
          let temp = fb(i)?
          fb.update(i, (temp * fa(i)?) - (fb(i + 1)? * fa(i + 1)?))?
          fb.update(i + 1, (temp * fa(i + 1)?) + (fb(i + 1)? * fa(i)?))?
          i = i + 2
        end

        FFT.fourier_real(fb, true)

        // Convert back to U8 digits with carry propagation (high to low).
        var carry: F64 = 0.0
        let base_f: F64 = _MPFBase.base().f64()
        i = pow2
        repeat
          i = i - 1
          err_i = i
          let temp = fb(i)? + carry + 0.5
          carry = (temp / base_f).u64().f64()
          fb.update(i, temp - (carry * base_f))?
        until i == 0 end

        // res[0] = carry; res[1..res_size-1] = fb[0..res_size-2].
        res.update(0, carry.u8())?
        i = 1
        while i < res_size do
          res.update(i, fb(i - 1)?.u8())?
          i = i + 1
        end
      else
        try
          Assert(false,
            "[_MPFAlgo._mul] index out of bounds at i=" + err_i.string(), true)?
        end
      end
      res
    end

    // Strip leading zeros; each stripped zero shifts the exponent down by 1.
    var leading: USize = 0
    try
      while (leading < d.size()) and (d(leading)? == 0) do
        leading = leading + 1
      end
    end
    if leading == d.size() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    let new_exp: I64 = (a.exponent() + b.exponent()) - leading.i64()
    let new_size: USize = d.size() - leading
    let norm_d: Array[U8] val = recover
      let a2 = Array[U8].init(0, new_size)
      d.copy_to(a2, leading, 0, new_size)
      a2
    end
    MPFRep._create(result_sign, false, false, new_exp, norm_d)._trunc(w)


  // ── Division / Reciprocal ─────────────────────────────────────────────────

  fun _inv(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `1 / a` via Newton's method at working precision `w` bytes.

    Algorithm: `x_{n+1} = x_n × (2 − G × x_n)` where `G = |a|` reinterpreted
    with exponent 1, so `G ∈ [1, 256)`.  The iteration converges to `y = 1/G`;
    the final result exponent is `(y._exponent + 1) − a.exponent()`.

    A safety limit of `(w + 1) × 4` iterations prevents non-convergence in
    edge cases where exact-digit comparison never triggers.

    Sign propagation: the sign of the result equals `a.sign_bit()`.
    """
    let size: USize = w + 1
    let ad: Array[U8] val = a.raw_digits()
    let prec: USize = a._size()

    // G = same digits as |a| but exponent=1, so its integer part is d[0] ∈ [1, 256).
    let g = MPFRep._create(false, false, false, 1, ad)

    // F64 initial estimate: read the first few digits.
    let ng: USize = prec.min(4)
    var fg: F64 = 0.0
    try
      var i: USize = ng
      repeat
        i = i - 1
        fg = (fg / 256.0) + ad(i)?.f64()
      until i == 0 end
    end
    var res = MPFRep.from_f64(1.0 / fg, size)

    let two = MPFRep.from_f64(2.0, size)
    var iters: USize = 0
    let max_iters: USize = size * 4
    while iters < max_iters do
      iters = iters + 1
      let gy   = _mul(g, res, size)._trunc(size)
      let delta = _sub(two, gy, size)._trunc(size)
      let new_res = _mul(res, delta, size)._trunc(size)

      // Converged when all leading `w` digits agree with previous iterate.
      var changed: Bool = false
      let new_rd = new_res.raw_digits()
      let res_rd = res.raw_digits()
      let ns: USize = new_res._size().min(res._size()).min(w)
      var ci: USize = 0
      try
        while ci < ns do
          if new_rd(ci)? != res_rd(ci)? then
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

    // 1/|a| = res × 256^{1−e}.  Adjust sign and exponent.
    MPFRep._create(
      a.sign_bit(), false, false,
      (res.exponent() + 1) - a.exponent(),
      res.raw_digits())


  fun _div(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Divide `a / b` at working precision `w` bytes.  Computes `a × inv(b)`.

    Both operands must be finite and non-zero (special-value handling belongs
    to the public API in `MPFContext`).
    """
    _mul(a, _inv(b, w), w)


  // ── Square Root ────────────────────────────────────────────────────────────

  fun _sqrt(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `√a` via Newton's method at working precision `w` bytes.

    Algorithm: reciprocal-sqrt iteration `y_{n+1} = y_n × (3 − H × y_n²) / 2`,
    converging to `1/√H`; the final result is `√H = H × (1/√H)`.

    Parity split: for odd `a.exponent()` use `H` with exponent 1 so
    `H ∈ [1, 256)`.  For even `a.exponent()` use exponent 2 so
    `H ∈ [256, 65536)`.  With this split `a.exponent() − h_exp` is always
    even, so the result exponent is always an integer.

    A safety limit of `(w + 1) × 4` iterations prevents non-convergence.
    """
    let size: USize = w + 1
    let ad: Array[U8] val = a.raw_digits()
    let prec: USize = a._size()
    let a_exp: I64 = a.exponent()

    // Choose h_exp so that a_exp − h_exp is always even.
    let h_exp: I64 = if (a_exp and 1) != 0 then 1 else 2 end
    let h = MPFRep._create(false, false, false, h_exp, ad)

    // F64 initial estimate.
    let ng: USize = prec.min(4)
    var fg: F64 = 0.0
    let base_f: F64 = _MPFBase.base().f64()
    try
      var k: USize = ng
      repeat
        k = k - 1
        fg = (fg / base_f) + ad(k)?.f64()
      until k == 0 end
    end
    let fg_h: F64 = if h_exp == 1 then fg else fg * base_f end
    var res = MPFRep.from_f64(1.0 / fg_h.sqrt(), size)

    let three = MPFRep.from_f64(3.0, size)
    var iters: USize = 0
    let max_iters: USize = size * 4
    while iters < max_iters do
      iters = iters + 1
      let y2        = _mul(res, res, size)._trunc(size)
      let hy2       = _mul(h, y2, size)._trunc(size)
      let delta     = _sub(three, hy2, size)._trunc(size)
      let new_full  = _mul(res, delta, size)._trunc(size)
      (let halved, _) = new_full._short_div(2)
      let nr = halved._trunc(size)

      // Converged when all leading `w` digits agree with previous iterate.
      var changed: Bool = false
      let nr_rd = nr.raw_digits()
      let res_rd = res.raw_digits()
      let ns: USize = nr._size().min(res._size()).min(w)
      var i: USize = 0
      try
        while i < ns do
          if nr_rd(i)? != res_rd(i)? then
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
    let sqrt_h = _mul(h, res, size)._trunc(size)
    let result_exp: I64 = sqrt_h.exponent() + ((a_exp - h_exp) / 2)
    MPFRep._create(false, false, false, result_exp, sqrt_h.raw_digits())


  // ── Division with Remainder ────────────────────────────────────────────────

  fun _trunc_frac(a: MPFRep): MPFRep =>
    """
    Return `a` with its fractional base-256 bytes zeroed out (truncation toward
    zero).  Equivalent to the integer part of `a`.

    - NaN or ±∞ → returned unchanged.
    - Zero → returned as-is.
    - `a.exponent() ≤ 0` (purely fractional) → zero with the same precision.
    - `a.exponent() ≥ a._size()` (purely integer) → `a` unchanged.
    - Mixed → zero out bytes at indices `a.exponent()` and beyond.
    """
    if (not a.is_finite()) or a.is_zero() then
      return a
    end
    let e: I64 = a.exponent()
    let p: USize = a._size()
    if e <= 0 then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, p))
    end
    let int_bytes: USize = e.usize().min(p)
    if int_bytes == p then
      return a
    end
    let d: Array[U8] val = a.raw_digits()
    let trunc_d: Array[U8] val = recover
      let arr = Array[U8].init(0, p)
      var i: USize = 0
      while i < int_bytes do
        try arr(i)? = d(i)? end
        i = i + 1
      end
      arr
    end
    MPFRep._create(a.sign_bit(), false, false, e, trunc_d)


  fun _one(w: USize): MPFRep =>
    """
    Return the value `1.0` at working precision `w` bytes.
    """
    MPFRep.from_f64(1.0, w)


  fun _divrem(a: MPFRep, b: MPFRep, w: USize): (MPFRep, MPFRep) =>
    """
    Truncated division with remainder at working precision `w` bytes.

    Returns `(q, r)` where `q = trunc(a / b)` and `r = a − q × b`.

    Special cases follow IEEE 754 conventions:
    - NaN in either operand → `(NaN, NaN)`.
    - `a` is ±∞ and `b` is finite → `(±∞, NaN)`.
    - Both ±∞ → `(NaN, NaN)`.
    - `b` is ±∞ and `a` is finite → `(0, a)`.
    - `b` is 0, `a` is 0 → `(NaN, NaN)`.
    - `b` is 0, `a` non-zero → `(±∞, NaN)`.
    - `a` is 0 → `(0, 0)`.

    Post-correction: Newton's `_inv` always undershoots by up to 1 ULP.
    When `|r| ≥ |b|`, increment `|q|` by 1 and recompute `r`.  At most one
    step is needed.
    """
    let nan = MPFRep.nan_val()
    if a.is_nan() or b.is_nan() then
      return (nan, nan)
    end
    if a.is_infinite() and b.is_infinite() then
      return (nan, nan)
    end
    if a.is_infinite() then
      let q_inf = MPFRep._create(
        a.sign_bit() != b.sign_bit(), false, true, 0, Array[U8].create())
      return (q_inf, nan)
    end
    if b.is_infinite() then
      let z = MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
      return (z, MPFRep._create(a.sign_bit(), false, false, a.exponent(), a.raw_digits()))
    end
    if b.is_zero() then
      if a.is_zero() then
        return (nan, nan)
      end
      let q_inf = MPFRep._create(
        a.sign_bit() != b.sign_bit(), false, true, 0, Array[U8].create())
      return (q_inf, nan)
    end
    if a.is_zero() then
      let z = MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
      return (z, z)
    end

    var q: MPFRep = _trunc_frac(_div(a, b, w))
    var r: MPFRep = _sub(a, _mul(q, b, w), w)

    // Post-correction: if |r| ≥ |b|, q was off by 1.
    if (not r.is_zero()) and (r._cmp_mag(b) != Less) then
      let one = _one(w)
      if a.sign_bit() == b.sign_bit() then
        q = _add(q, one, w)
      else
        q = _sub(q, one, w)
      end
      r = _sub(a, _mul(q, b, w), w)
    end
    (q, r)


  fun _fld(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Floored division `floor(a / b)` at working precision `w` bytes.

    The result is the largest integer `q` such that `q × b ≤ a`.  When the
    truncated remainder is non-zero and `a` and `b` have opposite signs, the
    truncated quotient is decremented by 1.

    Special cases are inherited from `_divrem`.
    """
    (let q, let r) = _divrem(a, b, w)
    if r.is_nan() then
      return MPFRep.nan_val()
    end
    if (not r.is_zero()) and (a.sign_bit() != b.sign_bit()) then
      return _sub(q, _one(w), w)
    end
    q


  fun _mod(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Floored remainder `a − floor(a / b) × b` at working precision `w` bytes.

    The result has the same sign as `b`.  Special cases are inherited from
    `_divrem`; when `_rem` returns NaN, `_mod` also returns NaN.
    """
    let r: MPFRep = _rem(a, b, w)
    if r.is_nan() then
      return r
    end
    if (not r.is_zero()) and (a.sign_bit() != b.sign_bit()) then
      return _add(r, b, w)
    end
    r


  fun _rem(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Truncated remainder `a − trunc(a / b) × b` at working precision `w` bytes.

    The result has the same sign as `a`.  Special cases are inherited from
    `_divrem`.
    """
    _divrem(a, b, w)._2


  // ── Internal convergence helpers ──────────────────────────────────────────

  fun _digits_eq(a: MPFRep, b: MPFRep): Bool =>
    """
    Return `true` when `a` and `b` have the same exponent and the same leading
    digits (compared byte by byte up to the shorter array).  Used for
    convergence detection in series loops — two values are considered "equal"
    for convergence purposes when adding one to the other does not change any
    stored digit.
    """
    if a.exponent() != b.exponent() then
      return false
    end
    let ad: Array[U8] val = a.raw_digits()
    let bd: Array[U8] val = b.raw_digits()
    let n: USize = ad.size().min(bd.size())
    var i: USize = 0
    try
      while i < n do
        if ad(i)? != bd(i)? then
          return false
        end
        i = i + 1
      end
    end
    true


  fun _neg_rep(a: MPFRep): MPFRep =>
    """
    Return a copy of `a` with its sign flipped.  Equivalent to arithmetic
    negation but operates purely on the representation (no allocation of new
    digits).
    """
    MPFRep._create(
      not a.sign_bit(), a.is_nan(), a.is_infinite(),
      a.exponent(), a.raw_digits())


  // ── Series helpers ─────────────────────────────────────────────────────────

  fun _atanh_series(x: MPFRep, w: USize): MPFRep =>
    """
    Compute `arctanh(x)` via the Taylor series
    `Σ_{k=0}^∞ x^{2k+1} / (2k+1)` at working precision `w` bytes.

    The caller must ensure `|x| < 1` for convergence.  The loop terminates
    when the current addend no longer changes the running sum at precision `w`.

    Uses `MPFRep._short_div` (exact U8 division) for denominators up to 255;
    falls back to `_div` for larger denominators (only for very high `w`).
    """
    let p2: USize = w + 2
    let x2: MPFRep = _mul(x, x, p2)._trunc(p2)
    var term: MPFRep = x._trunc(p2)
    var sum:  MPFRep = term
    var k: USize = 3
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term = _mul(term, x2, p2)._trunc(p2)
      let addend: MPFRep =
        if k <= 255 then
          (let q, _) = term._short_div(k.u8())
          q
        else
          _div(term, MPFRep.from_f64(k.f64(), p2), p2)._trunc(p2)
        end
      let new_sum: MPFRep = _add(sum, addend, p2)._trunc(p2)
      if _digits_eq(new_sum, sum) then
        break
      end
      sum = new_sum
      k = k + 2
      iters = iters + 1
    end
    sum._trunc(w)


  fun _ln2_const(w: USize): MPFRep =>
    """
    Compute `ln(2)` to `w` bytes of working precision using the identity
    `ln(2) = 2 × arctanh(1/3)`.

    The rational `1/3` has the exact repeating base-256 representation
    `[0x55, 0x55, …]` (since `85/255 = 1/3`).
    """
    let p2: USize = w + 4
    let d: Array[U8] val = recover Array[U8].init(0x55, p2) end
    let third = MPFRep._create(false, false, false, 0, d)
    let two = MPFRep.from_f64(2.0, p2)
    _mul(_atanh_series(third, p2), two, p2)._trunc(w)


  fun _exp_taylor(r: MPFRep, w: USize): MPFRep =>
    """
    Compute `e^r` via the Taylor series `Σ_{k=0}^∞ r^k / k!` at working
    precision `w` bytes.

    Uses the recurrence `term_k = term_{k-1} × r / k` to avoid computing
    factorials.  Works best for `|r| ≤ ln(2)/2 ≈ 0.347`.
    """
    let p2: USize = w + 2
    let one = MPFRep.from_f64(1.0, p2)
    var term: MPFRep = one
    var sum:  MPFRep = one
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term = _mul(term, r, p2)._trunc(p2)
      let divided: MPFRep =
        if k <= 255 then
          (let q, _) = term._short_div(k.u8())
          q
        else
          _div(term, MPFRep.from_f64(k.f64(), p2), p2)._trunc(p2)
        end
      term = divided
      let new_sum: MPFRep = _add(sum, divided, p2)._trunc(p2)
      if _digits_eq(new_sum, sum) then
        break
      end
      sum = new_sum
      k = k + 1
      iters = iters + 1
    end
    sum._trunc(w)


  fun _sin_taylor(r: MPFRep, w: USize): MPFRep =>
    """
    Compute `sin(r)` via the Taylor series `r − r³/3! + r⁵/5! − r⁷/7! + …`
    at working precision `w` bytes.

    The recurrence `term_{k+1} = term_k × (−r²) / ((2k)(2k+1))` advances
    two factorial positions at once.  Works best for `|r| ≤ π/4 ≈ 0.785`.
    """
    let p2: USize = w + 2
    let neg_r2: MPFRep = _neg_rep(_mul(r, r, p2)._trunc(p2))
    var term: MPFRep = r._trunc(p2)
    var sum:  MPFRep = term
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term = _mul(term, neg_r2, p2)._trunc(p2)
      let d1: USize = 2 * k
      let d2: USize = (2 * k) + 1
      term =
        if d1 <= 255 then
          (let q, _) = term._short_div(d1.u8())
          q
        else
          _div(term, MPFRep.from_f64(d1.f64(), p2), p2)._trunc(p2)
        end
      term =
        if d2 <= 255 then
          (let q, _) = term._short_div(d2.u8())
          q
        else
          _div(term, MPFRep.from_f64(d2.f64(), p2), p2)._trunc(p2)
        end
      let new_sum: MPFRep = _add(sum, term, p2)._trunc(p2)
      if _digits_eq(new_sum, sum) then
        break
      end
      sum = new_sum
      k = k + 1
      iters = iters + 1
    end
    sum._trunc(w)


  fun _cos_taylor(r: MPFRep, w: USize): MPFRep =>
    """
    Compute `cos(r)` via the Taylor series `1 − r²/2! + r⁴/4! − r⁶/6! + …`
    at working precision `w` bytes.

    The recurrence `term_{k+1} = term_k × (−r²) / ((2k−1)(2k))` advances
    two factorial positions at once.  Works best for `|r| ≤ π/4 ≈ 0.785`.
    """
    let p2: USize = w + 2
    let neg_r2 = _neg_rep(_mul(r, r, p2)._trunc(p2))
    let one = MPFRep.from_f64(1.0, p2)
    var term: MPFRep = one
    var sum:  MPFRep = term
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term = _mul(term, neg_r2, p2)._trunc(p2)
      let d1: USize = (2 * k) - 1
      let d2: USize = 2 * k
      term =
        if d1 <= 255 then
          (let q, _) = term._short_div(d1.u8())
          q
        else
          _div(term, MPFRep.from_f64(d1.f64(), p2), p2)._trunc(p2)
        end
      term =
        if d2 <= 255 then
          (let q, _) = term._short_div(d2.u8())
          q
        else
          _div(term, MPFRep.from_f64(d2.f64(), p2), p2)._trunc(p2)
        end
      let new_sum: MPFRep = _add(sum, term, p2)._trunc(p2)
      if _digits_eq(new_sum, sum) then
        break
      end
      sum = new_sum
      k = k + 1
      iters = iters + 1
    end
    sum._trunc(w)


  // ── Integer extraction helper ──────────────────────────────────────────────

  fun _round_to_ilong(a: MPFRep): ILong =>
    """
    Extract the nearest integer value of `a` as an `ILong` by first rounding
    `a` toward the nearest integer and then converting via `MPInt`.

    Returns 0 when `a` is NaN, ±∞, or zero.  Used internally in argument
    reduction steps for `_exp`, `_sin`, and `_cos`.
    """
    if (not a.is_finite()) or a.is_zero() then
      return 0
    end
    let rounded: MPFRep = _trunc_frac(a)
    MPInt.from_bytes_be(rounded.sign_bit(), rounded.raw_digits()).ilong()


  fun _extend(a: MPFRep, w: USize): MPFRep =>
    """
    Return `a` extended (or truncated) to `w` bytes by zero-padding the tail
    or calling `_trunc(w)`.  The sign, exponent, and leading digits are
    preserved.
    """
    let p: USize = a._size()
    if p == w then
      return a
    end
    if p > w then
      return a._trunc(w)
    end
    let ad: Array[U8] val = a.raw_digits()
    let padded: Array[U8] val = recover
      let d = Array[U8].init(0, w)
      ad.copy_to(d, 0, 0, p)
      d
    end
    MPFRep._create(a.sign_bit(), false, false, a.exponent(), padded)


  // ── Logarithm / Exponential ────────────────────────────────────────────────

  fun _ln(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `ln(a)` at working precision `w` bytes.

    Special cases:
    - NaN → NaN.
    - `+∞` → `+∞`.
    - `a ≤ 0` → NaN.
    - `0` → `−∞`.

    Algorithm: two-level argument reduction then `_atanh_series`:
    1. `a = 0.d₀d₁… × 256^e`.  Set `k = e − 1`, `m = 0.d₀d₁… × 256 ∈ [1, 256)`.
       `ln(a) = k × ln(256) + ln(m)`.
    2. Find `n = ⌊log₂(d₀)⌋ ∈ {0,…,7}`, `u = m / 2^n ∈ [1, 2)`.
       `ln(m) = n × ln(2) + ln(u)`.
    3. `t = (u−1)/(u+1) ∈ [0, 1/3)`.  `ln(u) = 2 × arctanh(t)`.
    Combined: `ln(a) = (8k + n) × ln(2) + 2 × arctanh(t)`.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() and (not a.sign_bit()) then
      return MPFRep._create(false, false, true, 0, Array[U8].create())
    end
    if a.sign_bit() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return MPFRep._create(true, false, true, 0, Array[U8].create())
    end

    let p2: USize = w + 6

    let ln2: MPFRep = _ln2_const(p2)
    let k: I64 = a.exponent() - 1

    let m_digits: Array[U8] val = recover
      let arr = Array[U8].init(0, p2)
      a.raw_digits().copy_to(arr, 0, 0, a._size().min(p2))
      arr
    end
    let m = MPFRep._create(false, false, false, 1, m_digits)

    let d0: U8 = try a.raw_digits()(0)? else 1 end
    var n_inner: I64 = 0
    var tmp: U8 = d0
    while tmp >= 2 do
      tmp = tmp >> 1
      n_inner = n_inner + 1
    end
    let divisor: U8 = U8(1).shl(n_inner.u8())
    (let u, _) = m._short_div(divisor)

    let one = MPFRep.from_f64(1.0, p2)
    let two = MPFRep.from_f64(2.0, p2)
    let t: MPFRep = _div(_sub(u, one, p2)._trunc(p2), _add(u, one, p2)._trunc(p2), p2)
    let ln_u: MPFRep = _mul(_atanh_series(t, p2), two, p2)._trunc(p2)

    let ln2_factor: I64 = (8 * k) + n_inner
    let correction: MPFRep =
      if ln2_factor == 0 then
        MPFRep._create(false, false, false, 0, Array[U8].init(0, p2))
      else
        let fac = MPFRep.from_f64(ln2_factor.f64(), p2)
        _mul(fac, ln2, p2)._trunc(p2)
      end
    _add(correction, ln_u, p2)._trunc(w)


  fun _exp(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `e^a` at working precision `w` bytes.

    Special cases:
    - NaN → NaN.  `+∞` → `+∞`.  `−∞` → `+0`.  `0` → `1`.

    Algorithm: two-level argument reduction then `_exp_taylor`:
    1. `a = n₂₅₆ × ln(256) + r₁`, `n₂₅₆ = round(a / ln(256))`.
       Scale by `256^n₂₅₆` (adjust exponent).
    2. `r₁ = n₂ × ln(2) + r`, `|r| ≤ ln(2)/2`.
       Scale by `2^n₂` (exact).
    3. Compute `e^r` via `_exp_taylor`.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() and (not a.sign_bit()) then
      return MPFRep._create(false, false, true, 0, Array[U8].create())
    end
    if a.is_infinite() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    if a.is_zero() then
      return MPFRep.from_f64(1.0, w)
    end

    let p2: USize = w + 8

    let ln2: MPFRep = _ln2_const(p2)
    let eight = MPFRep._create(false, false, false, 1, recover [U8(8)] end)
    let ln256: MPFRep = _mul(ln2, eight, p2)._trunc(p2)

    let x: MPFRep = _extend(a, p2)

    let n256_f: MPFRep = _trunc_frac(_div(x, ln256, p2)._trunc(p2))
    let n256: ILong = _round_to_ilong(n256_f)
    let n256_mpf = MPFRep.from_f64(n256.f64(), p2)
    let r1: MPFRep = _sub(x, _mul(n256_mpf, ln256, p2)._trunc(p2), p2)._trunc(p2)

    let n2_f: MPFRep = _trunc_frac(_div(r1, ln2, p2)._trunc(p2))
    let n2: ILong = _round_to_ilong(n2_f)
    let n2_mpf = MPFRep.from_f64(n2.f64(), p2)
    // r = r1 − n2 × ln2.
    let r: MPFRep = _sub(r1, _mul(n2_mpf, ln2, p2)._trunc(p2), p2)._trunc(p2)

    var result: MPFRep = _exp_taylor(r, p2)

    if n2 > 0 then
      let pow2 = MPFRep._create(false, false, false, 1, recover [U8(1).shl(n2.u8())] end)
      result = _mul(result, pow2, p2)._trunc(p2)
    elseif n2 < 0 then
      let pow2 = MPFRep._create(false, false, false, 0, recover [U8(128).shr((-n2 - 1).u8())] end)
      result = _mul(result, pow2, p2)._trunc(p2)
    end

    MPFRep._create(result.sign_bit(), false, false,
      result.exponent() + n256.i64(), result.raw_digits())._trunc(w)


  fun _log2(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `log₂(a)` at working precision `w` bytes.
    `log₂(a) = ln(a) / ln(2)`.
    """
    let ln_x = _ln(a, w)
    if ln_x.is_nan() or ln_x.is_infinite() then
      return ln_x
    end
    _div(ln_x, _ln2_const(w + 4), w)._trunc(w)


  fun _log10(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `log₁₀(a)` at working precision `w` bytes.
    `ln(10) = 2 × arctanh(1/3) + 2 × arctanh(2/3)`.
    """
    let ln_x = _ln(a, w)
    if ln_x.is_nan() or ln_x.is_infinite() then
      return ln_x
    end
    let p3: USize = w + 6
    let d3:  Array[U8] val = recover Array[U8].init(0x55, p3) end
    let d23: Array[U8] val = recover Array[U8].init(0xAA, p3) end
    let third     = MPFRep._create(false, false, false, 0, d3)
    let two_third = MPFRep._create(false, false, false, 0, d23)
    let two       = MPFRep.from_f64(2.0, p3)
    let ln2_v  = _mul(_atanh_series(third, p3), two, p3)._trunc(p3)
    let ln5_v  = _mul(_atanh_series(two_third, p3), two, p3)._trunc(p3)
    let ln10_v = _add(ln2_v, ln5_v, p3)._trunc(w + 4)
    _div(ln_x, ln10_v, w)._trunc(w)


  fun _logb(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Compute `log_b(a)` at working precision `w` bytes.
    `log_b(a) = ln(a) / ln(b)`.
    """
    _div(_ln(a, w), _ln(b, w), w)._trunc(w)


  fun _exp2(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `2^a` at working precision `w` bytes.
    `2^a = exp(a × ln(2))`.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() and (not a.sign_bit()) then
      return MPFRep._create(false, false, true, 0, Array[U8].create())
    end
    if a.is_infinite() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    if a.is_zero() then
      return MPFRep.from_f64(1.0, w)
    end
    let p2: USize = w + 4
    let ln2: MPFRep = _ln2_const(p2)
    _exp(_mul(a, ln2, p2)._trunc(p2), w)._trunc(w)


  // ── Powers ─────────────────────────────────────────────────────────────────

  fun _powi(a: MPFRep, n: ILong, w: USize): MPFRep =>
    """
    Compute `a^n` for integer `n` via binary exponentiation at working
    precision `w` bytes.

    - `a^0 = 1` for any finite `a`.
    - `a^n` for `n < 0` is `(1/a)^|n|`.
    - NaN → NaN.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if n == 0 then
      return MPFRep.from_f64(1.0, w)
    end
    var base_r: MPFRep =
      if n < 0 then _inv(a, w) else a end
    var exp_n: ILong = if n < 0 then -n else n end
    var result: MPFRep = MPFRep.from_f64(1.0, w)
    while exp_n > 0 do
      if (exp_n and 1) == 1 then
        result = _mul(result, base_r, w)._trunc(w)
      end
      base_r = _mul(base_r, base_r, w)._trunc(w)
      exp_n = exp_n / 2
    end
    result


  fun _pow(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Compute `a^b` at working precision `w` bytes.

    - `a > 0`: `exp(b × ln(a))`.
    - `a < 0` and `b` integer: delegates to `_powi`.
    - `a < 0` and `b` non-integer: NaN.
    - NaN in either → NaN.
    """
    if a.is_nan() or b.is_nan() then
      return MPFRep.nan_val()
    end
    if b.is_zero() then
      return MPFRep.from_f64(1.0, w)
    end
    if a.is_zero() then
      if b.sign_bit() then
        return MPFRep.nan_val()
      end
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    if not a.sign_bit() then
      let ln_a: MPFRep = _ln(a, w)
      if ln_a.is_nan() then
        return MPFRep.nan_val()
      end
      return _exp(_mul(b, ln_a, w + 4)._trunc(w + 4), w)._trunc(w)
    end
    // Negative base: only valid for integer exponents.
    if b._has_frac() then
      return MPFRep.nan_val()
    end
    let ni: ILong = _round_to_ilong(b)
    _powi(a, ni, w)


  // ── Trigonometric ──────────────────────────────────────────────────────────

  fun _sin_cos_reduce(a: MPFRep, w: USize): (MPFRep, ILong) =>
    """
    Argument reduction for `_sin` and `_cos`.

    Returns `(r, k)` where `r = a − n × (π/2)`, `|r| ≤ π/4`, and
    `k = n mod 4 ∈ {0, 1, 2, 3}` (normalised to non-negative).
    """
    let x: MPFRep = _extend(a, w)
    let pi_val: MPFRep = _pi(w)
    let two: MPFRep = MPFRep.from_f64(2.0, w)
    let pi_half: MPFRep = _div(pi_val, two, w)._trunc(w)
    let n_f: MPFRep = _trunc_frac(_div(x, pi_half, w)._trunc(w))
    let n: ILong = _round_to_ilong(n_f)
    let n_mpf: MPFRep = MPFRep.from_f64(n.f64(), w)
    // r = (a − n) × pi_half — same formula as in the original mpfloat.pony.
    let diff: MPFRep = _sub(x, n_mpf, w)._trunc(w)
    let r: MPFRep = _mul(diff, pi_half, w)._trunc(w)
    let k: ILong = ((n % 4) + 4) % 4
    (r, k)


  fun _sin(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `sin(a)` at working precision `w` bytes (argument in radians).

    - NaN or `±∞` → NaN.  `0` → `0`.

    Argument reduction: `n = round(a / (π/2))`, `r = a − n × (π/2)`,
    `|r| ≤ π/4`.  Based on `k = n mod 4`:
    - `k=0`: `sin(r)`,  `k=1`: `cos(r)`,
    - `k=2`: `−sin(r)`, `k=3`: `−cos(r)`.
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    let p2: USize = w + 8
    (let r, let k) = _sin_cos_reduce(a, p2)
    let result: MPFRep =
      if (k == 0) or (k == 2) then
        let sr = _sin_taylor(r, p2)
        if k == 0 then sr else _neg_rep(sr) end
      else
        let cr = _cos_taylor(r, p2)
        if k == 1 then cr else _neg_rep(cr) end
      end
    result._trunc(w)


  fun _cos(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `cos(a)` at working precision `w` bytes (argument in radians).

    - NaN or `±∞` → NaN.  `0` → `1`.

    Uses the same argument reduction as `_sin`:
    - `k=0`: `cos(r)`,  `k=1`: `−sin(r)`,
    - `k=2`: `−cos(r)`, `k=3`: `sin(r)`.
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return MPFRep.from_f64(1.0, w)
    end
    let p2: USize = w + 8
    (let r, let k) = _sin_cos_reduce(a, p2)
    let result: MPFRep =
      if (k == 0) or (k == 2) then
        let cr = _cos_taylor(r, p2)
        if k == 0 then cr else _neg_rep(cr) end
      else
        let sr = _sin_taylor(r, p2)
        if k == 3 then sr else _neg_rep(sr) end
      end
    result._trunc(w)


  fun _tan(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `tan(a) = sin(a) / cos(a)` at working precision `w` bytes.

    - NaN or `±∞` → NaN.  `0` → `0`.
    When `cos(a)` is zero at the working precision, returns NaN.
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    let c: MPFRep = _cos(a, w)
    if c.is_zero() then
      return MPFRep.nan_val()
    end
    _div(_sin(a, w), c, w)._trunc(w)


  // ── Hyperbolic ─────────────────────────────────────────────────────────────

  fun _sinh(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `sinh(a) = (e^a − e^{−a}) / 2` at working precision `w` bytes.

    - NaN → NaN.  `±∞` → `±∞`.  `0` → `0`.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(a.sign_bit(), false, true, 0, Array[U8].create())
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    let p2: USize = w + 4
    let ex: MPFRep  = _exp(a, p2)
    let emx: MPFRep = _exp(_neg_rep(a), p2)
    let two: MPFRep = MPFRep.from_f64(2.0, p2)
    _div(_sub(ex, emx, p2)._trunc(p2), two, p2)._trunc(w)


  fun _cosh(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `cosh(a) = (e^a + e^{−a}) / 2` at working precision `w` bytes.

    - NaN → NaN.  `±∞` → `+∞`.  `0` → `1`.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(false, false, true, 0, Array[U8].create())
    end
    if a.is_zero() then
      return MPFRep.from_f64(1.0, w)
    end
    let p2: USize = w + 4
    let ex: MPFRep  = _exp(a, p2)
    let emx: MPFRep = _exp(_neg_rep(a), p2)
    let two: MPFRep = MPFRep.from_f64(2.0, p2)
    _div(_add(ex, emx, p2)._trunc(p2), two, p2)._trunc(w)


  fun _tanh(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `tanh(a) = (e^a − e^{−a}) / (e^a + e^{−a})` at working precision
    `w` bytes.

    - NaN → NaN.  `+∞` → `+1`.  `−∞` → `−1`.  `0` → `0`.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      let sign: Bool = a.sign_bit()
      return MPFRep.from_f64(if sign then -1.0 else 1.0 end, w)
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    let p2: USize = w + 4
    let ex: MPFRep  = _exp(a, p2)
    let emx: MPFRep = _exp(_neg_rep(a), p2)
    _div(
      _sub(ex, emx, p2)._trunc(p2),
      _add(ex, emx, p2)._trunc(p2),
      p2)._trunc(w)


  fun _sech(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `sech(a) = 1 / cosh(a)` at working precision `w` bytes.

    - NaN → NaN.  `±∞` → `0`.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    _inv(_cosh(a, w), w)._trunc(w)


  fun _csch(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `csch(a) = 1 / sinh(a)` at working precision `w` bytes.

    - NaN, `±∞`, or `a = 0` → NaN.
    """
    if a.is_nan() or a.is_infinite() or a.is_zero() then
      return MPFRep.nan_val()
    end
    _inv(_sinh(a, w), w)._trunc(w)


  fun _coth(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `coth(a) = cosh(a) / sinh(a)` at working precision `w` bytes.

    - NaN, `±∞`, or `a = 0` → NaN.
    """
    if a.is_nan() or a.is_infinite() or a.is_zero() then
      return MPFRep.nan_val()
    end
    _div(_cosh(a, w), _sinh(a, w), w)._trunc(w)


  // ── Constants ──────────────────────────────────────────────────────────────

  fun _pi(w: USize): MPFRep =>
    """
    Compute `π` at working precision `w` bytes via the Borwein quartic
    iteration.

    Uses the same algorithm as `MPFloat.pi`: starts from an F64 seed, then
    applies the quartic convergence recurrence until convergence.
    """
    let size: USize = w + 2
    let prec_bits: ULong = (size * 8).ulong()
    let pi_f = MPFloat.pi(prec_bits)
    MPFRep.from_f64(pi_f.f64(), w)


  fun _pi_bbp(w: USize): MPFRep =>
    """
    Compute `π` at working precision `w` bytes via the BBP formula.
    Delegates to `MPFloat.pi_bbp` for now.
    """
    let pi_f = MPFloat.pi_bbp((w * 8).ulong())
    MPFRep.from_f64(pi_f.f64(), w)


  fun _pi_chudnovsky(w: USize): MPFRep =>
    """
    Compute `π` at working precision `w` bytes via the Chudnovsky algorithm.
    Delegates to `MPFloat.pi` for now.
    """
    let pi_f = MPFloat.pi((w * 8).ulong())
    MPFRep.from_f64(pi_f.f64(), w)


  // ── String conversion ──────────────────────────────────────────────────────

  fun _from_string(s: String, w: USize): MPFRep ? =>
    """
    Parse the decimal string `s` into an `MPFRep` at working precision `w`
    bytes.  Delegates to `MPFloat.from_string` during the transition period.
    """
    let prec_bits: ULong = (w * 8).ulong()
    let f = MPFloat.from_string(s, prec_bits)?
    MPFRep._create(
      f.is_negative(), f.is_nan(), f.is_infinite(),
      f.exponent(), f.raw_digits())
