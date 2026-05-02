// Multi-precision float algorithms

use "../assertx"
use "../formatx"

use "debug"
use "collections"


primitive _MPFAlgo
  """
  All arithmetic and transcendental algorithms for `MPFRep` values.

  ## Two layers of methods

  **Public (context-aware)**: `add`, `sub`, `mul`, `div`, `inv`, `sqrt`,
  `divrem`, `fld`, `rem`, `mod`, `ln`, `log`, `log2`, `log10`, `logb`,
  `exp`, `exp2`, `powi`, `pow`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`,
  `sech`, `csch`, `coth`, `pi`, `pi_bbp`, `pi_chudnovsky`, `from_string`.

  Each takes `ctx: MPFContext` as its first argument. It handles **all**
  special values (NaN, ±∞, ±0) with early returns before touching `_*`
  methods, runs the computation at working precision `ctx.working_bytes(op)`,
  then rounds to output precision via `ctx._round_to`. `MPFloat` calls these
  directly — there is no dispatch layer in `MPFContext`.

  **Internal (precision-explicit)**: `_add`, `_sub`, `_mul`, `_inv`, `_sqrt`,
  `_divrem`, `_ln`, `_exp`, etc., and series helpers `_atanh_series`,
  `_exp_taylor`, `_sin_taylor`, `_cos_taylor`, `_ln2_const`.

  Each takes `w: USize` (working precision in bytes), performs the computation,
  and returns an `MPFRep` at that precision (unrounded). Internal loop bodies
  call `MPFRep._trunc(w)` to bound intermediate growth; this is distinct from
  the final output rounding applied by the public layer.

  ## Invariant: special-value handling belongs exclusively to public methods

  Internal `_*` methods are pure computation kernels. They assume their
  inputs are **finite, non-zero, and non-NaN**. They make no guarantees for
  any other input and must not contain NaN/inf/zero guard clauses.

  When adding a new operation:
  1. Implement the algorithm in a `_foo(a, w)` method with no special-value
     checks. Document the finite/non-zero precondition in its docstring.
  2. Add a public `foo(ctx, a)` method that dispatches all special values
     (NaN, ±∞, ±0, sign edge cases) before calling `_foo`.
  3. Never add guards to `_foo` as a shortcut — that duplicates logic and
     silently breaks the contract for callers that rely on the kernel being
     a pure computation.
  """

  // ── Addition / Subtraction ─────────────────────────────────────────────────

  fun _add(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Add `a + b` at working precision `w` bytes.

    Assumes both operands are finite and non-zero. Special-value handling
    (NaN, ±∞, zero) belongs to the public `add` method.

    When operand signs agree, magnitudes are added (`_add_mag`). When they
    differ, the smaller magnitude is subtracted from the larger and the result
    takes the sign of the larger. The result is truncated to `w` bytes.
    """
    let res: MPFRep =
      if a.sign_bit() == b.sign_bit() then
        // Same sign: add magnitudes.
        if a.exponent() >= b.exponent() then
          a._add_mag(b, a.sign_bit())
        else
          b._add_mag(a, a.sign_bit())
        end
      else
        // Different signs: subtract smaller magnitude from larger.
        match a._cmp_mag(b)
        | Greater => a._sub_mag(b, a.sign_bit())
        | Less    => b._sub_mag(a, b.sign_bit())
        else
          MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
        end
      end
    res._trunc(w)


  fun _sub(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Subtract `a − b` at working precision `w` bytes. Delegates to `_add`
    after negating `b`, inheriting all sign and special-value handling.
    """
    _add(a, _neg_rep(b), w)


  // ── Multiplication ─────────────────────────────────────────────────────────

  fun _mul(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Multiply `a × b` at working precision `w` bytes via FFT convolution.

    The result sign is `a.sign_bit() XOR b.sign_bit()`. The raw convolution
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
        let base_f: F64 = 256.0
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
        Fail(Format("[_MPFAlgo._mul] Index out of bounds at i={}", err_i))
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
    with exponent 1, so `G ∈ [1, 256)`. The iteration converges to `y = 1/G`;
    the final result exponent is `(y._exponent + 1) − a.exponent()`.

    After the main loop converges at `w + 1` bytes (faithful, ≤ 1 ULP), one
    extra Newton step at `w + 3` bytes halves the error to < 0.5 ULP relative
    to `w`, so the caller's `_round_to(result, p)` is always unambiguous.

    Sign propagation: the sign of the result equals `a.sign_bit()`.
    """
    let size: USize = w + 1
    let size2: USize = w + 3  // refinement precision
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
      let gy = _mul(g, res, size)
      let delta = _sub(two, gy, size)
      let new_res = _mul(res, delta, size)

      let converged = _digits_eq(new_res, res)
      res = new_res
      if converged then
        break
      end
    end

    // One extra refinement step at size2 bytes to push error below 0.5 ULP.
    let two2 = MPFRep.from_f64(2.0, size2)
    let gy2 = _mul(g, res, size2)
    let delta2 = _sub(two2, gy2, size2)
    res = _mul(res, delta2, size2)

    // 1/|a| = res × 256^{1−e}. Adjust sign and exponent.
    MPFRep._create(
      a.sign_bit(), false, false,
      (res.exponent() + 1) - a.exponent(),
      res.raw_digits())


  fun _div(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Divide `a / b` at working precision `w` bytes. Computes `a × inv(b)`.

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
    `H ∈ [1, 256)`. For even `a.exponent()` use exponent 2 so
    `H ∈ [256, 65536)`. With this split `a.exponent() − h_exp` is always
    even, so the result exponent is always an integer.

    After the main loop converges at `w + 1` bytes (faithful, ≤ 1 ULP), one
    extra refinement step at `w + 3` bytes halves the error to < 0.5 ULP
    relative to `w`, so the caller's `_round_to(result, p)` is unambiguous.
    """
    let size: USize = w + 1
    let size2: USize = w + 3  // refinement precision
    let ad: Array[U8] val = a.raw_digits()
    let prec: USize = a._size()
    let a_exp: I64 = a.exponent()

    // Choose h_exp so that a_exp − h_exp is always even.
    let h_exp: I64 = if (a_exp and 1) != 0 then 1 else 2 end
    let h = MPFRep._create(false, false, false, h_exp, ad)

    // F64 initial estimate.
    let ng: USize = prec.min(4)
    var fg: F64 = 0.0
    let base_f: F64 = 256.0
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
      let y2 = _mul(res, res, size)
      let hy2 = _mul(h, y2, size)
      let delta = _sub(three, hy2, size)
      let new_full = _mul(res, delta, size)
      (let halved, _) = new_full._short_div(2)
      let nr = halved._trunc(size)

      let converged = _digits_eq(nr, res)
      res = nr
      if converged then
        break
      end
    end

    // One extra refinement step at size2 bytes to push error below 0.5 ULP.
    let three2 = MPFRep.from_f64(3.0, size2)
    let y2r = _mul(res, res, size2)
    let hy2r = _mul(h, y2r, size2)
    let delta2 = _sub(three2, hy2r, size2)
    let new_full2 = _mul(res, delta2, size2)
    (let halved2, _) = new_full2._short_div(2)
    res = halved2._trunc(size2)

    // √H = H × (1/√H); exponent adjusted for the parity split.
    let sqrt_h = _mul(h, res, size2)
    let result_exp: I64 = sqrt_h.exponent() + ((a_exp - h_exp) / 2)
    MPFRep._create(false, false, false, result_exp, sqrt_h.raw_digits())


  // ── Division with Remainder ────────────────────────────────────────────────

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
    When `|r| ≥ |b|`, increment `|q|` by 1 and recompute `r`. At most one
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

    var q: MPFRep = _div(a, b, w)._trunc_frac()
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

    The result is the largest integer `q` such that `q × b ≤ a`. When the
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

    The result has the same sign as `b`. Special cases are inherited from
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

    The result has the same sign as `a`. Special cases are inherited from
    `_divrem`.
    """
    _divrem(a, b, w)._2


  // ── Internal convergence helpers ──────────────────────────────────────────

  fun _digits_eq(a: MPFRep, b: MPFRep): Bool =>
    """
    Return `true` when `a` and `b` have the same exponent and the same leading
    digits (compared byte by byte up to the shorter array). Used for
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
    try
      for i in Range(0, n) do
        if ad(i)? != bd(i)? then
          return false
        end
      end
    end
    true


  fun _neg_rep(a: MPFRep): MPFRep =>
    """
    Return a copy of `a` with its sign flipped. Equivalent to arithmetic
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

    The caller must ensure `|x| < 1` for convergence. The loop terminates
    when the current addend no longer changes the running sum at precision `w`.

    Uses `MPFRep._short_div` (exact U8 division) for denominators up to 255;
    falls back to `_div` for larger denominators (only for very high `w`).
    """
    let p2: USize = w + 2
    let x2: MPFRep = _mul(x, x, p2)
    var term: MPFRep = x
    var sum:  MPFRep = term
    var k: USize = 3
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term = _mul(term, x2, p2)
      let addend: MPFRep = if k <= 255 then
          (let q, _) = term._short_div(k.u8())
          q
        else
          _div(term, MPFRep.from_f64(k.f64(), p2), p2)
        end
      let new_sum: MPFRep = _add(sum, addend, p2)
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
    factorials. Works best for `|r| ≤ ln(2)/2 ≈ 0.347`.
    """
    let p2: USize = w + 2
    let one = _one(p2)
    var term: MPFRep = one
    var sum:  MPFRep = one
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term = _mul(term, r, p2)
      let divided: MPFRep = if k <= 255 then
          (let q, _) = term._short_div(k.u8())
          q
        else
          _div(term, MPFRep.from_f64(k.f64(), p2), p2)
        end
      term = divided
      let new_sum: MPFRep = _add(sum, divided, p2)
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
    two factorial positions at once. Works best for `|r| ≤ π/4 ≈ 0.785`.
    """
    let p2: USize = w + 2
    let neg_r2: MPFRep = _neg_rep(_mul(r, r, p2))
    var term: MPFRep = r
    var sum:  MPFRep = term
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term = _mul(term, neg_r2, p2)
      let d1: USize = 2 * k
      let d2: USize = (2 * k) + 1
      term = if d1 <= 255 then
          (let q, _) = term._short_div(d1.u8())
          q
        else
          _div(term, MPFRep.from_f64(d1.f64(), p2), p2)
        end
      term =
        if d2 <= 255 then
          (let q, _) = term._short_div(d2.u8())
          q
        else
          _div(term, MPFRep.from_f64(d2.f64(), p2), p2)
        end
      let new_sum: MPFRep = _add(sum, term, p2)
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
    two factorial positions at once. Works best for `|r| ≤ π/4 ≈ 0.785`.
    """
    let p2: USize = w + 2
    let neg_r2 = _neg_rep(_mul(r, r, p2))
    let one = _one(p2)
    var term: MPFRep = one
    var sum:  MPFRep = term
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term = _mul(term, neg_r2, p2)
      let d1: USize = (2 * k) - 1
      let d2: USize = 2 * k
      term = if d1 <= 255 then
          (let q, _) = term._short_div(d1.u8())
          q
        else
          _div(term, MPFRep.from_f64(d1.f64(), p2), p2)
        end
      term = if d2 <= 255 then
          (let q, _) = term._short_div(d2.u8())
          q
        else
          _div(term, MPFRep.from_f64(d2.f64(), p2), p2)
        end
      let new_sum: MPFRep = _add(sum, term, p2)
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
    Extract the nearest integer value of `a` as an `ILong` by first truncating
    the fractional part and then converting via `MPInt`.

    Returns 0 when `a` is NaN, ±∞, or zero. Used internally in argument
    reduction steps for `_exp`, `_sin`, and `_cos`.

    Only the first `exponent()` bytes of the digit array are the integer part;
    the rest are fractional and must be excluded before passing to `MPInt`.
    """
    if (not a.is_finite()) or a.is_zero() then
      return 0
    end
    let rounded: MPFRep = a._trunc_frac()
    let e: USize = rounded.exponent().usize().min(rounded._size())
    if e == 0 then
      return 0
    end
    let int_bytes: Array[U8] val = recover
      let src = rounded.raw_digits()
      let dst = Array[U8].init(0, e)
      src.copy_to(dst, 0, 0, e)
      dst
    end
    MPInt.from_bytes_be(rounded.sign_bit(), int_bytes).ilong()


  fun _extend(a: MPFRep, w: USize): MPFRep =>
    """
    Return `a` extended (or truncated) to `w` bytes by zero-padding the tail
    or calling `_trunc(w)`. The sign, exponent, and leading digits are
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

    Assumes `a` is finite, positive, and non-zero. Special-value handling
    belongs to the public `ln` method.

    Algorithm: two-level argument reduction then `_atanh_series`:
    1. `a = 0.d₀d₁… × 256^e`. Set `k = e − 1`, `m = 0.d₀d₁… × 256 ∈ [1, 256)`.
       `ln(a) = k × ln(256) + ln(m)`.
    2. Find `n = ⌊log₂(d₀)⌋ ∈ {0,…,7}`, `u = m / 2^n ∈ [1, 2)`.
       `ln(m) = n × ln(2) + ln(u)`.
    3. `t = (u−1)/(u+1) ∈ [0, 1/3)`. `ln(u) = 2 × arctanh(t)`.
    Combined: `ln(a) = (8k + n) × ln(2) + 2 × arctanh(t)`.
    """
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

    let one = _one(p2)
    let two = MPFRep.from_f64(2.0, p2)
    let t: MPFRep = _div(_sub(u, one, p2), _add(u, one, p2), p2)
    let ln_u: MPFRep = _mul(_atanh_series(t, p2), two, p2)

    let ln2_factor: I64 = (8 * k) + n_inner
    let correction: MPFRep = if ln2_factor == 0 then
        MPFRep._create(false, false, false, 0, Array[U8].init(0, p2))
      else
        let fac = MPFRep.from_f64(ln2_factor.f64(), p2)
        _mul(fac, ln2, p2)
      end
    _add(correction, ln_u, p2)._trunc(w)


  fun _exp(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `e^a` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to the
    public `exp` method.

    Algorithm: two-level argument reduction then `_exp_taylor`:
    1. `a = n₂₅₆ × ln(256) + r₁`, `n₂₅₆ = round(a / ln(256))`.
       Scale by `256^n₂₅₆` (adjust exponent).
    2. `r₁ = n₂ × ln(2) + r`, `|r| ≤ ln(2)/2`.
       Scale by `2^n₂` (exact).
    3. Compute `e^r` via `_exp_taylor`.
    """
    let p2: USize = w + 8

    let ln2: MPFRep = _ln2_const(p2)
    let eight = MPFRep._create(false, false, false, 1, recover [U8(8)] end)
    let ln256: MPFRep = _mul(ln2, eight, p2)

    let x: MPFRep = _extend(a, p2)

    let n256_f: MPFRep = _div(x, ln256, p2)._trunc_frac()
    let n256: ILong = _round_to_ilong(n256_f)
    let n256_mpf = MPFRep.from_f64(n256.f64(), p2)
    let r1: MPFRep = _sub(x, _mul(n256_mpf, ln256, p2), p2)

    let n2_f: MPFRep = _div(r1, ln2, p2)._trunc_frac()
    let n2: ILong = _round_to_ilong(n2_f)
    let n2_mpf = MPFRep.from_f64(n2.f64(), p2)
    // r = r1 − n2 × ln2.
    let r: MPFRep = _sub(r1, _mul(n2_mpf, ln2, p2), p2)

    var result: MPFRep = _exp_taylor(r, p2)

    if n2 > 0 then
      let pow2 = MPFRep._create(false, false, false, 1, recover [U8(1).shl(n2.u8())] end)
      result = _mul(result, pow2, p2)
    elseif n2 < 0 then
      let pow2 = MPFRep._create(false, false, false, 0, recover [U8(128).shr((-n2 - 1).u8())] end)
      result = _mul(result, pow2, p2)
    end

    MPFRep._create(result.sign_bit(), false, false,
      result.exponent() + n256.i64(), result.raw_digits())._trunc(w)


  fun _log2(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `log₂(a)` at working precision `w` bytes.
    `log₂(a) = ln(a) / ln(2)`.
    Assumes `a` is finite, positive, and non-zero.
    """
    _div(_ln(a, w), _ln2_const(w + 4), w)


  fun _log10(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `log₁₀(a)` at working precision `w` bytes.
    `ln(10) = 2 × arctanh(1/3) + 2 × arctanh(2/3)`.
    Assumes `a` is finite, positive, and non-zero.
    """
    let ln_x = _ln(a, w)
    let p3: USize = w + 6
    let d3:  Array[U8] val = recover Array[U8].init(0x55, p3) end
    let d23: Array[U8] val = recover Array[U8].init(0xAA, p3) end
    let third     = MPFRep._create(false, false, false, 0, d3)
    let two_third = MPFRep._create(false, false, false, 0, d23)
    let two       = MPFRep.from_f64(2.0, p3)
    let ln2_v  = _mul(_atanh_series(third, p3), two, p3)
    let ln5_v  = _mul(_atanh_series(two_third, p3), two, p3)
    let ln10_v = _add(ln2_v, ln5_v, p3)._trunc(w + 4)
    _div(ln_x, ln10_v, w)


  fun _logb(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Compute `log_b(a)` at working precision `w` bytes.
    `log_b(a) = ln(a) / ln(b)`.
    """
    _div(_ln(a, w), _ln(b, w), w)


  fun _exp2(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `2^a` at working precision `w` bytes.
    `2^a = exp(a × ln(2))`.
    Assumes `a` is finite and non-zero.
    """
    let p2: USize = w + 4
    let ln2: MPFRep = _ln2_const(p2)
    _exp(_mul(a, ln2, p2), w)


  // ── Powers ─────────────────────────────────────────────────────────────────

  fun _powi(a: MPFRep, n: ILong, w: USize): MPFRep =>
    """
    Compute `a^n` for integer `n` via binary exponentiation at working
    precision `w` bytes.

    Assumes `a` is finite and non-NaN, and `n ≠ 0`. Special-value handling
    belongs to the public `powi` method.
    `a^n` for `n < 0` is `(1/a)^|n|`.
    """
    var base_r: MPFRep = if n < 0 then _inv(a, w) else a end
    var exp_n: ILong = if n < 0 then -n else n end
    var result: MPFRep = _one(w)
    while exp_n > 0 do
      if (exp_n and 1) == 1 then
        result = _mul(result, base_r, w)
      end
      base_r = _mul(base_r, base_r, w)
      exp_n = exp_n / 2
    end
    result


  fun _pow(a: MPFRep, b: MPFRep, w: USize): MPFRep =>
    """
    Compute `a^b` at working precision `w` bytes.

    Assumes both operands are finite and non-NaN, `b ≠ 0`, `a ≠ 0`, and
    that negative-base/non-integer-exponent has already been rejected.
    Special-value handling belongs to the public `pow` method.
    `a > 0`: `exp(b × ln(a))`. `a < 0` with integer `b`: delegates to `_powi`.
    """
    if not a.sign_bit() then
      return _exp(_mul(b, _ln(a, w + 4), w + 4), w)
    end
    // Negative base with integer exponent (non-integer already rejected by caller).
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
    let half: MPFRep = MPFRep.from_f64(0.5, w)
    let pi_half: MPFRep = _div(pi_val, two, w)
    // round-to-nearest: n = floor(x/pi_half + 0.5)
    let ratio: MPFRep = _div(x, pi_half, w)
    let n_f: MPFRep = _add(ratio, half, w)._trunc_frac()
    let n: ILong = _round_to_ilong(n_f)
    let n_mpf: MPFRep = MPFRep.from_f64(n.f64(), w)
    // r = x − n × pi_half  (NOT (x−n) × pi_half — that is a bug).
    let r: MPFRep = _sub(x, _mul(n_mpf, pi_half, w), w)
    let k: ILong = ((n % 4) + 4) % 4
    (r, k)


  fun _sin(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `sin(a)` at working precision `w` bytes (argument in radians).

    Assumes `a` is finite and non-zero. Special-value handling belongs to `sin`.

    Argument reduction: `n = round(a / (π/2))`, `r = a − n × (π/2)`,
    `|r| ≤ π/4`. Based on `k = n mod 4`:
    - `k=0`: `sin(r)`,  `k=1`: `cos(r)`,
    - `k=2`: `−sin(r)`, `k=3`: `−cos(r)`.
    """
    let p2: USize = w + 8
    (let r, let k) = _sin_cos_reduce(a, p2)
    let result: MPFRep = if (k == 0) or (k == 2) then
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

    Assumes `a` is finite and non-zero. Special-value handling belongs to `cos`.

    Argument reduction: `n = round(a / (π/2))`, `r = a − n × (π/2)`,
    `|r| ≤ π/4`. Based on `k = n mod 4`:
    - `k=0`: `cos(r)`,  `k=1`: `−sin(r)`,
    - `k=2`: `−cos(r)`, `k=3`: `sin(r)`.
    """
    let p2: USize = w + 8
    (let r, let k) = _sin_cos_reduce(a, p2)
    let result: MPFRep = if (k == 0) or (k == 2) then
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

    Assumes `a` is finite and non-zero. Special-value handling belongs to `tan`.
    When `cos(a)` is zero at the working precision, returns NaN.
    """
    let c: MPFRep = _cos(a, w)
    if c.is_zero() then
      return MPFRep.nan_val()
    end
    _div(_sin(a, w), c, w)


  // ── Hyperbolic ─────────────────────────────────────────────────────────────

  fun _sinh(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `sinh(a) = (e^a − e^{−a}) / 2` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `sinh`.
    """
    let p2: USize = w + 4
    let ex: MPFRep  = _exp(a, p2)
    let emx: MPFRep = _exp(_neg_rep(a), p2)
    let two: MPFRep = MPFRep.from_f64(2.0, p2)
    _div(_sub(ex, emx, p2), two, p2)._trunc(w)


  fun _cosh(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `cosh(a) = (e^a + e^{−a}) / 2` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `cosh`.
    """
    let p2: USize = w + 4
    let ex: MPFRep  = _exp(a, p2)
    let emx: MPFRep = _exp(_neg_rep(a), p2)
    let two: MPFRep = MPFRep.from_f64(2.0, p2)
    _div(_add(ex, emx, p2), two, p2)._trunc(w)


  fun _tanh(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `tanh(a) = (e^a − e^{−a}) / (e^a + e^{−a})` at working precision
    `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `tanh`.
    """
    let p2: USize = w + 4
    let ex: MPFRep  = _exp(a, p2)
    let emx: MPFRep = _exp(_neg_rep(a), p2)
    _div(_sub(ex, emx, p2), _add(ex, emx, p2), p2)._trunc(w)


  fun _sech(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `sech(a) = 1 / cosh(a)` at working precision `w` bytes.

    Assumes `a` is finite. Special-value handling belongs to `sech`.
    """
    _inv(_cosh(a, w), w)


  fun _csch(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `csch(a) = 1 / sinh(a)` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `csch`.
    """
    _inv(_sinh(a, w), w)


  fun _coth(a: MPFRep, w: USize): MPFRep =>
    """
    Compute `coth(a) = cosh(a) / sinh(a)` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `coth`.
    """
    _div(_cosh(a, w), _sinh(a, w), w)


  // ── Constants ──────────────────────────────────────────────────────────────

  fun _pi(w: USize): MPFRep =>
    """
    Compute `π` at working precision `w` bytes via the Machin formula
    `π = 16·arctan(1/5) − 4·arctan(1/239)`.

    Uses Taylor `arctan(x) = x − x³/3 + x⁵/5 − ...`
    """
    let p: USize = w + 4
    let one: MPFRep = _one(p)
    let five: MPFRep = MPFRep.from_f64(5.0, p)
    let two39: MPFRep = MPFRep.from_f64(239.0, p)
    let x5: MPFRep = _inv(five, p)
    let x239: MPFRep = _inv(two39, p)

    // arctan(1/5) via Taylor series
    let neg_x52: MPFRep = _neg_rep(_mul(x5, x5, p))
    var pow5: MPFRep = x5
    var atan5: MPFRep = x5
    var k5: USize = 1
    var iters5: USize = 0
    let max_iters: USize = p * 30
    while iters5 < max_iters do
      pow5 = _mul(pow5, neg_x52, p)
      let d5: USize = (2 * k5) + 1
      let term5: MPFRep = if d5 <= 255 then
          (let q5, _) = pow5._short_div(d5.u8())
          q5
        else
          _div(pow5, MPFRep.from_f64(d5.f64(), p), p)
        end
      let new_a5: MPFRep = _add(atan5, term5, p)
      if _digits_eq(new_a5, atan5) then break end
      atan5 = new_a5
      k5 = k5 + 1
      iters5 = iters5 + 1
    end

    // arctan(1/239) via Taylor series
    let neg_x2392: MPFRep = _neg_rep(_mul(x239, x239, p))
    var pow239: MPFRep = x239
    var atan239: MPFRep = x239
    var k239: USize = 1
    var iters239: USize = 0
    while iters239 < max_iters do
      pow239 = _mul(pow239, neg_x2392, p)
      let d239: USize = (2 * k239) + 1
      let term239: MPFRep = if d239 <= 255 then
          (let q239, _) = pow239._short_div(d239.u8())
          q239
        else
          _div(pow239, MPFRep.from_f64(d239.f64(), p), p)
        end
      let new_a239: MPFRep = _add(atan239, term239, p)
      if _digits_eq(new_a239, atan239) then break end
      atan239 = new_a239
      k239 = k239 + 1
      iters239 = iters239 + 1
    end

    // π = 16·arctan(1/5) − 4·arctan(1/239)
    let sixteen: MPFRep = MPFRep.from_f64(16.0, p)
    let four: MPFRep = MPFRep.from_f64(4.0, p)
    _sub(_mul(sixteen, atan5, p), _mul(four, atan239, p), p)._trunc(w)


  fun _pi_bbp(w: USize): MPFRep =>
    """
    Compute `π` at working precision `w` bytes via the BBP formula:
    `π = Σ_{k=0}^∞ (1/16^k) × [4/(8k+1) − 2/(8k+4) − 1/(8k+5) − 1/(8k+6)]`.

    The `1/16^k` factor is maintained as a running product to avoid `powi`.
    """
    let p: USize = w + 4
    let k_1: MPFRep = _one(p)
    let k_2: MPFRep = MPFRep.from_f64(2.0, p)
    let k_4: MPFRep = MPFRep.from_f64(4.0, p)
    let k_5: MPFRep = MPFRep.from_f64(5.0, p)
    let k_6: MPFRep = MPFRep.from_f64(6.0, p)
    let k_8: MPFRep = MPFRep.from_f64(8.0, p)
    let inv16: MPFRep = _inv(MPFRep.from_f64(16.0, p), p)
    var k: USize = 0
    var result: MPFRep = MPFRep._create(false, false, false, 0, Array[U8].init(0, p))
    var prev_res: MPFRep = result
    var pow16k: MPFRep = k_1
    repeat
      let t0: MPFRep = _mul(k_8, MPFRep.from_ulong(k.ulong(), p), p)
      let t1: MPFRep = _div(k_4, _add(t0, k_1, p), p)
      let t2: MPFRep = _div(k_2, _add(t0, k_4, p), p)
      let t3: MPFRep = _div(k_1, _add(t0, k_5, p), p)
      let t4: MPFRep = _div(k_1, _add(t0, k_6, p), p)
      let inner: MPFRep = _sub(_sub(_sub(t1, t2, p), t3, p), t4, p)
      let term: MPFRep = _mul(pow16k, inner, p)
      prev_res = result
      result = _add(result, term, p)
      pow16k = _mul(pow16k, inv16, p)
      k = k + 1
    until _digits_eq(result, prev_res) or (k > p) end
    result._trunc(w)


  fun _pi_chudnovsky(w: USize): MPFRep =>
    """
    Compute `π` at working precision `w` bytes via the Chudnovsky algorithm.

    TODO: Correct following bug,
    ⚠ Warning: gives only ~4 correct digits — the algorithm has a known bug in
    `640320^(3k+1.5)` computed via `pow` (transcendental path) losing precision.
    """
    let p: USize = w + 4
    let one: MPFRep = _one(p)
    let minus_one: MPFRep = _neg_rep(one)
    let k_3: MPFRep = MPFRep.from_f64(3.0, p)
    let k_1_5: MPFRep = _div(k_3, MPFRep.from_f64(2.0, p), p)
    let k_640320: MPFRep = MPFRep.from_ulong(640320, p)
    let k_13591409: MPFRep = MPFRep.from_ulong(13591409, p)
    let k_545140134: MPFRep = MPFRep.from_ulong(545140134, p)
    let k_12: MPFRep = MPFRep.from_f64(12.0, p)
    var result: MPFRep = MPFRep._create(false, false, false, 0, Array[U8].init(0, p))
    var prev_res: MPFRep = result
    var k: USize = 0
    repeat
      (let fact_k, let fact_3k, let fact_6k) = _fact_chudnovsky(k.ulong(), p)
      let sgn: MPFRep = if (k %% 2) == 0 then one else minus_one end
      let num: MPFRep = _add(k_13591409,
        _mul(k_545140134, MPFRep.from_ulong(k.ulong(), p), p), p)
      let exp_arg: MPFRep = _add(
        _mul(k_3, MPFRep.from_ulong(k.ulong(), p), p), k_1_5, p)
      let den_pow: MPFRep = _pow(k_640320, exp_arg, p)
      let fact_k3: MPFRep = _powi(fact_k, 3, p)
      let den: MPFRep = _mul(_mul(fact_3k, fact_k3, p), den_pow, p)
      let term: MPFRep = _mul(_mul(sgn, fact_6k, p), _div(num, den, p), p)
      prev_res = result
      result = _add(result, term, p)
      k = k + 1
    until _digits_eq(result, prev_res) or (k > p) end
    _inv(_mul(k_12, result, p), p)._trunc(w)


  fun _fact_chudnovsky(n: ULong, p: USize): (MPFRep, MPFRep, MPFRep) =>
    """
    Compute `n!`, `(3n)!`, `(6n)!` as `MPFRep` at precision `p` bytes for the
    Chudnovsky algorithm.
    """
    var fact_n: MPInt = MPInt.from[ULong](1)
    for i in Range[ULong](1, n + 1) do
      fact_n = fact_n * MPInt.from[ULong](i)
    end
    var fact_3n: MPInt = fact_n
    for i in Range[ULong](n + 1, (3 * n) + 1) do
      fact_3n = fact_3n * MPInt.from[ULong](i)
    end
    var fact_6n: MPInt = fact_3n
    for i in Range[ULong]((3 * n) + 1, (6 * n) + 1) do
      fact_6n = fact_6n * MPInt.from[ULong](i)
    end
    ( MPFRep.from_mpint(fact_n, p),
      MPFRep.from_mpint(fact_3n, p),
      MPFRep.from_mpint(fact_6n, p) )


  // ── String conversion ──────────────────────────────────────────────────────

  fun _from_string(s: String, w: USize): MPFRep ? =>
    """
    Parse the decimal string `s` into an `MPFRep` at working precision `w`
    bytes using the same multi-precision decimal algorithm as `MPFloat.from_string`.
    """
    let p_digits: USize = w
    let st: String = s.clone() .> strip()

    if st.size() == 0 then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, p_digits))
    end
    if (st == "nan") or (st == "NaN") or (st == "@NaN@") then
      return MPFRep.nan_val()
    end
    if (st == "+inf") or (st == "@Inf@") or (st == "inf") then
      return MPFRep.inf_val(true)
    end
    if (st == "-inf") or (st == "-@Inf@") then
      return MPFRep.inf_val(false)
    end
    if (st == "0") or (st == "0.0") or (st == "+0") or (st == "+0.0") then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, p_digits))
    end
    if (st == "-0") or (st == "-0.0") then
      return MPFRep._create(true, false, false, 0, Array[U8].init(0, p_digits))
    end

    let sig_limit: USize = ((p_digits.f64() * 2.41).usize() + 4).max(1)
    var pos: USize = 0
    let sz: USize = st.size()
    var fsign: Bool = false
    if st(pos)? == '-' then
      fsign = true ; pos = pos + 1
    elseif st(pos)? == '+' then
      pos = pos + 1
    end

    let ten_int: MPInt = MPInt.from[ILong](10)
    var n_int: MPInt = MPInt.from[ILong](0)
    var sig_count: USize = 0
    var int_exp: I64 = 0
    var frac_count: I64 = 0
    var has_digit: Bool = false
    var sep_seen: Bool = false

    while (pos < sz) and
          (((st(pos)? >= '0') and (st(pos)? <= '9')) or (st(pos)? == '_'))
    do
      if st(pos)? == '_' then
        if not sep_seen then
          sep_seen = true
          pos = pos + 1
          continue
        else
          Debug("[_MPFAlgo._from_string] consecutive '_' separators")
        end
      end
      sep_seen = false
      let d: U8 = st(pos)? - '0'
      has_digit = true
      if (d != 0) or (sig_count > 0) then
        if sig_count < sig_limit then
          n_int = (n_int * ten_int) + MPInt.from[ILong](d.ilong())
          sig_count = sig_count + 1
        else
          int_exp = int_exp + 1
        end
      end
      pos = pos + 1
    end

    if (pos < sz) and (st(pos)? == '.') then
      pos = pos + 1
      sep_seen = false
      while (pos < sz) and
            (((st(pos)? >= '0') and (st(pos)? <= '9')) or (st(pos)? == '_'))
      do
        if st(pos)? == '_' then
          if not sep_seen then
            sep_seen = true
            pos = pos + 1
            continue
          else
            Debug("[_MPFAlgo._from_string] consecutive '_' separators in fractional part")
          end
        end
        sep_seen = false
        let d: U8 = st(pos)? - '0'
        has_digit = true
        if (d != 0) or (sig_count > 0) then
          if sig_count < sig_limit then
            n_int = (n_int * ten_int) + MPInt.from[ILong](d.ilong())
            sig_count = sig_count + 1
            frac_count = frac_count + 1
          end
        else
          frac_count = frac_count + 1
        end
        pos = pos + 1
      end
    end

    if not has_digit then
      error
    end

    var str_exp: I64 = 0
    if (pos < sz) and
       ((st(pos)? == 'e') or (st(pos)? == 'E') or (st(pos)? == '@'))
    then
      pos = pos + 1
      var esign: Bool = false
      if (pos < sz) and (st(pos)? == '-') then
        esign = true ; pos = pos + 1
      elseif (pos < sz) and (st(pos)? == '+') then
        pos = pos + 1
      end
      var ehas_digit: Bool = false
      sep_seen = false
      while (pos < sz) and
            (((st(pos)? >= '0') and (st(pos)? <= '9')) or (st(pos)? == '_'))
      do
        if st(pos)? == '_' then
          if not sep_seen then
            sep_seen = true
            pos = pos + 1
            continue
          else
            Debug("[_MPFAlgo._from_string] consecutive '_' in exponent")
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
      if esign then str_exp = -str_exp end
    end

    if pos != sz then
      error
    end

    let dec_exp: I64 = (int_exp - frac_count) + str_exp
    if n_int.is_zero() then
      return MPFRep._create(fsign, false, false, 0, Array[U8].init(0, p_digits))
    end

    let p2: USize = p_digits + 2
    let n_max_exact: USize = ((p2.f64() * 14.0).usize() + 50).min(300)
    let ten_mp: MPFRep = MPFRep.from_f64(10.0, p2)

    let scaled: MPFRep =
      if dec_exp >= 0 then
        var n_mp: MPFRep = MPFRep.from_mpint(n_int, p2)
        if dec_exp > 0 then
          var scale: MPFRep = _one(p2)
          var sbase: MPFRep = ten_mp
          var sn: I64 = dec_exp
          while sn > 0 do
            if (sn and 1) == 1 then
              scale = _mul(scale, sbase, p2)
            end
            sbase = _mul(sbase, sbase, p2)
            sn = sn / 2
          end
          n_mp = _mul(n_mp, scale, p2)
        end
        n_mp._trunc(p_digits)
      else
        let n: USize = (-dec_exp).usize()
        if n <= n_max_exact then
          let extra: USize = (((n * 5) + 11) / 12) + p_digits + 4
          let shift_bits: USize = (extra * 8) - n
          let five_pow: MPInt = MPInt.from[ILong](5).pow(MPInt.from[ILong](n.ilong()))
          let n_shifted: MPInt = n_int.shl(MPInt.from[ILong](shift_bits.ilong()))
          (let q, let r) = n_shifted.divrem(five_pow)
          let q_rounded: MPInt = if (r + r) >= five_pow then
              q + MPInt.from[ILong](1)
            else
              q
            end
          let qbytes: Array[U8] val = q_rounded.raw_digits()
          let qlen: USize = qbytes.size()
          let new_exp: I64 = qlen.i64() - extra.i64()
          let keep: USize = p_digits.min(qlen)
          let new_digits: Array[U8] val = recover
            let d = Array[U8].create(keep)
            for i in Range(0, keep) do
              try d.push(qbytes(i)?) end
            end
            d
          end
          MPFRep._create(false, false, false, new_exp, new_digits)
        else
          var n_mp: MPFRep = MPFRep.from_mpint(n_int, p2)
          var scale: MPFRep = _one(p2)
          var sbase: MPFRep = ten_mp
          var sn: I64 = (-dec_exp)
          while sn > 0 do
            if (sn and 1) == 1 then
              scale = _mul(scale, sbase, p2)
            end
            sbase = _mul(sbase, sbase, p2)
            sn = sn / 2
          end
          _div(n_mp, scale, p2)._trunc(p_digits)
        end
      end

    MPFRep._create(fsign, scaled.is_nan(), scaled.is_infinite(),
      scaled.exponent(), scaled.raw_digits())


  // ── Ziv's rounding-determination test ────────────────────────────────────

  fun _ziv_ok(result: MPFRep, p: USize, mode: RoundingMode): Bool =>
    """
    Return `true` when the rounding direction of `result` at output precision
    `p` bytes can be determined unambiguously from the bits already computed.

    This is the core predicate of Ziv's iteration: if the first discarded byte
    is well away from every rounding boundary for the requested mode, we know
    that computing more guard bits cannot change the rounded result.

    Boundary locations by mode (in the first discarded byte):
    - `RoundingNearest` / `RoundingFaithful`: boundary at `0x80` (exact
      halfway). Ambiguous when `first_discarded == 0x80` AND all remaining
      discarded bytes are zero (exact tie) — extremely rare in practice.
    - `RoundingZero`: boundary at `0x00`. Ambiguous only when ALL discarded
      bytes are zero (result is exactly representable).
    - `RoundingPosInf` / `RoundingNegInf` / `RoundingAwayZ`: same boundary
      at `0x00`.

    We conservatively declare ambiguity when the first discarded byte is in
    `[0x7C, 0x84]` — within 4 of the nearest boundary — giving 3 extra bits
    of safety margin beyond the hard boundary. This avoids spurious retries
    while still catching all cases where the rounded result could change.
    """
    let digits = result.raw_digits()
    if digits.size() <= p then
      return true  // result already fits exactly; no rounding needed
    end
    let first_disc: U8 = try digits(p)? else 0 end
    var sticky: Bool = false
    for i in Range(p + 1, digits.size()) do
      if (try digits(i)? else 0 end) != 0 then
        sticky = true
        break
      end
    end
    match mode
    | RoundingNearest | RoundingFaithful =>
      // Ambiguous only at exact tie: first_disc == 0x80 and no sticky bits.
      not ((first_disc == 0x80) and not sticky)
    else
      // For directed modes, ambiguous only when ALL discarded bytes are zero
      // (result is exactly representable — no increment needed but cannot
      // distinguish from a result infinitesimally above a boundary).
      not ((first_disc == 0) and not sticky)
    end


  // ── Guard-byte table ──────────────────────────────────────────────────────

  fun _guard_bytes(ctx: MPFContext, op: String): USize =>
    """
    Working precision in bytes for the named operation `op` under context `ctx`.

    Returns `ctx.p_bytes() + guard`, where `guard` is a per-operation constant
    chosen to keep the rounding error of the full computation ≤ 1 ULP at
    output precision.

    The guard values are derived from the analysis in the MPFR thesis
    (C. Lauter, https://www.christoph-lauter.org/these.pdf) for exact rounding
    at F64 precision, and generalised conservatively to arbitrary precision:

    | `op`                        | guard | Reason                                  |
    |-----------------------------|-------|-----------------------------------------|
    | `"add"`, `"sub"`, `"mul"`   |   2   | Cancellation / FFT rounding ≤ 1 ULP     |
    | `"inv"`, `"sqrt"`           |   2   | Newton error ≤ 1 ULP at convergence     |
    | `"ln"`                      |   6   | Argument reduction cancellation         |
    | `"exp"`, `"trig"`, `"pi"`   |   8   | Reduction + Taylor / Machin series      |
    | (default)                   |   4   | Conservative fallback                   |

    This method is the single authoritative source for guard bytes. Adding a
    new operation requires only extending this match — no other file changes.
    `MPFContext.working_bytes` delegates here.
    """
    let guard: USize = match op
      | "add" | "sub" | "mul" | "inv" | "sqrt" => 2
      | "ln" => 6
      | "exp" | "trig" | "pi" => 8
      else 4
    end
    ctx.p_bytes() + guard


  // ── Public context-aware API ───────────────────────────────────────────────
  //
  // Each method below takes `ctx: MPFContext` as its first argument, handles
  // special values (NaN, ±∞, ±0), calls the corresponding internal `_*` method
  // at working precision `ctx.working_bytes(op)`, then rounds the result to
  // output precision via `ctx._round_to`.
  //
  // `MPFloat` calls these directly — there is no dispatch layer in `MPFContext`.

  fun add(ctx: MPFContext, a: MPFRep, b: MPFRep): MPFRep =>
    """
    Compute `a + b` at output precision `ctx.p_bytes()` using rounding
    `ctx.rounding`.

    Special value rules:
    - NaN propagates: `NaN + x = NaN`.
    - `+∞ + (−∞) = NaN`.
    - `±∞ + finite = ±∞`.
    - Zero is the additive identity.

    Finite non-zero operands are added via `_add` at working precision, then
    rounded to output precision via `ctx._round_to`.

    **Rounding**: correctly rounded — IEEE 754 §5.4. The intermediate sum at
    working precision is exact (no cancellation beyond digit alignment);
    a single `_round_to` step produces the nearest representable value.
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
      return ctx._round_to(b, ctx.p_bytes(), ctx.rounding)
    end
    if b.is_zero() then
      return ctx._round_to(a, ctx.p_bytes(), ctx.rounding)
    end
    let w = ctx.working_bytes("add")
    ctx._round_to(_add(a, b, w), ctx.p_bytes(), ctx.rounding)


  fun sub(ctx: MPFContext, a: MPFRep, b: MPFRep): MPFRep =>
    """
    Compute `a − b` at output precision `ctx.p_bytes()` using rounding
    `ctx.rounding`.

    Delegates to `add` after negating `b`, inheriting all special-value
    handling from `add`.
    """
    add(ctx, a, _neg_rep(b))


  fun mul(ctx: MPFContext, a: MPFRep, b: MPFRep): MPFRep =>
    """
    Compute `a × b` at output precision `ctx.p_bytes()` using rounding
    `ctx.rounding`.

    Special value rules: NaN propagates; `∞ × 0 = NaN`; `∞ × finite = ±∞`;
    `0 × finite = ±0`. Finite non-zero operands are multiplied via FFT
    convolution at working precision, then rounded to output precision.

    **Rounding**: correctly rounded — IEEE 754 §5.4. FFT convolution produces
    the exact full-length product; a single `_round_to` step rounds to output
    precision.
    """
    if a.is_nan() or b.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() or b.is_infinite() then
      if a.is_zero() or b.is_zero() then
        return MPFRep.nan_val()
      end
      return MPFRep._create(a.sign_bit() != b.sign_bit(), false, true, 0,
        Array[U8].create())
    end
    if a.is_zero() or b.is_zero() then
      return MPFRep._create(a.sign_bit() != b.sign_bit(), false, false, 0,
        Array[U8].init(0, ctx.p_bytes()))
    end
    let w = ctx.working_bytes("mul")
    ctx._round_to(_mul(a, b, w), ctx.p_bytes(), ctx.rounding)


  fun inv(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `1 / a` at output precision `ctx.p_bytes()`.

    NaN → NaN. ±∞ → ±0. 0 → ±∞.

    **Rounding**: correctly rounded — IEEE 754 §5.4. `_inv` runs Newton to
    convergence at `p + 2` bytes, then one extra refinement step at `p + 4`
    bytes halves the error below 0.5 ULP, making `_round_to` unambiguous.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(a.sign_bit(), false, false, 0,
        Array[U8].init(0, ctx.p_bytes()))
    end
    if a.is_zero() then
      return MPFRep._create(a.sign_bit(), false, true, 0, Array[U8].create())
    end
    let w = ctx.working_bytes("inv")
    ctx._round_to(_inv(a, w), ctx.p_bytes(), ctx.rounding)


  fun div(ctx: MPFContext, a: MPFRep, b: MPFRep): MPFRep =>
    """
    Compute `a / b` at output precision `ctx.p_bytes()`.

    Handles all IEEE 754 special cases: NaN propagation, ±∞/±∞ = NaN,
    finite/±∞ = ±0, ±∞/finite = ±∞, finite/0 = ±∞, 0/0 = NaN.

    **Rounding**: correctly rounded — IEEE 754 §5.4. Computed as `a × inv(b)`;
    `_inv`'s refinement step ensures the combined error is below 0.5 ULP,
    making the subsequent `_round_to` unambiguous.
    """
    if a.is_nan() or b.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() and b.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(a.sign_bit() != b.sign_bit(), false, true, 0,
        Array[U8].create())
    end
    if b.is_infinite() then
      return MPFRep._create(a.sign_bit() != b.sign_bit(), false, false, 0,
        Array[U8].init(0, ctx.p_bytes()))
    end
    if b.is_zero() then
      if a.is_zero() then
        return MPFRep.nan_val()
      end
      return MPFRep._create(a.sign_bit() != b.sign_bit(), false, true, 0,
        Array[U8].create())
    end
    if a.is_zero() then
      return MPFRep._create(a.sign_bit() != b.sign_bit(), false, false, 0,
        Array[U8].init(0, ctx.p_bytes()))
    end
    let w = ctx.working_bytes("inv")
    ctx._round_to(_div(a, b, w), ctx.p_bytes(), ctx.rounding)


  fun sqrt(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `√a` at output precision `ctx.p_bytes()`.

    NaN → NaN. Negative finite → NaN. +∞ → +∞. ±0 → ±0.

    **Rounding**: correctly rounded — IEEE 754 §5.4. `_sqrt` runs the
    reciprocal-sqrt Newton loop to convergence at `p + 2` bytes, then one
    extra refinement step at `p + 4` bytes halves the error below 0.5 ULP,
    making `_round_to` unambiguous.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      if a.sign_bit() then
        return MPFRep.nan_val()
      end
      return MPFRep.inf_val(true)
    end
    if a.is_zero() then
      return MPFRep._create(a.sign_bit(), false, false, 0,
        Array[U8].init(0, ctx.p_bytes()))
    end
    if a.sign_bit() then
      return MPFRep.nan_val()
    end
    let w = ctx.working_bytes("sqrt")
    ctx._round_to(_sqrt(a, w), ctx.p_bytes(), ctx.rounding)


  fun divrem(ctx: MPFContext, a: MPFRep, b: MPFRep): (MPFRep, MPFRep) =>
    """
    Truncated division with remainder: `(q, r)` such that
    `a = q × b + r`, `q = trunc(a / b)`, `r` has the same sign as `a`.

    All results are calculated at working precision and rounded at `ctx.p_bytes()`
    using `ctx.rounding`. Special cases follow IEEE 754 conventions.
    """
    let p = ctx.p_bytes()
    let w = ctx.working_bytes("inv")
    let nan = MPFRep.nan_val()
    if a.is_nan() or b.is_nan() then
      return (nan, nan)
    end
    if a.is_infinite() and b.is_infinite() then
      return (nan, nan)
    end
    if a.is_infinite() then
      return (MPFRep._create(a.sign_bit() != b.sign_bit(), false, true, 0, Array[U8].create()),
              nan)
    end
    if b.is_infinite() then
      return (MPFRep._create(false, false, false, 0, Array[U8].init(0, p)),
        MPFRep._create(a.sign_bit(), false, false, a.exponent(), a.raw_digits()))
    end
    if b.is_zero() then
      if a.is_zero() then
        return (nan, nan)
      end
      return (MPFRep._create(a.sign_bit() != b.sign_bit(), false, true, 0, Array[U8].create()),
              nan)
    end
    if a.is_zero() then
      let z = MPFRep._create(false, false, false, 0, Array[U8].init(0, p))
      return (z, z)
    end
    (let q, let r) = _divrem(a, b, w)
    (ctx._round_to(q, p, ctx.rounding), ctx._round_to(r, p, ctx.rounding))


  fun fld(ctx: MPFContext, a: MPFRep, b: MPFRep): MPFRep =>
    """
    Floored division: `floor(a / b)` at output precision `ctx.p_bytes()`.

    The result is the largest integer `q` such that `q × b ≤ a`. When the
    truncated remainder is non-zero and `a` and `b` have opposite signs, the
    truncated quotient is decremented by 1.

    Special cases are inherited from `divrem`.
    """
    let w = ctx.working_bytes("inv")
    let p = ctx.p_bytes()
    (let q, let r) = divrem(ctx, a, b)
    if r.is_nan() then
      return MPFRep.nan_val()
    end
    if r._has_frac() or (not r.is_zero()) then
      if a.sign_bit() != b.sign_bit() then
        let one = _one(p)
        return ctx._round_to(_sub(q, one, w), p, ctx.rounding)
      end
    end
    q


  fun rem(ctx: MPFContext, a: MPFRep, b: MPFRep): MPFRep =>
    """
    Truncated remainder: `a − trunc(a/b) × b`. Sign of result = sign of `a`.

    Special cases are inherited from `divrem`.
    """
    divrem(ctx, a, b)._2


  fun mod(ctx: MPFContext, a: MPFRep, b: MPFRep): MPFRep =>
    """
    Floored remainder: `a − fld(a,b) × b`. Sign of result = sign of `b`.

    Special cases are inherited from `rem`; when `rem` returns NaN the result
    is NaN.
    """
    let r = rem(ctx, a, b)
    if r.is_nan() then
      return MPFRep.nan_val()
    end
    if (not r.is_zero()) and (a.sign_bit() != b.sign_bit()) then
      let w = ctx.working_bytes("add")
      return ctx._round_to(_add(r, b, w), ctx.p_bytes(), ctx.rounding)
    end
    r


  fun ln(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `ln(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +∞. 0 → −∞. Negative finite → NaN.

    **Rounding**: correctly rounded — the result is the nearest representable
    value to the true `ln(a)` under `ctx.rounding`. Ziv's iteration
    automatically increases working precision until the rounding direction is
    unambiguous, so the rounded result is identical to what infinite-precision
    arithmetic followed by rounding would produce.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      if a.sign_bit() then
        return MPFRep.nan_val()
      end
      return MPFRep.inf_val(true)
    end
    if a.is_zero() then
      return MPFRep.inf_val(false)
    end
    if a.sign_bit() then
      return MPFRep.nan_val()
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("ln")
    var result = _ln(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _ln(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun log(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Natural logarithm alias — same as `ln`. Correctly rounded.
    """
    ln(ctx, a)


  fun log2(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `log₂(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +∞. 0 → −∞. Negative finite → NaN.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      if a.sign_bit() then
        return MPFRep.nan_val()
      end
      return MPFRep.inf_val(true)
    end
    if a.is_zero() then
      return MPFRep.inf_val(false)
    end
    if a.sign_bit() then
      return MPFRep.nan_val()
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("ln")
    var result = _log2(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _log2(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun log10(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `log₁₀(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +∞. 0 → −∞. Negative finite → NaN.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      if a.sign_bit() then
        return MPFRep.nan_val()
      end
      return MPFRep.inf_val(true)
    end
    if a.is_zero() then
      return MPFRep.inf_val(false)
    end
    if a.sign_bit() then
      return MPFRep.nan_val()
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("ln")
    var result = _log10(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _log10(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun logb(ctx: MPFContext, a: MPFRep, base: MPFRep): MPFRep =>
    """
    Compute `log_base(a)` at output precision `ctx.p_bytes()`.

    NaN in either operand → NaN.

    **Rounding**: faithfully rounded (≤ 1 ULP). Implemented as `ln(a)/ln(base)`;
    each `ln` call is correctly rounded but the final division introduces an
    additional rounding step, so the overall result may differ from the
    correctly-rounded value by at most 1 ULP.
    """
    if a.is_nan() or base.is_nan() then
      return MPFRep.nan_val()
    end
    let w = ctx.working_bytes("ln")
    ctx._round_to(_logb(a, base, w), ctx.p_bytes(), ctx.rounding)


  fun exp(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `eᵃ` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +∞. −∞ → +0. 0 → 1.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      if a.sign_bit() then
        return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
      end
      return MPFRep.inf_val(true)
    end
    if a.is_zero() then
      return _one(ctx.p_bytes())
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("exp")
    var result = _exp(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _exp(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun exp2(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `2ᵃ` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +∞. −∞ → +0. 0 → 1.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      if a.sign_bit() then
        return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
      end
      return MPFRep.inf_val(true)
    end
    if a.is_zero() then
      return _one(ctx.p_bytes())
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("exp")
    var result = _exp2(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _exp2(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun powi(ctx: MPFContext, a: MPFRep, n: ILong): MPFRep =>
    """
    Compute `aⁿ` for integer exponent `n` at output precision `ctx.p_bytes()`.

    - NaN → NaN.
    - ±∞, n > 0 → ±∞ (sign = sign of a when n is odd, positive otherwise).
    - ±∞, n < 0 → ±0.
    - ±∞, n = 0 → 1 (by convention).

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    Implemented via repeated squaring; each squaring is a multiply whose
    individual rounding is faithfully rounded, but Ziv's outer loop ensures
    the final result is correctly rounded.
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      if n == 0 then
        return _one(ctx.p_bytes())
      elseif n > 0 then
        let result_sign = a.sign_bit() and ((n and 1) != 0)
        return MPFRep._create(result_sign, false, true, 0, Array[U8].create())
      else
        return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
      end
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("exp")
    var result = _powi(a, n, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _powi(a, n, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun pow(ctx: MPFContext, a: MPFRep, b: MPFRep): MPFRep =>
    """
    Compute `aᵇ` at output precision `ctx.p_bytes()`.

    NaN in either → NaN. b = 0 → 1. a = 0 and b > 0 → 0. a = 0 and b ≤ 0 → NaN.
    Negative base with non-integer exponent → NaN.
    For positive base: `exp(b × ln(a))`. For negative integer base: delegates to `powi`.

    **Rounding**: correctly rounded via Ziv's iteration (see `ln`).
    """
    if a.is_nan() or b.is_nan() then
      return MPFRep.nan_val()
    end
    if b.is_zero() then
      return _one(ctx.p_bytes())
    end
    if a.is_zero() then
      if b.sign_bit() then
        return MPFRep.nan_val()
      end
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    if a.sign_bit() and b._has_frac() then
      return MPFRep.nan_val()
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("exp")
    var result = _pow(a, b, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _pow(a, b, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun sin(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `sin(a)` at output precision `ctx.p_bytes()`.

    NaN or ±∞ → NaN. 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_sin(a, w), ctx.p_bytes(), ctx.rounding)


  fun cos(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `cos(a)` at output precision `ctx.p_bytes()`.

    NaN or ±∞ → NaN. 0 → 1.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return _one(ctx.p_bytes())
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_cos(a, w), ctx.p_bytes(), ctx.rounding)


  fun tan(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `tan(a)` at output precision `ctx.p_bytes()`.

    NaN or ±∞ → NaN. 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_tan(a, w), ctx.p_bytes(), ctx.rounding)


  fun sinh(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `sinh(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. ±∞ → ±∞. 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(a.sign_bit(), false, true, 0, Array[U8].create())
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_sinh(a, w), ctx.p_bytes(), ctx.rounding)


  fun cosh(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `cosh(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. ±∞ → +∞. 0 → 1.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep.inf_val(true)
    end
    if a.is_zero() then
      return _one(ctx.p_bytes())
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_cosh(a, w), ctx.p_bytes(), ctx.rounding)


  fun tanh(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `tanh(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +1. −∞ → −1. 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(a.sign_bit(), false, false, 1,
        Array[U8].init(1, 1))
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_tanh(a, w), ctx.p_bytes(), ctx.rounding)


  fun sech(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `sech(a) = 1/cosh(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. ±∞ → +0.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_sech(a, w), ctx.p_bytes(), ctx.rounding)


  fun csch(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `csch(a) = 1/sinh(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. ±∞ → +0. 0 → NaN.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    if a.is_zero() then
      return MPFRep.nan_val()
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_csch(a, w), ctx.p_bytes(), ctx.rounding)


  fun coth(ctx: MPFContext, a: MPFRep): MPFRep =>
    """
    Compute `coth(a) = cosh(a)/sinh(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +1. −∞ → −1. 0 → NaN.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      return MPFRep._create(a.sign_bit(), false, false, 1,
        Array[U8].init(1, 1))
    end
    if a.is_zero() then
      return MPFRep.nan_val()
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_coth(a, w), ctx.p_bytes(), ctx.rounding)


  fun pi(ctx: MPFContext): MPFRep =>
    """
    Compute π at output precision `ctx.p_bytes()`.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    let w = ctx.working_bytes("pi")
    ctx._round_to(_pi(w), ctx.p_bytes(), ctx.rounding)


  fun pi_bbp(ctx: MPFContext): MPFRep =>
    """
    Compute π using the Bailey–Borwein–Plouffe formula at output precision
    `ctx.p_bytes()`.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    let w = ctx.working_bytes("pi")
    ctx._round_to(_pi_bbp(w), ctx.p_bytes(), ctx.rounding)


  fun pi_chudnovsky(ctx: MPFContext): MPFRep =>
    """
    Compute π using the Chudnovsky algorithm at output precision `ctx.p_bytes()`.

    ⚠ Warning: this implementation currently gives only ~4 correct digits —
    see `_pi_chudnovsky` for the known bug description.
    """
    let w = ctx.working_bytes("pi")
    ctx._round_to(_pi_chudnovsky(w), ctx.p_bytes(), ctx.rounding)


  fun from_string(ctx: MPFContext, s: String): MPFRep ? =>
    """
    Parse `s` as a floating-point number at output precision `ctx.p_bytes()`.

    Raises an error if `s` is not a valid floating-point representation.
    Only base 10 is currently implemented.
    """
    let w = ctx.working_bytes("add")
    ctx._round_to(_from_string(s, w)?, ctx.p_bytes(), ctx.rounding)
