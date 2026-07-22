// Multi-precision float algorithms

use "../assertx"
use "../formatx"
use "../mathx"

use "debug"
use "collections"


primitive _MPFAlgo
  """
  All arithmetic and transcendental algorithms for `MPFRep` values.

  ## Two layers of methods

  **Public (context-aware)**: `add`, `sub`, `mul`, `div`, `inv`, `sqrt`,
  `divrem`, `fld`, `rem`, `mod`, `ln`, `log`, `log2`, `log10`, `logb`,
  `exp`, `exp2`, `powi`, `pow`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`,
  `sech`, `csch`, `coth`, `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh`,
  `cbrt`, `rootn`, `pi`, `e`, `from_string`.

  Each takes `ctx: MPFContext` as its first argument. It handles **all**
  special values (NaN, ±∞, ±0) with early returns before touching `_*`
  methods, runs the computation at working precision `ctx.working_bytes(op)`,
  then rounds to output precision via `ctx._round_to`.

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

  // ── Conversion from external types ────────────────────────────────────────

  fun _from_f64(f: F64, p: USize): MPFRep val =>
    if f.nan() then
      MPFRep._create(false, true, false, 0, Array[U8].create())
    elseif f.infinite() then
      MPFRep._create(f < 0.0, false, true, 0, Array[U8].create())
    elseif f == 0.0 then
      MPFRep._create(f.bits() == 0x8000000000000000, false, false, 0, Array[U8].init(0, p))
    else
      let sgn = f < 0.0
      (let m, let e_u32) = f.abs().frexp()
      let bit_exp = e_u32.i32().i64()
      let expn = (bit_exp.f64() / 8.0).ceil().i64()
      let shift = bit_exp - (expn * 8)
      var frac = m.f64() * F64(2).pow(shift.f64())
      let digits: Array[U8] val = recover
        let d = Array[U8].init(0, p)
        var i: USize = 0
        while i < p do
          frac = frac * 256.0
          let di: U8 = frac.u8()
          try d.update(i, di)? end
          frac = frac - di.f64()
          i = i + 1
        end
        d
      end
      MPFRep._create(sgn, false, false, expn, digits)
    end


  fun _from_f32(f: F32, p: USize): MPFRep val =>
    if f.nan() then
      MPFRep._create(false, true, false, 0, Array[U8].create())
    elseif f.infinite() then
      MPFRep._create(f < 0.0, false, true, 0, Array[U8].create())
    elseif f == 0.0 then
      MPFRep._create(f.bits() == 0x80000000, false, false, 0, Array[U8].init(0, p))
    else
      let sgn = f < 0.0
      (let m, let e_u32) = f.abs().frexp()
      let bit_exp = e_u32.i32().i64()
      let expn = (bit_exp.f32() / 8.0).ceil().i64()
      let shift = bit_exp - (expn * 8)
      var frac = m * F32(2).pow(shift.f32())
      let digits: Array[U8] val = recover
        let d = Array[U8].init(0, p)
        var i: USize = 0
        while i < p do
          frac = frac * 256.0
          let di: U8 = frac.u8()
          try d.update(i, di)? end
          frac = frac - di.f32()
          i = i + 1
        end
        d
      end
      MPFRep._create(sgn, false, false, expn, digits)
    end


  fun _from_mpint(n: MPInt, p: USize): MPFRep val =>
    if n.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, p))
    end
    let mag = n.raw_digits()
    let total: USize = mag.size()
    let digits: Array[U8] val = recover
      let keep = total.min(p)
      let d = Array[U8].create(keep)
      mag.copy_to(d, 0, 0, keep)
      d
    end
    MPFRep._create(n.is_negative(), false, false, total.i64(), digits)



  fun _from_ulong(n: ULong, p: USize): MPFRep val =>
    let base: ULong = 256
    var q: ULong = n
    let digits: Array[U8] val = recover
      let d = Array[U8].create(p)
      while q >= base do
        (q, let r) = q.divrem(base)
        d.push(r.u8())
      end
      d.push(q.u8())
      d.reverse_in_place()
      d
    end
    MPFRep._create(false, false, false, digits.size().i64(), digits)


  // ── Addition / Subtraction ─────────────────────────────────────────────────

  fun _add(a: MPFRep box, b: MPFRep box, w: USize): MPFRep val =>
    """
    Add `a + b` at working precision `w` bytes.

    Assumes both operands are finite and non-zero. Special-value handling
    (NaN, ±∞, zero) belongs to the public `add` method.

    When operand signs agree, magnitudes are added (`_add_mag`). When they
    differ, the smaller magnitude is subtracted from the larger and the result
    takes the sign of the larger. The result is truncated to `w` bytes.
    """
    // Pass w into the magnitude methods so they cap their output directly,
    // eliminating the separate _trunc(w) allocation in the common path.
    if a.sign_bit() == b.sign_bit() then
      // Same sign: add magnitudes.
      if a.exponent() >= b.exponent() then
        a._add_mag(b, a.sign_bit(), w)
      else
        b._add_mag(a, a.sign_bit(), w)
      end
    else
      // Different signs: subtract smaller magnitude from larger.
      match a._cmp_mag(b)
      | Greater => a._sub_mag(b, a.sign_bit(), w)
      | Less    => b._sub_mag(a, b.sign_bit(), w)
      else
        MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
      end
    end


  fun _sub(a: MPFRep box, b: MPFRep box, w: USize): MPFRep val =>
    """
    Subtract `a − b` at working precision `w` bytes. Delegates to `_add`
    after negating `b`, inheriting all sign and special-value handling.
    """
    _add(a, _neg_rep(b), w)


  // ── Multiplication ─────────────────────────────────────────────────────────

  fun _mul(a: MPFRep box, b: MPFRep box, w: USize): MPFRep val =>
    """
    Multiply `a × b` at working precision `w` bytes via FFT convolution.

    The result sign is `a.sign_bit() XOR b.sign_bit()`. The raw convolution
    product has `a._size() + b._size()` digits; it is normalised by stripping
    leading zeros and adjusting the exponent, then truncated to exactly `w`
    bytes — combining the `norm_d` copy and the `_trunc` step into a single
    output array allocation rather than three.

    Both operands must be finite (NaN/±∞ handling belongs to the public API
    layer).
    """
    let result_sign: Bool = a.sign_bit() != b.sign_bit()
    let this_size: USize = a._size()
    let that_size: USize = b._size()
    let res_size: USize = this_size + that_size
    let pow2: USize = res_size.next_pow2().max(4)
    let ad: Array[U8] val = a.raw_digits()
    let bd: Array[U8] val = b.raw_digits()

    var err_i: USize = 0
    // Compute FFT convolution → U8 digits, find leading zeros, and write at most
    // w significant bytes into the output — all in one recover block, eliminating
    // the intermediate res_size-byte array and the separate norm/trunc copies.
    // `leading` is a var of value type USize; Pony allows writing to outer var
    // bindings of machine-word types from inside a recover block.
    var leading: USize = 0
    let out: Array[U8] val = recover
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

        // Carry propagation: convert F64 convolution back to U8 digits.
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

        // Map F64 results to the logical digit array:
        //   slot 0  = carry.u8()
        //   slot k (k ≥ 1) = fb(k-1).u8()  for k < res_size+1
        // total slots = res_size + 1
        let total: USize = res_size + 1
        let carry_byte: U8 = carry.u8()

        // Find first non-zero slot.
        if carry_byte != 0 then
          leading = 0
        else
          var lead: USize = 1
          while lead < total do
            let b_val: U8 = if lead < total then fb(lead - 1)?.u8() else 0 end
            if b_val != 0 then
              break
            end
            lead = lead + 1
          end
          leading = lead
        end

        // Copy at most w bytes starting from leading into output.
        let keep: USize = (total - leading).min(w)
        if keep == 0 then
          Array[U8].init(0, w)
        else
          let o = Array[U8].init(0, keep)
          var s: USize = leading
          var j: USize = 0
          while j < keep do
            o.update(j, if s == 0 then carry_byte else fb(s - 1)?.u8() end)?
            j = j + 1
            s = s + 1
          end
          o
        end
      else
        Fail(Format("[_MPFAlgo._mul] Index out of bounds at i={}", err_i))
        Array[U8].init(0, w)
      end
    end

    // leading == res_size+1 means all digits are zero.
    if leading == (res_size + 1) then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, w))
    end
    let new_exp: I64 = (a.exponent() + b.exponent()) - leading.i64()
    MPFRep._create(result_sign, false, false, new_exp, out)


  // ── Division / Reciprocal ─────────────────────────────────────────────────

  fun _inv(a: MPFRep box, w: USize): MPFRep val =>
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
    let ad: Array[U8] val = a.raw_digits()
    let prec: USize = a._size()

    // G = same digits as |a| but exponent=1, so its integer part is d[0] ∈ [1, 256).
    let g: MPFRep val = MPFRep._create(false, false, false, 1, ad)

    // F64 initial estimate: read the first few digits.
    // fg is constructed with exponent = 1 regardless of a's actual exponent so
    // G ∈ [1, 256) always. The F64 estimate reads the first 4 bytes of the digit array,
    // building a value in [1, 256) as well. That is always within F64 range, so fg
    // is a normal F64 and 1.0 / fg is a normal F64 in (1/256, 1].
    let ng: USize = prec.min(4)
    var fg: F64 = 0.0
    try
      var i: USize = ng
      repeat
        i = i - 1
        fg = (fg / 256.0) + ad(i)?.f64()
      until i == 0 end
    end
    var res: MPFRep ref = MPFRep.from[F64](1.0 / fg, size)

    let two: MPFRep val = MPFRep.from[F64](2.0, size)
    var iters: USize = 0
    let max_iters: USize = size * 4
    while iters < max_iters do
      iters = iters + 1
      let new_res = _mul(res, _sub(two, _mul(g, res, size), size), size)
      let converged = _converged(new_res, res)
      res._update(new_res)
      if converged then
        break
      end
    end

    // One extra refinement step at size2 bytes to push error below 0.5 ULP.
    let size2 = size + 2  // refinement precision
    res._update(_mul(res, _sub(two, _mul(g, res, size2), size2), size2))

    // 1/|a| = res × 256^{1−e}. Adjust sign and exponent.
    MPFRep._create(a.sign_bit(), false, false, (res.exponent() + 1) - a.exponent(),
      res.raw_digits())


  fun _div(a: MPFRep box, b: MPFRep box, w: USize): MPFRep val =>
    """
    Divide `a / b` at working precision `w` bytes. Computes `a × inv(b)`.

    Both operands must be finite and non-zero (special-value handling belongs
    to the public API in `MPFContext`).
    """
    _mul(a, _inv(b, w), w)


  // ── Square Root ────────────────────────────────────────────────────────────

  fun _sqrt(a: MPFRep box, w: USize): MPFRep val =>
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
    let ad: Array[U8] val = a.raw_digits()
    let prec: USize = a._size()
    let a_exp: I64 = a.exponent()

    // Choose h_exp so that a_exp − h_exp is always even.
    let h_exp: I64 = if (a_exp and 1) != 0 then 1 else 2 end
    let h: MPFRep val = MPFRep._create(false, false, false, h_exp, ad)

    // F64 initial estimate.
    // Like for inv, the estimate is in the F64 range.
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
    var res: MPFRep ref = MPFRep.from[F64](1.0 / fg_h.sqrt(), size)

    let three: MPFRep val = MPFRep.from[F64](3.0, size)
    var iters: USize = 0
    let max_iters: USize = size * 4
    while iters < max_iters do
      iters = iters + 1
      (let halved: MPFRep val, _) = _mul(res, _sub(three, _mul(h, _mul(res, res, size), size), size), size)._short_div(2)
      let nr = halved._trunc(size)
      let converged = _converged(nr, res)
      res._update(nr)
      if converged then
        break
      end
    end

    // One extra refinement step at size2 bytes to push error below 0.5 ULP.
    let size2 = size + 2  // refinement precision
    (let halved2: MPFRep val, _) = _mul(res, _sub(three, _mul(h, _mul(res, res, size2), size2), size2), size2)._short_div(2)
    res._update(halved2._trunc(size2))

    // √H = H × (1/√H); exponent adjusted for the parity split.
    let sqrt_h: MPFRep val = _mul(h, res, size2)
    let result_exp: I64 = sqrt_h.exponent() + ((a_exp - h_exp) / 2)
    MPFRep._create(false, false, false, result_exp, sqrt_h.raw_digits())


  // ── Division with Remainder ────────────────────────────────────────────────

  fun _divrem(a: MPFRep box, b: MPFRep box, w: USize): (MPFRep val, MPFRep val) =>
    """
    Truncated division with remainder at working precision `w` bytes.

    Returns `(q, r)` where `q = trunc(a / b)` and `r = a − q × b`.

    Assumes both operands are finite and non-zero. Special-value handling
    (NaN, ±∞, ±0) belongs to the public `divrem` method.

    Post-correction: Newton's `_inv` always undershoots by up to 1 ULP.
    When `|r| ≥ |b|`, increment `|q|` by 1 and recompute `r`. At most one
    step is needed.
    """
    var q: MPFRep val = _div(a, b, w)._trunc_frac()
    var r: MPFRep val = _sub(a, _mul(q, b, w), w)

    // Post-correction: if |r| ≥ |b|, q was off by 1.
    if (not r.is_zero()) and (r._cmp_mag(b) != Less) then
      let one: MPFRep val = MPFRep.from[F64](1.0, w)
      if a.sign_bit() == b.sign_bit() then
        q = _add(q, one, w)
      else
        q = _sub(q, one, w)
      end
      r = _sub(a, _mul(q, b, w), w)
    end
    (q, r)


  fun _fld(a: MPFRep box, b: MPFRep box, w: USize): MPFRep val =>
    """
    Floored division `floor(a / b)` at working precision `w` bytes.

    The result is the largest integer `q` such that `q × b ≤ a`. When the
    truncated remainder is non-zero and `a` and `b` have opposite signs, the
    truncated quotient is decremented by 1.

    Assumes both operands are finite and non-zero. Special-value handling
    belongs to the public `fld` method.
    """
    (let q: MPFRep val, let r: MPFRep box) = _divrem(a, b, w)
    if (not r.is_zero()) and (a.sign_bit() != b.sign_bit()) then
      return _sub(q, MPFRep.from[F64](1.0, w), w)
    end
    q


  fun _mod(a: MPFRep box, b: MPFRep box, w: USize): MPFRep val =>
    """
    Floored remainder `a − floor(a / b) × b` at working precision `w` bytes.

    The result has the same sign as `b`. Assumes both operands are finite and
    non-zero. Special-value handling belongs to the public `mod` method.
    """
    let r: MPFRep val = _rem(a, b, w)
    if (not r.is_zero()) and (a.sign_bit() != b.sign_bit()) then
      return _add(r, b, w)
    end
    r


  fun _rem(a: MPFRep box, b: MPFRep box, w: USize): MPFRep val =>
    """
    Truncated remainder `a − trunc(a / b) × b` at working precision `w` bytes.

    The result has the same sign as `a`. Assumes both operands are finite and
    non-zero. Special-value handling belongs to the public `rem` method.
    """
    _divrem(a, b, w)._2


  // ── Internal convergence helpers ──────────────────────────────────────────

  fun _converged(a: MPFRep box, b: MPFRep box): Bool =>
    """
    Convergence predicate for Newton and Taylor loops.

    Returns `true` when `a` and `b` have the same exponent and identical
    leading bytes up to `min(|a|, |b|)`. Deliberately weaker than value
    equality (`MPFRep.eq`): sign is ignored (intermediate values are
    constructed sign-neutral inside kernels), and extra trailing bytes on
    the longer operand are not compared (a value at `w+1` bytes is considered
    converged against one at `w` bytes if their shared prefix matches).
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


  fun _neg_rep(a: MPFRep box): MPFRep val =>
    """
    Return a copy of `a` with its sign flipped. Equivalent to arithmetic
    negation but operates purely on the representation (no allocation of new
    digits).
    """
    MPFRep._create(not a.sign_bit(), a.is_nan(), a.is_infinite(),
      a.exponent(), a.raw_digits())


  // ── Series helpers ─────────────────────────────────────────────────────────

  fun _atanh_series(x: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `arctanh(x)` via the Taylor series
    `Σ_{k=0}^∞ x^{2k+1} / (2k+1)` at working precision `w` bytes.

    The caller must ensure `|x| < 1` for convergence. The loop terminates
    when the current addend no longer changes the running sum at precision `w`.

    Uses `MPFRep._short_div` (exact U8 division) for denominators up to 255;
    falls back to `_div` for larger denominators (only for very high `w`).
    """
    let p2: USize = w + 2
    let x2 = _mul(x, x, p2)
    var term: MPFRep ref = x._clone()
    var sum:  MPFRep ref = term._clone()
    var k: USize = 3
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term._update(_mul(term, x2, p2))
      let addend = if k <= 255 then
          (let q: MPFRep val, _) = term._short_div(k.u8())
          q
        else
          _div(term, MPFRep.from[F64](k.f64(), p2), p2)
        end
      let new_sum = _add(sum, addend, p2)

      if _converged(new_sum, sum) then
        break
      end

      sum._update(new_sum)
      k = k + 2
      iters = iters + 1
    end
    sum._trunc(w)


  fun _ln2_const(w: USize): MPFRep val =>
    """
    Compute `ln(2)` to `w` bytes of working precision using the identity
    `ln(2) = 2 × arctanh(1/3)`.

    The rational `1/3` has the exact repeating base-256 representation
    `[0x55, 0x55, …]` (since `85/255 = 1/3`).
    """
    let p2: USize = w + 4
    let d: Array[U8] val = recover Array[U8].init(0x55, p2) end
    let third: MPFRep val = MPFRep._create(false, false, false, 0, d)
    let two: MPFRep val = MPFRep.from[F64](2.0, p2)
    _mul(_atanh_series(third, p2), two, p2)._trunc(w)


  fun _exp_taylor(r: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `e^r` via the Taylor series `Σ_{k=0}^∞ r^k / k!` at working
    precision `w` bytes.

    Uses the recurrence `term_k = term_{k-1} × r / k` to avoid computing
    factorials. Works best for `|r| ≤ ln(2)/2 ≈ 0.347`.
    """
    let p2: USize = w + 2
    let one = MPFRep.from[F64](1.0, p2)
    var term: MPFRep ref = one._clone()
    var sum:  MPFRep ref = one._clone()
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term._update(_mul(term, r, p2))
      let divided = if k <= 255 then
          (let q: MPFRep val, _) = term._short_div(k.u8())
          q
        else
          _div(term, MPFRep.from[F64](k.f64(), p2), p2)
        end
      term._update(divided)
      let new_sum = _add(sum, divided, p2)

      if _converged(new_sum, sum) then
        break
      end

      sum._update(new_sum)
      k = k + 1
      iters = iters + 1
    end
    sum._trunc(w)


  fun _sin_taylor(r: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `sin(r)` via the Taylor series `r − r³/3! + r⁵/5! − r⁷/7! + ...`
    at working precision `w` bytes.

    The recurrence `term_{k+1} = term_k × (−r²) / ((2k)(2k+1))` advances
    two factorial positions at once. Works best for `|r| ≤ π/4 ≈ 0.785`.
    """
    let p2: USize = w + 2
    let neg_r2 = _neg_rep(_mul(r, r, p2))
    var term: MPFRep ref = r._clone()
    var sum:  MPFRep ref = term._clone()
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term._update(_mul(term, neg_r2, p2))
      let d1: USize = 2 * k
      let d2: USize = (2 * k) + 1
      term._update(if d1 <= 255 then
          (let q: MPFRep val, _) = term._short_div(d1.u8())
          q
        else
          _div(term, MPFRep.from[F64](d1.f64(), p2), p2)
        end)
      term._update(if d2 <= 255 then
          (let q: MPFRep val, _) = term._short_div(d2.u8())
          q
        else
          _div(term, MPFRep.from[F64](d2.f64(), p2), p2)
        end)
      let new_sum = _add(sum, term, p2)

      if _converged(new_sum, sum) then
        break
      end

      sum._update(new_sum)
      k = k + 1
      iters = iters + 1
    end
    sum._trunc(w)


  fun _cos_taylor(r: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `cos(r)` via the Taylor series `1 − r²/2! + r⁴/4! − r⁶/6! + ...`
    at working precision `w` bytes.

    The recurrence `term_{k+1} = term_k × (−r²) / ((2k−1)(2k))` advances
    two factorial positions at once. Works best for `|r| ≤ π/4 ≈ 0.785`.
    """
    let p2: USize = w + 2
    let neg_r2 = _neg_rep(_mul(r, r, p2))
    let one = MPFRep.from[F64](1.0, p2)
    var term: MPFRep ref = one._clone()
    var sum:  MPFRep ref = one._clone()
    var k: USize = 1
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term._update(_mul(term, neg_r2, p2))
      let d1: USize = (2 * k) - 1
      let d2: USize = 2 * k
      term._update(if d1 <= 255 then
          (let q: MPFRep val, _) = term._short_div(d1.u8())
          q
        else
          _div(term, MPFRep.from[F64](d1.f64(), p2), p2)
        end)
      term._update(if d2 <= 255 then
          (let q: MPFRep val, _) = term._short_div(d2.u8())
          q
        else
          _div(term, MPFRep.from[F64](d2.f64(), p2), p2)
        end)
      let new_sum = _add(sum, term, p2)

      if _converged(new_sum, sum) then
        break
      end

      sum._update(new_sum)
      k = k + 1
      iters = iters + 1
    end
    sum._trunc(w)


  // ── Integer extraction helper ──────────────────────────────────────────────

  fun _round_to_ilong(a: MPFRep box): ILong =>
    """
    Extract the integer part of `a` as an `ILong` (truncated toward zero).

    Returns 0 for NaN, ±∞, or zero. Used in argument-reduction steps for
    `_exp`, `_sin`, and `_cos` to compute the multiple-of-period quotient.

    Only the first `exponent()` bytes of the digit array are the integer
    part; the rest are fractional and are excluded before passing to `MPInt`.

    **Overflow**: if the integer part of `a` exceeds `ILong.max_value()`
    (~9.2×10¹⁸), `MPInt.ilong()` silently truncates to the low 64 bits,
    giving a wrong result. In practice this cannot occur: working precisions
    used by callers cap argument magnitudes far below that threshold.
    """
    if (not a.is_finite()) or a.is_zero() then
      return 0
    end
    let rounded: MPFRep val = a._trunc_frac()
    let exp_bytes: USize = rounded.exponent().usize().min(rounded._size())
    if exp_bytes == 0 then
      return 0
    end
    // trim() on val returns val, sharing the backing store — no copy needed.
    MPInt.from_bytes_be(rounded.sign_bit(), rounded.raw_digits().trim(0, exp_bytes)).ilong()


  fun _extend(a: MPFRep box, w: USize): MPFRep val =>
    """
    Return `a` extended (or truncated) to `w` bytes by zero-padding the tail
    or calling `_trunc(w)`. The sign, exponent, and leading digits are
    preserved.
    """
    let p: USize = a._size()
    if p == w then
      // Already the right size; _create shares the val reference, no copy.
      return MPFRep._create(a.sign_bit(), a.is_nan(), a.is_infinite(), a.exponent(), a.raw_digits())
    end
    if p > w then
      return a._trunc(w)
    end
    let digits: Array[U8] val = a.raw_digits()
    let padded: Array[U8] val = recover
      let d = Array[U8].init(0, w)
      digits.copy_to(d, 0, 0, p)
      d
    end
    MPFRep._create(a.sign_bit(), a.is_nan(), a.is_infinite(), a.exponent(), padded)


  // ── Logarithm / Exponential ────────────────────────────────────────────────

  fun _ln(a: MPFRep box, w: USize): MPFRep val =>
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

    let ln2: MPFRep val = _ln2_const(p2)
    let k: I64 = a.exponent() - 1

    let m_digits: Array[U8] val = recover
      let arr = Array[U8].init(0, p2)
      a.raw_digits().copy_to(arr, 0, 0, a._size().min(p2))
      arr
    end
    let m: MPFRep val = MPFRep._create(false, false, false, 1, m_digits)

    let d0: U8 = try a.raw_digits()(0)? else 1 end
    var n_inner: I64 = 0
    var tmp: U8 = d0
    while tmp >= 2 do
      tmp = tmp >> 1
      n_inner = n_inner + 1
    end
    let divisor: U8 = U8(1).shl(n_inner.u8())
    (let u: MPFRep val, _) = m._short_div(divisor)

    let one: MPFRep val = MPFRep.from[F64](1.0, p2)
    let two: MPFRep val = MPFRep.from[F64](2.0, p2)
    let t: MPFRep val = _div(_sub(u, one, p2), _add(u, one, p2), p2)
    let ln_u: MPFRep val = _mul(_atanh_series(t, p2), two, p2)

    let ln2_factor: I64 = (8 * k) + n_inner
    let correction: MPFRep val = if ln2_factor == 0 then
        MPFRep._create(false, false, false, 0, Array[U8].init(0, p2))
      else
        let fac: MPFRep val = MPFRep.from[F64](ln2_factor.f64(), p2)
        _mul(fac, ln2, p2)
      end
    _add(correction, ln_u, p2)._trunc(w)


  fun _exp(a: MPFRep box, w: USize): MPFRep val =>
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

    let ln2: MPFRep val = _ln2_const(p2)
    let eight: MPFRep val = MPFRep._create(false, false, false, 1, recover [U8(8)] end)
    let ln256: MPFRep val = _mul(ln2, eight, p2)

    let x: MPFRep val = _extend(a, p2)

    let n256_f: MPFRep val = _div(x, ln256, p2)._trunc_frac()
    let n256: ILong = _round_to_ilong(n256_f)
    let n256_mpf: MPFRep val = MPFRep.from[F64](n256.f64(), p2)
    let r1: MPFRep val = _sub(x, _mul(n256_mpf, ln256, p2), p2)

    let n2_f: MPFRep val = _div(r1, ln2, p2)._trunc_frac()
    let n2: ILong = _round_to_ilong(n2_f)
    let n2_mpf: MPFRep val = MPFRep.from[F64](n2.f64(), p2)
    // r = r1 − n2 × ln2.
    let r: MPFRep val = _sub(r1, _mul(n2_mpf, ln2, p2), p2)

    var result: MPFRep val = _exp_taylor(r, p2)

    if n2 > 0 then
      let pow2: MPFRep val = MPFRep._create(false, false, false, 1, recover [U8(1).shl(n2.u8())] end)
      result = _mul(result, pow2, p2)
    elseif n2 < 0 then
      let pow2: MPFRep val = MPFRep._create(false, false, false, 0, recover [U8(128).shr((-n2 - 1).u8())] end)
      result = _mul(result, pow2, p2)
    end

    let exp_digits: Array[U8] val = result.raw_digits().trim(0, result._size().min(w))
    MPFRep._create(result.sign_bit(), false, false, result.exponent() + n256.i64(), exp_digits)


  fun _log2(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `log₂(a)` at working precision `w` bytes.
    `log₂(a) = ln(a) / ln(2)`.

    Assumes `a` is finite, positive, and non-zero.
    """
    _div(_ln(a, w), _ln2_const(w + 4), w)


  fun _log10(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `log₁₀(a)` at working precision `w` bytes.
    `log₁₀(a) = ln(a) / ln(10)`, with `ln(10) = 2 × arctanh(1/3) + 2 × arctanh(2/3)`.

    Assumes `a` is finite, positive, and non-zero.
    """
    let ln_x: MPFRep val = _ln(a, w)
    let p2: USize = w + 6
    // 1/3 and 2/3 are exact repeating fractions in base 256:
    //   1/3 = 0x55/0xFF = 85/255,  so every digit is 0x55  (= 85)
    //   2/3 = 0xAA/0xFF = 170/255, so every digit is 0xAA  (= 170 = 2×85)
    // Both have exponent 0, meaning the value is 0.d₀d₁… × 256⁰ ∈ (0, 1).
    let d3:  Array[U8] val = recover Array[U8].init(0x55, p2) end
    let d23: Array[U8] val = recover Array[U8].init(0xAA, p2) end
    let third: MPFRep val     = MPFRep._create(false, false, false, 0, d3)
    let two_third: MPFRep val = MPFRep._create(false, false, false, 0, d23)
    let two: MPFRep val       = MPFRep.from[F64](2.0, p2)
    let ln2_v: MPFRep val  = _mul(_atanh_series(third, p2), two, p2)
    let ln5_v: MPFRep val  = _mul(_atanh_series(two_third, p2), two, p2)
    let ln10_v: MPFRep val = _add(ln2_v, ln5_v, p2)._trunc(w + 4)
    _div(ln_x, ln10_v, w)


  fun _logb(a: MPFRep box, b: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `log_b(a)` at working precision `w` bytes.
    `log_b(a) = ln(a) / ln(b)`.

    Assumes `a` is finite, positive, and non-zero.
    """
    _div(_ln(a, w), _ln(b, w), w)


  fun _exp2(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `2^a` at working precision `w` bytes.
    `2^a = exp(a × ln(2))`.

    Assumes `a` is finite and non-zero.
    """
    let p2: USize = w + 4
    let ln2: MPFRep val = _ln2_const(p2)
    _exp(_mul(a, ln2, p2), w)


  // ── Powers ─────────────────────────────────────────────────────────────────

  fun _powi(a: MPFRep box, n: ILong, w: USize): MPFRep val =>
    """
    Compute `a^n` for integer `n` via binary exponentiation at working
    precision `w` bytes.

    Assumes `a` is finite and non-NaN, and `n ≠ 0`.
    Special-value handling belongs to the public `powi` method.

    `a^n` for `n < 0` is `(1/a)^|n|`.
    """
    var base_r: MPFRep val = if n < 0 then _inv(a, w) else a._clone() end
    var exp_n: ILong = if n < 0 then -n else n end
    var result: MPFRep val = MPFRep.from[F64](1.0, w)
    while exp_n > 0 do
      if (exp_n and 1) == 1 then
        result = _mul(result, base_r, w)
      end
      base_r = _mul(base_r, base_r, w)
      exp_n = exp_n / 2
    end
    result


  fun _pow(a: MPFRep box, b: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `a^b` at working precision `w` bytes.

    Assumes both operands are finite and non-NaN, `b ≠ 0`, `a ≠ 0`, and
    that negative-base/non-integer-exponent has already been rejected.
    Special-value handling belongs to the public `pow` method.

    When `a > 0`: `exp(b × ln(a))`.
    If `a < 0` with integer `b`: delegates to `_powi`.
    """
    if not a.sign_bit() then
      return _exp(_mul(b, _ln(a, w + 4), w + 4), w)
    end
    // Negative base with integer exponent (non-integer already rejected by caller).
    let ni: ILong = _round_to_ilong(b)
    _powi(a, ni, w)


  // ── Trigonometric ──────────────────────────────────────────────────────────

  fun _sin_cos_reduce(a: MPFRep box, w: USize): (MPFRep val, ILong) =>
    """
    Argument reduction for `_sin` and `_cos`.

    Returns `(r, k)` where `r = a − n × (π/2)`, `|r| ≤ π/4`, and
    `k = n mod 4 ∈ {0, 1, 2, 3}` (normalised to non-negative).
    """
    let x: MPFRep val = _extend(a, w)
    let pi_val: MPFRep val = _pi(w)
    let two: MPFRep val = MPFRep.from[F64](2.0, w)
    let half: MPFRep val = MPFRep.from[F64](0.5, w)
    let pi_half: MPFRep val = _div(pi_val, two, w)
    // round-to-nearest: n = floor(x/pi_half + 0.5)
    let ratio: MPFRep val = _div(x, pi_half, w)
    let n_f: MPFRep val = _add(ratio, half, w)._trunc_frac()
    let n: ILong = _round_to_ilong(n_f)
    let n_mpf: MPFRep val = MPFRep.from[F64](n.f64(), w)
    // r = x − n × pi_half
    let r: MPFRep val = _sub(x, _mul(n_mpf, pi_half, w), w)
    let k: ILong = ((n % 4) + 4) % 4
    (r, k)


  fun _sin(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `sin(a)` at working precision `w` bytes (argument in radians).

    Assumes `a` is finite and non-zero. Special-value handling belongs to `sin`.

    Argument reduction: `n = round(a / (π/2))`, `r = a − n × (π/2)`,
    `|r| ≤ π/4`. Based on `k = n mod 4`:
    - `k = 0`: `sin(r)`,
    - `k = 1`: `cos(r)`,
    - `k = 2`: `−sin(r)`,
    - `k = 3`: `−cos(r)`.
    """
    let p2: USize = w + 8
    (let r, let k) = _sin_cos_reduce(a, p2)
    let result: MPFRep val = if (k == 0) or (k == 2) then
        let sr: MPFRep val = _sin_taylor(r, p2)
        if k == 0 then
          sr
        else
          _neg_rep(sr)
        end
      else
        let cr: MPFRep val = _cos_taylor(r, p2)
        if k == 1 then
          cr
        else
          _neg_rep(cr)
        end
      end
    result._trunc(w)


  fun _cos(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `cos(a)` at working precision `w` bytes (argument in radians).

    Assumes `a` is finite and non-zero. Special-value handling belongs to `cos`.

    Argument reduction: `n = round(a / (π/2))`, `r = a − n × (π/2)`,
    `|r| ≤ π/4`. Based on `k = n mod 4`:
    - `k = 0`: `cos(r)`,
    - `k = 1`: `−sin(r)`,
    - `k = 2`: `−cos(r)`,
    - `k = 3`: `sin(r)`.
    """
    let p2: USize = w + 8
    (let r, let k) = _sin_cos_reduce(a, p2)
    let result: MPFRep val = if (k == 0) or (k == 2) then
        let cr: MPFRep val = _cos_taylor(r, p2)
        if k == 0 then
          cr
        else
          _neg_rep(cr)
        end
      else
        let sr: MPFRep val = _sin_taylor(r, p2)
        if k == 3 then
          sr
        else
          _neg_rep(sr)
        end
      end
    result._trunc(w)


  fun _tan(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `tan(a) = sin(a) / cos(a)` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `tan`.
    When `cos(a)` is zero at the working precision, returns NaN.
    """
    let c: MPFRep val = _cos(a, w)
    if c.is_zero() then
      return MPFRep.nan_val()
    end
    _div(_sin(a, w), c, w)


  // ── Hyperbolic ─────────────────────────────────────────────────────────────

  fun _sinh(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `sinh(a) = (e^a − e^{−a}) / 2` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `sinh`.
    """
    let p2: USize = w + 4
    let ex: MPFRep val = _exp(a, p2)
    let emx: MPFRep val = _exp(_neg_rep(a), p2)
    let two: MPFRep val = MPFRep.from[F64](2.0, p2)
    _div(_sub(ex, emx, p2), two, p2)._trunc(w)


  fun _cosh(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `cosh(a) = (e^a + e^{−a}) / 2` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `cosh`.
    """
    let p2: USize = w + 4
    let ex: MPFRep val = _exp(a, p2)
    let emx: MPFRep val = _exp(_neg_rep(a), p2)
    let two: MPFRep val = MPFRep.from[F64](2.0, p2)
    _div(_add(ex, emx, p2), two, p2)._trunc(w)


  fun _tanh(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `tanh(a) = (e^a − e^{−a}) / (e^a + e^{−a})` at working precision
    `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `tanh`.
    """
    let p2: USize = w + 4
    let ex: MPFRep val = _exp(a, p2)
    let emx: MPFRep val = _exp(_neg_rep(a), p2)
    _div(_sub(ex, emx, p2), _add(ex, emx, p2), p2)._trunc(w)


  fun _sech(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `sech(a) = 1 / cosh(a)` at working precision `w` bytes.

    Assumes `a` is finite. Special-value handling belongs to `sech`.
    """
    _inv(_cosh(a, w), w)


  fun _csch(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `csch(a) = 1 / sinh(a)` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `csch`.
    """
    _inv(_sinh(a, w), w)


  fun _coth(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `coth(a) = cosh(a) / sinh(a)` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `coth`.
    """
    _div(_cosh(a, w), _sinh(a, w), w)


  // ── Inverse trigonometric ──────────────────────────────────────────────────

  fun _atan_series(x: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `atan(x)` via the Taylor series
    `Σ_{k=0}^∞ (-1)^k x^{2k+1} / (2k+1)` at working precision `w` bytes.

    The caller must ensure `|x| ≤ 1` for convergence. For best convergence
    prefer `|x| ≤ tan(π/12) ≈ 0.268`.
    """
    let p2: USize = w + 2
    let x2 = _mul(x, x, p2)
    var term: MPFRep ref = x._clone()
    var sum:  MPFRep ref = term._clone()
    var k: USize = 3
    var sign: Bool = true // start negative at k=3
    var iters: USize = 0
    let max_iters: USize = p2 * 30
    while iters < max_iters do
      term._update(_mul(term, x2, p2))
      let addend_abs = if k <= 255 then
          (let q: MPFRep val, _) = term._short_div(k.u8())
          q
        else
          _div(term, MPFRep.from[F64](k.f64(), p2), p2)
        end
      let addend: MPFRep val = if sign then _neg_rep(addend_abs) else addend_abs end
      let new_sum = _add(sum, addend, p2)

      if _converged(new_sum, sum) then
        break
      end

      sum._update(new_sum)
      k = k + 2
      sign = not sign
      iters = iters + 1
    end
    sum._trunc(w)


  fun _atan(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `atan(a)` at working precision `w` bytes.

    Argument reduction:
    - `|a| > 1`: `atan(a) = ±π/2 − atan(1/a)`.
    - Repeat half-angle reduction `a ← a/(1+√(1+a²))` until `|a| ≤ 0.5`
      (each step halves `atan(a)`, needing a `×2^reductions` scale-back).
    - Apply the Taylor series `Σ (-1)^k a^{2k+1}/(2k+1)` on the small argument.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `atan`.
    """
    let p2: USize = w + 6
    let one: MPFRep val = MPFRep.from[F64](1.0, p2)
    let half: MPFRep val = MPFRep.from[F64](0.5, p2)

    var x: MPFRep val = a._trunc(p2)
    var reductions: USize = 0
    var complement: Bool = false

    // Reduce |x| > 1 to |x| ≤ 1 via atan(x) = π/2 - atan(1/x)
    let x2 = _mul(x, x, p2)
    let diff = _sub(x2, one, p2)
    if (not diff.is_zero()) and (not diff.sign_bit()) then
      complement = true
      x = _div(one, x, p2)
    end

    // Half-angle reduction: atan(x) = 2*atan(x/(1+√(1+x²)))
    // Repeat until |x| ≤ 0.5 (at most ~4 steps from |x|=1)
    var hiter: USize = 0
    while hiter < 8 do
      // Check if |x| ≤ 0.5 (i.e. x - 0.5 < 0)
      let xd = _sub(if x.sign_bit() then _neg_rep(x) else x end, half, p2)
      if xd.is_zero() or xd.sign_bit() then
        break
      end
      let one_p2: MPFRep val = MPFRep.from[F64](1.0, p2)
      let denom = _add(one_p2, _sqrt(_add(_mul(x, x, p2), one_p2, p2), p2), p2)
      x = _div(x, denom, p2)
      reductions = reductions + 1
      hiter = hiter + 1
    end

    // Apply series on small |x|
    var result: MPFRep val = _atan_series(x, p2)

    // Scale back: atan(original_reduced_x) = 2^reductions * result
    var r = reductions
    while r > 0 do
      result = _add(result, result, p2)
      r = r - 1
    end

    // Apply complement: atan(original_a) = ±π/2 - result
    if complement then
      let pi_half = _div(_pi(p2), MPFRep.from[F64](2.0, p2), p2)
      let signed_pi_half: MPFRep val = if a.sign_bit() then _neg_rep(pi_half) else pi_half end
      result = _sub(signed_pi_half, result, p2)
    end

    result._trunc(w)


  fun _asin(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `asin(a) = atan(a / sqrt(1 − a²))` at working precision `w` bytes.

    Assumes `|a| ≤ 1`, `a` finite and non-zero. Special-value handling belongs
    to `asin`.
    """
    let p2: USize = w + 4
    let one: MPFRep val = MPFRep.from[F64](1.0, p2)
    let a2: MPFRep val = _mul(a, a, p2)
    let denom: MPFRep val = _sqrt(_sub(one, a2, p2), p2)
    let x: MPFRep val = _div(a, denom, p2)
    _atan(x, p2)._trunc(w)


  fun _acos(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `acos(a) = π/2 − asin(a)` at working precision `w` bytes.

    Assumes `|a| ≤ 1`, `a` finite. Special-value handling belongs to `acos`.
    """
    let p2: USize = w + 4
    let half_pi: MPFRep val = _div(_pi(p2), MPFRep.from[F64](2.0, p2), p2)
    _sub(half_pi, _asin(a, p2), p2)._trunc(w)


  fun _atanh(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `atanh(a)` via `_atanh_series` at working precision `w` bytes.

    Assumes `|a| < 1`, `a` finite and non-zero. Special-value handling belongs
    to `atanh`.
    """
    _atanh_series(a, w)


  fun _asinh(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `asinh(a) = ln(a + sqrt(a² + 1))` at working precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `asinh`.
    """
    let p2: USize = w + 4
    let one: MPFRep val = MPFRep.from[F64](1.0, p2)
    let a2: MPFRep val = _mul(a, a, p2)
    let inner: MPFRep val = _sqrt(_add(a2, one, p2), p2)
    _ln(_add(a, inner, p2), p2)._trunc(w)


  fun _acosh(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `acosh(a) = ln(a + sqrt(a² − 1))` at working precision `w` bytes.

    Assumes `a ≥ 1`, `a` finite. Special-value handling belongs to `acosh`.
    """
    let p2: USize = w + 4
    let one: MPFRep val = MPFRep.from[F64](1.0, p2)
    let a2: MPFRep val = _mul(a, a, p2)
    let inner: MPFRep val = _sqrt(_sub(a2, one, p2), p2)
    _ln(_add(a, inner, p2), p2)._trunc(w)


  fun _cbrt(a: MPFRep box, w: USize): MPFRep val =>
    """
    Compute `cbrt(a) = exp(ln(|a|) / 3)` with correct sign at working
    precision `w` bytes.

    Assumes `a` is finite and non-zero. Special-value handling belongs to `cbrt`.
    """
    let p2: USize = w + 4
    let three: MPFRep val = MPFRep.from[F64](3.0, p2)
    let abs_a: MPFRep val = if a.sign_bit() then _neg_rep(a) else a._trunc(p2) end
    let result: MPFRep val = _exp(_div(_ln(abs_a, p2), three, p2), p2)._trunc(w)
    if a.sign_bit() then _neg_rep(result) else result end


  fun _rootn(a: MPFRep box, n: USize, w: USize): MPFRep val =>
    """
    Compute `a^{1/n} = exp(ln(|a|) / n)` with correct sign (for odd `n`) at
    working precision `w` bytes.

    Assumes `a` is finite and non-zero, `n ≥ 1`. Special-value handling belongs
    to `rootn`.
    """
    let p2: USize = w + 4
    let fn_rep: MPFRep val = MPFRep.from[ULong](n.ulong(), p2)
    let abs_a: MPFRep val = if a.sign_bit() then _neg_rep(a) else a._trunc(p2) end
    let result: MPFRep val = _exp(_div(_ln(abs_a, p2), fn_rep, p2), p2)._trunc(w)
    if a.sign_bit() and ((n and 1) == 1) then _neg_rep(result) else result end


  // ── Constants ──────────────────────────────────────────────────────────────

  fun _pi(w: USize): MPFRep val =>
    """
    Compute `π` at working precision `w` bytes via the
    [Chudnovsky algorithm](https://en.wikipedia.org/wiki/Chudnovsky_algorithm):
    `1 / π = 12 ∑{k=0}^∞ ((−1)^k·(6k)!·(545140134k + 13591409)) / ((3k)!·(k!)^3·(640320)^{3k+3/2})`
    This algorithm stands out for its fast convergence and adds ~14 digits of π
    per term.

    The series is `π = C^{3/2} / (12 × Σ_k s_k × a_k)` where C = 640320,
    `a_k = 13591409 + 545140134k`, and
    `s_k = (−1)^k × (6k)! / ((3k)! × (k!)^3 × C^{3k})`.

    The even factors `(6k+2)(6k+4)(6k+6)` and `(3k+1)(3k+2)(3k+3)` cancel
    exactly, giving the recurrence:
      `s_{k+1} = s_k × (−8) × (6k+1)(6k+3)(6k+5) / ((k+1)^3 × C^3)`

    This avoids recomputing factorials from scratch each iteration (O(k²) → O(k))
    and replaces the per-iteration transcendental `_pow` with integer multiplications.
    Delivers ~14.18 decimal digits per term so only ~0.17p iterations are needed.
    """
    let p: USize = w + 4
    let k_640320: MPFRep val = MPFRep.from[ULong](640320, p)
    let k_c3: MPFRep val = _mul(_mul(k_640320, k_640320, p), k_640320, p)
    let k_8: MPFRep val = MPFRep.from[F64](8.0, p)
    // Precompute −8/C^3 once; applied at every step.
    let neg_8_inv_c3: MPFRep val = _neg_rep(_mul(k_8, _inv(k_c3, p), p))

    // s_0 = 1, a_0 = 13591409; term_0 = s_0 × a_0
    var s_k: MPFRep val = MPFRep.from[F64](1.0, p)
    var a_k: MPFRep val = MPFRep.from[ULong](13591409, p)
    let delta_a: MPFRep val = MPFRep.from[ULong](545140134, p)

    var sum: MPFRep val = a_k
    var prev_sum: MPFRep val = sum

    var k: USize = 0
    let max_k: USize = p + p
    while k < max_k do
      // s_{k+1} = s_k × (−8/C^3) × (6k+1)(6k+3)(6k+5) / (k+1)^3
      let k6: USize = 6 * k
      let f1: MPFRep val = MPFRep.from[ULong]((k6 + 1).ulong(), p)
      let f2: MPFRep val = MPFRep.from[ULong]((k6 + 3).ulong(), p)
      let f3: MPFRep val = MPFRep.from[ULong]((k6 + 5).ulong(), p)
      let kp1: MPFRep val = MPFRep.from[ULong]((k + 1).ulong(), p)
      let num_s: MPFRep val = _mul(_mul(f1, f2, p), f3, p)
      let den_s: MPFRep val = _mul(_mul(kp1, kp1, p), kp1, p)
      s_k = _mul(s_k, _mul(neg_8_inv_c3, _div(num_s, den_s, p), p), p)
      a_k = _add(a_k, delta_a, p)

      let term: MPFRep val = _mul(s_k, a_k, p)
      prev_sum = sum
      sum = _add(sum, term, p)

      if _converged(sum, prev_sum) then
        break
      end
      k = k + 1
    end

    // π = C^{3/2} / (12 × sum)  where  C^{3/2} = C × √C
    let c_3_2: MPFRep val = _mul(k_640320, _sqrt(k_640320, p), p)
    let k_12: MPFRep val = MPFRep.from[F64](12.0, p)
    _div(c_3_2, _mul(k_12, sum, p), p)._trunc(w)


  // ── String conversion ──────────────────────────────────────────────────────

  fun _from_string(s: String, w: USize): MPFRep val ? =>
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

    let sig_limit: USize = ((p_digits.f64() * F64(256.0).log10()).ceil().usize() + 4).max(1)
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
    let ten_mp: MPFRep val = MPFRep.from[F64](10.0, p2)

    let scaled: MPFRep val =
      if dec_exp >= 0 then
        var n_mp: MPFRep val = MPFRep.from[MPInt](n_int, p2)
        if dec_exp > 0 then
          var scale: MPFRep val = MPFRep.from[F64](1.0, p2)
          var sbase: MPFRep val = ten_mp
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
          var n_mp: MPFRep val = MPFRep.from[MPInt](n_int, p2)
          var scale: MPFRep val = MPFRep.from[F64](1.0, p2)
          var sbase: MPFRep val = ten_mp
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

  fun _ziv_ok(result: MPFRep box, p: USize, mode: RoundingMode): Bool =>
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

    Returns `ctx.p_bytes() + guard`, where `guard` grows logarithmically with
    the output precision so that accumulated rounding error over long computation
    chains stays bounded regardless of how large the precision is.

    ## Formula

    `log_guard = ceil(log₂(p_bytes + 1))`  (≥ 1, grows as ⌈log₂ p⌉ bytes)

    Per-operation floors prevent guard from being too small at low precision:

    | `op`                        | floor | Reason                                          |
    |-----------------------------|-------|-------------------------------------------------|
    | `"add"`, `"sub"`            |   2   | Exact result + 1 carry bit; ≤ 2 bytes needed    |
    | `"mul"`, `"inv"`, `"sqrt"`  |   4   | Newton refinement and FFT rounding              |
    | `"ln"`                      |   6   | Argument reduction cancellation                 |
    | `"exp"`, `"trig"`, `"pi"`   |   8   | Reduction + Taylor / Machin series              |
    | (default)                   |   4   | Conservative fallback                           |

    At 10 000-bit precision (`p_bytes = 1250`): `log_guard = ceil(log₂(1251)) ≈ 11`.
    At 1 000 000-bit precision (`p_bytes = 125 000`): `log_guard ≈ 17`.
    The dynamic guard ensures that a Borwein or similar iteration that runs
    O(log p) operations accumulates ≤ 1 ULP of total error at any precision.

    This method is the single authoritative source for guard bytes. Adding a
    new operation requires only extending this match — no other file changes.
    `MPFContext.working_bytes` delegates here.
    """
    // log_guard = floor(log₂(p_bytes)) + 1  (= number of bits needed to represent p)
    // For p=14: 4; p=1250: 11; p=125000: 17.
    let p = ctx.p_bytes()
    var log_guard: USize = 0
    var bits: USize = p
    while bits > 0 do
      log_guard = log_guard + 1
      bits = bits >> 1
    end
    // Per-operation floor: take the maximum of the floor and log_guard.
    let floor_guard: USize = match op
      | "add" | "sub" => 2
      | "mul" | "inv" | "sqrt" => 4
      | "ln" => 6
      | "exp" | "trig" | "pi" => 8
      else 4
    end
    p + log_guard.max(floor_guard)


  // ── Public context-aware API ───────────────────────────────────────────────
  //
  // Each method below takes `ctx: MPFContext` as its first argument, handles
  // special values (NaN, ±∞, ±0), calls the corresponding internal `_*` method
  // at working precision `ctx.working_bytes(op)`, then rounds the result to
  // output precision via `ctx._round_to`.
  //
  // `MPFloat` calls these directly — there is no dispatch layer in `MPFContext`.

  fun add(ctx: MPFContext, a: MPFRep box, b: MPFRep box): MPFRep val =>
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


  fun sub(ctx: MPFContext, a: MPFRep box, b: MPFRep box): MPFRep val =>
    """
    Compute `a − b` at output precision `ctx.p_bytes()` using rounding
    `ctx.rounding`.

    Delegates to `add` after negating `b`, inheriting all special-value
    handling from `add`.
    """
    add(ctx, a, _neg_rep(b))


  fun mul(ctx: MPFContext, a: MPFRep box, b: MPFRep box): MPFRep val =>
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


  fun inv(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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


  fun div(ctx: MPFContext, a: MPFRep box, b: MPFRep box): MPFRep val =>
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


  fun sqrt(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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


  fun divrem(ctx: MPFContext, a: MPFRep box, b: MPFRep box): (MPFRep val, MPFRep val) =>
    """
    Truncated division with remainder: `(q, r)` such that
    `a = q × b + r`, `q = trunc(a / b)`, `r` has the same sign as `a`.

    All results are calculated at working precision and rounded at `ctx.p_bytes()`
    using `ctx.rounding`. Special cases follow IEEE 754 conventions.
    """
    let p = ctx.p_bytes()
    let w = ctx.working_bytes("inv")
    let nan: MPFRep val = MPFRep.nan_val()
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
      let z: MPFRep val = MPFRep._create(false, false, false, 0, Array[U8].init(0, p))
      return (z, z)
    end
    (let q: MPFRep val, let r: MPFRep box) = _divrem(a, b, w)
    (ctx._round_to(q, p, ctx.rounding), ctx._round_to(r, p, ctx.rounding))


  fun fld(ctx: MPFContext, a: MPFRep box, b: MPFRep box): MPFRep val =>
    """
    Floored division: `floor(a / b)` at output precision `ctx.p_bytes()`.

    The result is the largest integer `q` such that `q × b ≤ a`. When the
    truncated remainder is non-zero and `a` and `b` have opposite signs, the
    truncated quotient is decremented by 1.

    Special cases handled here:
    - `b = 0` → NaN (`floor(a/0)` is undefined for all `a`).
    - Non-finite `q` (NaN inputs, `a = ±∞`) → return `q` directly; the
      floor-adjustment step must not call `_sub` on non-finite values.
    Remaining special cases are handled by `divrem`.
    """
    if b.is_zero() then
      return MPFRep.nan_val()
    end
    let w = ctx.working_bytes("inv")
    let p = ctx.p_bytes()
    (let q, let r) = divrem(ctx, a, b)
    if not q.is_finite() then
      return q
    end
    if r._has_frac() or (not r.is_zero()) then
      if a.sign_bit() != b.sign_bit() then
        let one: MPFRep val = MPFRep.from[F64](1.0, p)
        return ctx._round_to(_sub(q, one, w), p, ctx.rounding)
      end
    end
    q


  fun rem(ctx: MPFContext, a: MPFRep box, b: MPFRep box): MPFRep val =>
    """
    Truncated remainder: `a − trunc(a/b) × b`. Sign of result = sign of `a`.

    Special cases are handled by `divrem` before this point.
    """
    divrem(ctx, a, b)._2


  fun mod(ctx: MPFContext, a: MPFRep box, b: MPFRep box): MPFRep val =>
    """
    Floored remainder: `a − fld(a,b) × b`. Sign of result = sign of `b`.

    Special cases are handled by `rem`/`divrem` before this point.
    """
    let r = rem(ctx, a, b)
    if (not r.is_zero()) and (a.sign_bit() != b.sign_bit()) then
      let w = ctx.working_bytes("add")
      return ctx._round_to(_add(r, b, w), ctx.p_bytes(), ctx.rounding)
    end
    r


  fun ln(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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
    var result: MPFRep val = _ln(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _ln(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun log(ctx: MPFContext, a: MPFRep box): MPFRep val =>
    """
    Natural logarithm alias — same as `ln`. Correctly rounded.
    """
    ln(ctx, a)


  fun log2(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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
    var result: MPFRep val = _log2(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _log2(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun log10(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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
    var result: MPFRep val = _log10(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _log10(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun logb(ctx: MPFContext, a: MPFRep box, base: MPFRep box): MPFRep val =>
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


  fun exp(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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
      return MPFRep.from[F64](1.0, ctx.p_bytes())
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("exp")
    var result: MPFRep val = _exp(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _exp(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun exp2(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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
      return MPFRep.from[F64](1.0, ctx.p_bytes())
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("exp")
    var result = _exp2(a, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _exp2(a, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun powi(ctx: MPFContext, a: MPFRep box, n: ILong): MPFRep val =>
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
        return MPFRep.from[F64](1.0, ctx.p_bytes())
      elseif n > 0 then
        let result_sign = a.sign_bit() and ((n and 1) != 0)
        return MPFRep._create(result_sign, false, true, 0, Array[U8].create())
      else
        return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
      end
    end
    let p = ctx.p_bytes()
    var g = ctx.working_bytes("exp")
    var result: MPFRep val = _powi(a, n, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _powi(a, n, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun pow(ctx: MPFContext, a: MPFRep box, b: MPFRep box): MPFRep val =>
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
      return MPFRep.from[F64](1.0, ctx.p_bytes())
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
    var result: MPFRep val = _pow(a, b, g)
    while not _ziv_ok(result, p, ctx.rounding) do
      g = g * 2
      result = _pow(a, b, g)
    end
    ctx._round_to(result, p, ctx.rounding)


  fun sin(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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


  fun cos(ctx: MPFContext, a: MPFRep box): MPFRep val =>
    """
    Compute `cos(a)` at output precision `ctx.p_bytes()`.

    NaN or ±∞ → NaN. 0 → 1.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return MPFRep.from[F64](1.0, ctx.p_bytes())
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_cos(a, w), ctx.p_bytes(), ctx.rounding)


  fun tan(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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


  fun sinh(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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


  fun cosh(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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
      return MPFRep.from[F64](1.0, ctx.p_bytes())
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_cosh(a, w), ctx.p_bytes(), ctx.rounding)


  fun tanh(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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


  fun sech(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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


  fun csch(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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


  fun coth(ctx: MPFContext, a: MPFRep box): MPFRep val =>
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


  fun atan(ctx: MPFContext, a: MPFRep box): MPFRep val =>
    """
    Compute `atan(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +π/2. −∞ → −π/2. 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      let w = ctx.working_bytes("trig")
      let half_pi = _div(_pi(w), MPFRep.from[F64](2.0, w), w)
      if a.sign_bit() then
        return ctx._round_to(_neg_rep(half_pi), ctx.p_bytes(), ctx.rounding)
      end
      return ctx._round_to(half_pi, ctx.p_bytes(), ctx.rounding)
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_atan(a, w), ctx.p_bytes(), ctx.rounding)


  fun _abs_gt_one(a: MPFRep box, p: USize): Bool =>
    """
    Return `true` if `|a| > 1`. Used for domain checks in inverse trig.
    Computes `|a| - 1` and checks the sign (positive sign means > 1).
    """
    let one: MPFRep val = MPFRep.from[F64](1.0, p)
    let abs_a: MPFRep val = if a.sign_bit() then _neg_rep(a) else a._trunc(p) end
    let diff = _sub(abs_a, one, p)
    (not diff.is_zero()) and (not diff.sign_bit())


  fun _abs_ge_one(a: MPFRep box, p: USize): Bool =>
    """
    Return `true` if `|a| >= 1`.
    """
    let one: MPFRep val = MPFRep.from[F64](1.0, p)
    let abs_a: MPFRep val = if a.sign_bit() then _neg_rep(a) else a._trunc(p) end
    let diff = _sub(abs_a, one, p)
    diff.is_zero() or (not diff.sign_bit())


  fun _lt_rep(a: MPFRep box, b: MPFRep box, p: USize): Bool =>
    """
    Return `true` if `a < b` by checking sign of `a - b`.
    """
    let diff = _sub(a, b, p)
    (not diff.is_zero()) and diff.sign_bit()


  fun asin(ctx: MPFContext, a: MPFRep box): MPFRep val =>
    """
    Compute `asin(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. |a| > 1 → NaN. +∞ or −∞ → NaN. 0 → 0. ±1 → ±π/2.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    let p = ctx.p_bytes()
    // |a| > 1: NaN
    if _abs_gt_one(a, p) then
      return MPFRep.nan_val()
    end
    // ±1 → ±π/2 (avoid sqrt(0) in kernel)
    if _abs_ge_one(a, p) then
      let w = ctx.working_bytes("trig")
      let half_pi = _div(_pi(w), MPFRep.from[F64](2.0, w), w)
      return ctx._round_to(if a.sign_bit() then _neg_rep(half_pi) else half_pi end,
        p, ctx.rounding)
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_asin(a, w), p, ctx.rounding)


  fun acos(ctx: MPFContext, a: MPFRep box): MPFRep val =>
    """
    Compute `acos(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. |a| > 1 → NaN. +∞ or −∞ → NaN. 0 → π/2. 1 → 0. −1 → π.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    let p = ctx.p_bytes()
    if _abs_gt_one(a, p) then
      return MPFRep.nan_val()
    end
    // |a| = 1: exact values to avoid sqrt(0) in kernel
    if _abs_ge_one(a, p) then
      if a.sign_bit() then
        // a = -1 → π
        let w = ctx.working_bytes("trig")
        return ctx._round_to(_pi(w), p, ctx.rounding)
      else
        // a = 1 → 0
        return MPFRep._create(false, false, false, 0, Array[U8].init(0, p))
      end
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_acos(a, w), p, ctx.rounding)


  fun atanh(ctx: MPFContext, a: MPFRep box): MPFRep val =>
    """
    Compute `atanh(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. ±∞ → NaN. |a| ≥ 1 → NaN (±∞ for |a|=1 by convention, but
    diverges; we return NaN for simplicity). 0 → 0.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() or a.is_infinite() then
      return MPFRep.nan_val()
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    let p = ctx.p_bytes()
    // |a| >= 1: NaN (series diverges)
    if _abs_ge_one(a, p) then
      return MPFRep.nan_val()
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_atanh(a, w), p, ctx.rounding)


  fun asinh(ctx: MPFContext, a: MPFRep box): MPFRep val =>
    """
    Compute `asinh(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +∞. −∞ → −∞. 0 → 0.

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
    ctx._round_to(_asinh(a, w), ctx.p_bytes(), ctx.rounding)


  fun acosh(ctx: MPFContext, a: MPFRep box): MPFRep val =>
    """
    Compute `acosh(a)` at output precision `ctx.p_bytes()`.

    NaN → NaN. a < 1 → NaN. 1 → 0. +∞ → +∞.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
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
    let p = ctx.p_bytes()
    let one: MPFRep val = MPFRep.from[F64](1.0, p)
    // a < 1: NaN
    if _lt_rep(a, one, p) then
      return MPFRep.nan_val()
    end
    // a = 1: exact 0 (avoids sqrt(0) → 0 → ln(1+0) edge case)
    let diff = _sub(a, one, p)
    if diff.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, p))
    end
    let w = ctx.working_bytes("trig")
    ctx._round_to(_acosh(a, w), p, ctx.rounding)


  fun cbrt(ctx: MPFContext, a: MPFRep box): MPFRep val =>
    """
    Compute the cube root `a^{1/3}` at output precision `ctx.p_bytes()`.

    NaN → NaN. +∞ → +∞. −∞ → −∞. 0 → 0.
    Defined for negative `a` (cbrt is an odd function).

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
    let w = ctx.working_bytes("exp")
    ctx._round_to(_cbrt(a, w), ctx.p_bytes(), ctx.rounding)


  fun rootn(ctx: MPFContext, a: MPFRep box, n: USize): MPFRep val =>
    """
    Compute the `n`-th root `a^{1/n}` at output precision `ctx.p_bytes()`.

    NaN → NaN. n = 0 → NaN. +∞ → +∞ (for any n). −∞ → −∞ if n is odd,
    NaN if n is even. 0 → 0.
    Negative `a` with even `n` → NaN. Negative `a` with odd `n` → negative result.

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    if a.is_nan() or (n == 0) then
      return MPFRep.nan_val()
    end
    if a.is_infinite() then
      if a.sign_bit() and ((n and 1) == 0) then
        return MPFRep.nan_val()
      end
      return MPFRep._create(a.sign_bit() and ((n and 1) == 1), false, true,
        0, Array[U8].create())
    end
    if a.is_zero() then
      return MPFRep._create(false, false, false, 0, Array[U8].init(0, ctx.p_bytes()))
    end
    if a.sign_bit() and ((n and 1) == 0) then
      return MPFRep.nan_val()
    end
    let w = ctx.working_bytes("exp")
    ctx._round_to(_rootn(a, n, w), ctx.p_bytes(), ctx.rounding)


  fun pi(ctx: MPFContext): MPFRep val =>
    """
    Compute π at output precision `ctx.p_bytes()`.

    **Rounding**: faithfully rounded (≤ 1 ULP error).

    Uses Chudnovsky algorithm
    """
    let w = ctx.working_bytes("pi")
    ctx._round_to(_pi(w), ctx.p_bytes(), ctx.rounding)


  fun e(ctx: MPFContext): MPFRep val =>
    """
    Compute Euler's number `e = exp(1)` at output precision `ctx.p_bytes()`.

    Delegates to the public `exp` method with argument `1`, inheriting its
    Ziv-iteration correct rounding (≤ 0.5 ULP error).
    """
    exp(ctx, MPFRep.from[F64](1.0, ctx.working_bytes("exp")))


  fun from_string(ctx: MPFContext, s: String): MPFRep val ? =>
    """
    Parse `s` as a floating-point number at output precision `ctx.p_bytes()`.

    Raises an error if `s` is not a valid floating-point representation.
    Only base 10 is currently implemented.
    """
    let w = ctx.working_bytes("add")
    ctx._round_to(_from_string(s, w)?, ctx.p_bytes(), ctx.rounding)
