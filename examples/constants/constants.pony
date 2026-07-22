// Calculation of some constants

use "collections"

use "../../mathx"
use "../../bignum"
use gmp = "../../bignum/gmp"

use @mpfr_snprintf[I32](buf: Pointer[U8] tag, buf_size: U64, fmt: Pointer[U8] tag, ...)
use @pony_ctx[Pointer[None]]()
use @pony_alloc[Pointer[U8]](ctx: Pointer[None], size: USize)


actor Main
  """
  Calculates various mathematical constants using multiple precision floating point numbers.
  """
  new create(env: Env) =>
    """
    Calculate a few *well-known* mathematical constants using big floats.
    """
    env.out.print("Calculate some well-known constants using formulas and" +
                  " infinite-precision floats.")
    env.out.print("")
    env.out.print("Pi is calculated with 10_000 decimals. e calculation uses 1000 decimals")
    env.out.print("")
    e(env)
    pi_mpfr(env)
    pi_mpf(env)
    pi_gmpfloat(env)
    pi_mpfloat(env)
    pi_compare(env)


  fun pi_mpfr(env: Env) =>
    """
    The Pi constant calculated by MPFR
    """
    env.out.print("Using MPFR")
    env.out.print("----------")

    let pi_const = gmp.MPFloat.pi(10000)
    env.out.print("Pi constant from MPFR =")
    env.out.print(pi_const.string())
    env.out.print("")


  fun pi_mpfloat(env: Env) =>
    """
    Calculate Pi = 3.14149... ($$\pi$$) using the pure pony `MPFloat` class.

    Translation of the mathematical formulas is almost immediate.

    This algorithm uses [Borwein quadratic](https://en.wikipedia.org/wiki/Borwein%27s_algorithm#Quadratic_convergence_(1984))
    formula:

    Use the serie:
      * x_0 = sqrt(2)
      * y_0 = 0
      * pi_0 = 2 + sqrt(2)
    
      * x_n+1 = (sqrt(x_n) + 1/sqrt(x_n)) / 2
      * y_n+1 = (1 + y_n) * sqrt(x_n) / (x_n + y_n)
      * pi_n+1 = pi_n * (1 + x_n+1) * y_n+1 / (1 + y_n+1)
    """
    var i: USize = 0
    let precision: USize = 10000 // binary digits for the calculations, and result
    let zero = MPFloat.from[F64](0.0, precision)
    let one = MPFloat.from[F64](1.0, precision)
    let two = MPFloat.from[F64](2.0, precision)
    let two_sqrt = two.sqrt()

    let x_0 = two_sqrt
    let y_0 = zero
    let pi_0 = two + two_sqrt

    var x_n = x_0
    var y_n = y_0
    var pi_n = pi_0
    var old_pi = pi_0
    repeat
      old_pi = pi_n
      let x_n1 = (x_n.sqrt() + (one / x_n.sqrt())) / two
      let y_n1 = ((one + y_n) * x_n.sqrt()) / (x_n + y_n)
      let pi_n1 = (pi_n * (one + x_n1) * y_n1) / (one + y_n1)

      // Loop to next n
      x_n = x_n1
      y_n = y_n1
      pi_n = pi_n1
      i = i + 1
      // We stop if i > precision: there's probably a bug
    until (old_pi == pi_n) or (i > precision) end

    // Show result
    env.out.print("Using MPFloat")
    env.out.print("-------------")
    env.out.print("pi (after " + i.string() + " loops) = ")
    env.out.print(pi_n.string())
    env.out.print("")

    // Compare with pi constant
    env.out.print("Compare result with MPFloat constant")
    let pi_const = MPFloat.pi(precision)
    if pi_const == pi_n then
      env.out.print("Calculation is exact")
    else
      env.out.print("Calculation is inexact")
      let delta = pi_n - pi_const
      //env.out.print("Delta = " + delta.string())
    end
    env.out.print("")


  fun pi_gmpfloat(env: Env) =>
    """
    Calculate Pi = 3.14149... ($$\pi$$) using the GMP `MPFloat` class wrapper.
    The code is the same as the one of `pi_mpfloat`.

    Translation of the mathematical formulas is almost immediate. From a
    performance point of view, this version allocates more memory than
    `pi_mpf`.

    This algorithm uses [Borwein quadratic](https://en.wikipedia.org/wiki/Borwein%27s_algorithm#Quadratic_convergence_(1984))
    formula:

    Use the serie:
      * x_0 = sqrt(2)
      * y_0 = 0
      * pi_0 = 2 + sqrt(2)
    
      * x_n+1 = (sqrt(x_n) + 1/sqrt(x_n)) / 2
      * y_n+1 = (1 + y_n) * sqrt(x_n) / (x_n + y_n)
      * pi_n+1 = pi_n * (1 + x_n+1) * y_n+1 / (1 + y_n+1)
    """
    var i: USize = 0
    let precision: USize = 10000 // binary digits for the calculations, and result
    let zero = gmp.MPFloat.from[F64](0.0, precision)
    let one = gmp.MPFloat.from[F64](1.0, precision)
    let two = gmp.MPFloat.from[F64](2.0, precision)
    let two_sqrt = two.sqrt()

    let x_0 = two_sqrt
    let y_0 = zero
    let pi_0 = two + two_sqrt

    var x_n = x_0
    var y_n = y_0
    var pi_n = pi_0
    var old_pi = pi_0
    repeat
      old_pi = pi_n
      let x_n1 = (x_n.sqrt() + (one / x_n.sqrt())) / two
      let y_n1 = ((one + y_n) * x_n.sqrt()) / (x_n + y_n)
      let pi_n1 = (pi_n * (one + x_n1) * y_n1) / (one + y_n1)

      // Loop to next n
      x_n = x_n1
      y_n = y_n1
      pi_n = pi_n1
      i = i + 1
      // We stop if i > precision: there's probably a bug
    until (old_pi == pi_n) or (i > precision) end

    // Show result
    env.out.print("Using GMP MPFloat")
    env.out.print("-----------------")
    env.out.print("pi (after " + i.string() + " loops) = ")
    env.out.print(pi_n.string())
    env.out.print("")

    // Compare with pi constant
    env.out.print("Compare result with MPFloat constant")
    let pi_const = MPFloat.pi(precision.usize())
    if pi_const.string() == pi_n.string() then
      env.out.print("Calculation is exact")
    else
      env.out.print("Calculation is inexact")
      //let delta = pi_n - pi_const
      //env.out.print("Delta = " + delta.string())
    end
    env.out.print("")


  fun pi_mpf(env: Env) =>
    """
    This version of the calculation of Pi = 3.14159 ($$\pi$$) uses also the
    Borwein's quadratic serie but implement it using the `MPF` class that is
    a raw binding to [MPF](https://www.mpfr.org/) library. This class limits
    memory use by reusing variables, but the developer has to manage memory
    allocations and deallocations. Compared with `pi` implementation, this
    version feels less *natural* to mathematicians and is much longer as
    each line of code is a step of calculation in a formula. It looks more
    like machine code...

    Use the serie:
      * x_0 = sqrt(2)
      * y_0 = 0
      * pi_0 = 2 + sqrt(2)
    
      * x_n+1 = (sqrt(x_n) + 1/sqrt(x_n)) / 2
      * y_n+1 = (1 + y_n) * sqrt(x_n) / (x_n + y_n)
      * pi_n+1 = pi_n * (1 + x_n+1) * y_n+1 / (1 + y_n+1)
    """
    let prec: ILong = 10000 // Precision of calculation
    var i: ILong = 0
    var x: gmp.SMPFr = gmp.SMPFr
    var y: gmp.SMPFr = gmp.SMPFr
    var p: gmp.SMPFr = gmp.SMPFr
    var d: gmp.SMPFr = gmp.SMPFr
    var t: gmp.SMPFr = gmp.SMPFr
    var s: gmp.SMPFr = gmp.SMPFr
    var u: gmp.SMPFr = gmp.SMPFr

    // x_0 = sqrt(2)
    gmp.MPF.init2(x, prec)
    gmp.MPF.set_d(x, 2.0, RoundingNearest)
    gmp.MPF.sqrt(x, x, RoundingNearest)

    // y_0 = 0
    gmp.MPF.init2(y, prec)
    gmp.MPF.set_d(y, 0.0, RoundingNearest)

    // pi_0 = 2 + sqrt(2)
    gmp.MPF.init2(p, prec)
    gmp.MPF.set_d(p, 2.0, RoundingNearest)
    gmp.MPF.add(p, p, x, RoundingNearest)

    // Temporaries for calculations
    gmp.MPF.init2(d, prec) // Previous pi value
    gmp.MPF.init2(t, prec)
    gmp.MPF.init2(s, prec)
    gmp.MPF.init2(u, prec)

    repeat
      // Keep previous value of pi to determine when to stop
      gmp.MPF.set(d, p, RoundingNearest)

      //   x_n+1 = (sqrt(x_n) + 1/sqrt(x_n)) / 2
      gmp.MPF.sqrt(t, x, RoundingNearest)             // t = sqrt(x_n)
      gmp.MPF.d_div(s, 1.0, t, RoundingNearest)       // s = 1/sqrt(x_n)
      gmp.MPF.add(s, t, s, RoundingNearest)           // s = sqrt(x_n) + 1/sqrt(x_n)
      gmp.MPF.div_d(s, s, 2.0, RoundingNearest)       // s = s / 2 --> x_n+1

      //   y_n+1 = (1 + y_n) * sqrt(x_n) / (x_n + y_n)
      gmp.MPF.add_d(u, y, 1.0, RoundingNearest)       // u = 1 + y_n
      gmp.MPF.mul(u, u, t, RoundingNearest)           // u = (1 + y_n) * sqrt(x_n)
      gmp.MPF.add(t, x, y, RoundingNearest)           // t = x_n + y_n
      gmp.MPF.div(t, u, t, RoundingNearest)           // t = u / t --> y_n+1

      //   pi_n+1 = pi_n * (1 + x_n+1) * y_n+1 / (1 + y_n+1)
      gmp.MPF.add_d(u, s, 1.0, RoundingNearest)       // u = 1 + x_n+1
      gmp.MPF.mul(p, p, u, RoundingNearest)           // p = pi_n * (1 + x_n+1)
      gmp.MPF.mul(p, p, t, RoundingNearest)           // p = pi_n * (1 + x_n+1) * y_n+1
      gmp.MPF.add_d(u, t, 1.0, RoundingNearest)       // u = 1 + y_n+1
      gmp.MPF.div(p, p, u, RoundingNearest)           // p = p / u --> pi_n+1
      
      // Now n -> n+1
      gmp.MPF.set(x, s, RoundingNearest)              // x_n+1 = s
      gmp.MPF.set(y, t, RoundingNearest)              // y_n+1 = t

      // Stop when new value has not changed compared with old value
      i = i + 1
      // We stop if i > prec: there's a bug
    until gmp.MPF.equal_p(p, d) or (i > prec) end

    env.out.print("Using MPF")
    env.out.print("---------")
    env.out.print("Pi (after " + i.string() + " loops) = ")
    env.out.print(_string(p))
    env.out.print("")

    // Compare result with constant
    env.out.print("Compare result with MPF constant")
    gmp.MPF.const_pi(d, RoundingNearest)
    if gmp.MPF.equal_p(d, p) then
      env.out.print("Calculation is exact")
    else
      env.out.print("Calculation is inexact")
      gmp.MPF.sub(d, p, d, RoundingNearest)
      //env.out.print("Delta = " + _string(d))
    end
    env.out.print("")

    // Clean memory and cache
    gmp.MPF.clear(x)
    gmp.MPF.clear(y)
    gmp.MPF.clear(p)
    gmp.MPF.clear(d)
    gmp.MPF.clear(t)
    gmp.MPF.clear(s)
    gmp.MPF.clear(u)
    gmp.MPF.free_cache()


  fun _string(x: gmp.SMPFr): String iso^ =>
    """
    Print the result using MPFR's formatting with full precision as the `MPF`
    does not have an output function yet..

    The format string is built dynamically from the number's precision (in bits),
    converted to decimal digits via ceil(prec * log2 / log10). Without this,
    a bare `"%.RNe"` format would have precision=0 (the `.` with no digit count
    before `R` means 0 significant digits), producing output like `"3e+00"`.
    """
    let prec = gmp.MPF.get_prec(x)
    let p = ((prec.f64() * F64(2.0).log()) / F64(10.0).log()).ceil().i32()
    let fmt = ("%." + p.string() + "RNe").cstring()
    let buf_size: I32 = @mpfr_snprintf(Pointer[U8], 0, fmt, x) + 1
    let result: String iso = recover iso
      let buffer: Pointer[U8] = @pony_alloc(@pony_ctx(), buf_size.usize())
      @mpfr_snprintf(buffer, buf_size.u64(), fmt, x)
      String.copy_cstring(buffer)
    end
    consume result


  fun _window(s: String, agree: USize, tail: USize): String =>
    """
    Extract a `tail`-wide window from `s` anchored at the first-divergence
    position `agree`.

    The window covers positions `[start .. start+tail)` where
    `start = agree - (tail - lookahead)` and `lookahead = tail / 5` (≈ 20%
    of the window after the split point, so context before and after are both
    visible).  If `s` is shorter than `start + tail`, the window is
    right-padded with spaces so every caller always returns exactly `tail`
    characters, keeping all lines the same width.
    """
    let lookahead: USize = tail / 5
    let start: USize =
      if agree >= (tail - lookahead) then agree - (tail - lookahead) else 0 end
    let result: String iso = recover
      let r = String.create(tail)
      var i: USize = start
      while i < (start + tail) do
        r.push(try s(i)? else ' ' end)
        i = i + 1
      end
      r
    end
    consume result


  fun _agree(a: String box, b: String box): USize =>
    """Count how many characters agree between `a` and `b` (from position 0)."""
    var i: USize = 0
    while i < a.size().min(b.size()) do
      if (try a(i)? else 0 end) != (try b(i)? else 0 end) then
        break
      end
      i = i + 1
    end
    i


  fun pi_compare(env: Env) =>
    """
    Run the Borwein quadratic iteration in lockstep for both `MPFloat` and
    `gmp.MPFloat` at 10 000-bit precision, printing the last 35 significant
    digits of each intermediate per iteration, right-aligned so differences
    are visible at a glance.  An `agree=N` count shows how many leading
    characters the two results share.
    """
    env.out.print("=== Borwein step-by-step comparison: MPFloat vs gmp.MPFloat ===")
    env.out.print("")

    // 200 bits ≈ 60 significant decimal digits; we print the last 35 of those.
    let precision: USize = 10000
    let tail: USize = 35

    // ── MPFloat setup ──────────────────────────────────────────────────────────
    let mzero = MPFloat.from[F64](0.0, precision)
    let mone  = MPFloat.from[F64](1.0, precision)
    let mtwo  = MPFloat.from[F64](2.0, precision)
    var mx_n  = mtwo.sqrt()
    var my_n  = mzero
    var mpi_n = mtwo + mx_n

    // ── gmp.MPFloat setup ──────────────────────────────────────────────────────
    let gzero = gmp.MPFloat.from[F64](0.0, precision)
    let gone  = gmp.MPFloat.from[F64](1.0, precision)
    let gtwo  = gmp.MPFloat.from[F64](2.0, precision)
    var gx_n  = gtwo.sqrt()
    var gy_n  = gzero
    var gpi_n = gtwo + gx_n

    // Reference: MPFR π at same precision
    let pi_ref: String = gmp.MPFloat.pi(precision).string()

    let max_iter: USize = 15
    var iter: USize = 0
    while iter < max_iter do
      env.out.print("── Iteration " + iter.string() + " ──────────────────────────────────")

      // ── MPFloat intermediates ─────────────────────────────────────────────────
      let msqx   = mx_n.sqrt()
      let mx_n1  = (msqx + (mone / msqx)) / mtwo
      let my_n1  = ((mone + my_n) * msqx) / (mx_n + my_n)
      let mpi_n1 = (mpi_n * (mone + mx_n1) * my_n1) / (mone + my_n1)

      // ── gmp.MPFloat intermediates ─────────────────────────────────────────────
      let gsqx   = gx_n.sqrt()
      let gx_n1  = (gsqx + (gone / gsqx)) / gtwo
      let gy_n1  = ((gone + gy_n) * gsqx) / (gx_n + gy_n)
      let gpi_n1 = (gpi_n * (gone + gx_n1) * gy_n1) / (gone + gy_n1)

      // ── Collect strings (iso^ bound to iso var, then consumed to val) ─────────
      let tms_sqx: String iso = msqx.string()
      let ms_sqx: String val = consume tms_sqx
      let tgs_sqx: String iso = gsqx.string()
      let gs_sqx:String val = consume tgs_sqx
      let tms_x: String iso = mx_n1.string()
      let ms_x: String val = consume tms_x
      let tgs_x: String iso = gx_n1.string()
      let gs_x: String val = consume tgs_x
      let tms_y: String iso = my_n1.string()
      let ms_y: String val = consume tms_y
      let tgs_y: String iso = gy_n1.string()
      let gs_y: String val = consume tgs_y
      let tms_pi: String iso = mpi_n1.string()
      let ms_pi: String val = consume tms_pi
      let tgs_pi: String iso = gpi_n1.string()
      let gs_pi: String val = consume tgs_pi

      let ag_sqx = _agree(ms_sqx, gs_sqx)
      let ag_x   = _agree(ms_x,   gs_x)
      let ag_y   = _agree(ms_y,   gs_y)
      let ag_mpi_ref = _agree(ms_pi, pi_ref)
      let ag_gpi_ref = _agree(gs_pi, pi_ref)
      let ag_pi  = _agree(ms_pi,  gs_pi)
      // For pi ref line, anchor on the tighter of the two agree counts.
      let ag_ref_win = ag_mpi_ref.min(ag_gpi_ref)

      env.out.print("  sqrt(x) M: ..." + _window(ms_sqx, ag_sqx, tail) + "  agree=" + ag_sqx.string())
      env.out.print("  sqrt(x) G: ..." + _window(gs_sqx, ag_sqx, tail))
      env.out.print("  x_n+1   M: ..." + _window(ms_x,   ag_x,   tail) + "  agree=" + ag_x.string())
      env.out.print("  x_n+1   G: ..." + _window(gs_x,   ag_x,   tail))
      env.out.print("  y_n+1   M: ..." + _window(ms_y,   ag_y,   tail) + "  agree=" + ag_y.string())
      env.out.print("  y_n+1   G: ..." + _window(gs_y,   ag_y,   tail))
      env.out.print("  pi_n+1  M: ..." + _window(ms_pi,  ag_ref_win, tail) + "  agree=" + ag_mpi_ref.string() + " vs ref  |  " + ag_pi.string() + " vs GMP")
      env.out.print("  pi_n+1  G: ..." + _window(gs_pi,  ag_ref_win, tail) + "  agree=" + ag_gpi_ref.string() + " vs ref")
      env.out.print("  pi ref:    ..." + _window(pi_ref,  ag_ref_win, tail))
      env.out.print("")

      mx_n = mx_n1
      my_n = my_n1
      gx_n = gx_n1
      gy_n = gy_n1

      // Stop once pi stops changing (same string as previous iteration)
      let old_mpi_s: String iso = mpi_n.string()
      let old_ms: String val = consume old_mpi_s
      let old_gpi_s: String iso = gpi_n.string()
      let old_gs: String val = consume old_gpi_s
      mpi_n = mpi_n1
      gpi_n = gpi_n1
      if (old_ms == ms_pi) and (old_gs == gs_pi) then
        break
      end

      iter = iter + 1
    end
    env.out.print("=== end comparison ===")
    env.out.print("")


  fun e(env: Env) =>
    """
    Sample calculation of e, base of Neper Logarithm using the series:
    e = 1 + 1/1! + 1/2! + ... + 1/n!
    """
    // Sample calculation of e, base of Neper Logarithm
    //
    // e = 1 + 1/1! + 1/2! + ... + 1/1000!
    //
    let one = MPFloat.from[F64](1.0, 1000)
    var res = one
    var prev_fact = one
    var old_res = res

    env.out.print("Using MPFloat")
    env.out.print("-------------")

    for i in Range(1, 1001) do
      let fact = prev_fact * MPFloat.from[F64](F64.from[USize](i), 1000)
      res = res + (one / fact)
      prev_fact = fact

      // If value does not change, the sum has converged and there's no gain
      // to continue calculation
      if res == old_res then
        env.out.print("Stopping after " + i.string() + " iterations...")
        break
      end
      old_res = res
    end
    env.out.print("e = " + res.string())
    env.out.print("")


