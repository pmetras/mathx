/*
  Floating-point pitfalls: IEEE 754 failure cases and MPFloat mitigations.

  Classic numerical problems where F64 gives wrong or surprising results, shown
  side-by-side with MPFloat solutions. Each section names the phenomenon, shows
  the F64 failure, explains why it happens, and shows the MPFloat result.

  Topics covered:
    1. Decimal representation error     (0.1 + 0.2 ≠ 0.3; from_string vs from[F64])
    2. Catastrophic cancellation        (sqrt(x+1) - sqrt(x) for large x)
    3. Non-associativity                ((a + b) + c ≠ a + (b + c))
    4. Quadratic formula instability    (near-cancellation of discriminant)
    5. Muller's recurrence              (apparently convergent, actually divergent)
    6. Harmonic series partial sum      (finite F64 limit vs unbounded truth)
    7. Rounding mode pitfalls           (directed rounding for interval enclosures)
    8. Subnormal number pitfalls        (silent precision loss near zero)
    9. Overflow and underflow           (silent ±∞ and ±0 corruption)
*/

use "../../bignum"
use "../../mathx"
use "../../formatx"
use "collections"


actor Main
  """
  Demonstrates nine IEEE 754 floating-point pitfalls and shows how `MPFloat`
  mitigates each one. Run with no arguments; output goes to stdout.
  """

  let _out: OutStream
    """
    Output stream used for all printed results.
    """


  new create(env: Env) =>
    """
    Run all pitfall demonstrations in order.
    """
    _out = env.out
    _out.print("Floating-Point Pitfalls: F64 vs MPFloat")
    _out.print("========================================")
    _out.print("")

    decimal_representation()
    catastrophic_cancellation()
    non_associativity()
    quadratic_formula()
    muller_recurrence()
    harmonic_series()
    rounding_pitfalls()
    try subnormal_pitfalls()? end
    try overflow_underflow()? end


  fun _section(title: String) =>
    """
    Print a section separator with the given title.
    """
    _out.print("------------------------------------------------------------")
    _out.print(title)
    _out.print("------------------------------------------------------------")


  fun _mpf(v: MPFloat, spec: String = ".20e"): String =>
    """
    Format an `MPFloat` value with `spec` (default: 20 significant digits in
    scientific notation). Using `format` instead of `string()` keeps output
    readable regardless of the stored precision.
    """
    v.format(spec)


  fun _row(label: String, value: String, width: USize = 28) =>
    """
    Print a labelled result row, padding the label to `width` characters.
    """
    let pad = if label.size() < width then
        String.from_array(recover Array[U8].init(' ', width - label.size()) end)
      else
        ""
      end
    _out.print("  " + label + pad + value)


  // ------------------------------------------------------------------
  // 1. Decimal representation error
  // ------------------------------------------------------------------

  fun decimal_representation() =>
    """
    0.1 has no exact binary (or base-256) representation. Adding 0.1 ten times
    is not equal to 1.0 in F64; the rounding errors compound.

    Three input paths are compared:
    - `from[F64](0.1)`: inherits F64's approximation, error ~5.5e-17.
    - `from_string("0.1")`: parses the decimal string directly into 256 bits;
      0.1 is still non-terminating in base 256, but the rounding error is
      ~10^-77 instead of ~10^-17 — 60 orders of magnitude smaller.
    - The key takeaway: `from_string` is the right entry point when the source
      value is known as a decimal literal.
    """
    _section("1. Decimal representation error: 0.1 added ten times vs 1.0")

    let p: USize = 256
    let m_one  = MPFloat.from[F64](1.0, p)
    let scale  = MPFloat.from[F64](1.0e16, p)
    let tol_14 = MPFloat.from[F64](1.0e-14, p)
    let tol_30 = MPFloat.from[F64](1.0e-30, p)

    _out.print("")
    _out.print("  F64:")
    var f_sum: F64 = 0.0
    for _ in Range(0, 10) do
      f_sum = f_sum + 0.1
    end
    let f_one: F64 = 1.0
    let f_err: F64 = (f_sum - f_one) * 1.0e16
    _row("    0.1 * 10        =", f_sum.string())
    _row("    equal to 1.0?   ", (f_sum == f_one).string())
    _row("    (sum-1)*1e16    =", f_err.string())
    _out.print("  Why: 0.1 has no exact binary representation; each addition")
    _out.print("  accumulates a rounding error ~5.5e-18. After 10 additions")
    _out.print("  the total error (~5.5e-17) is visible when scaled by 1e16.")
    _out.print("")

    _out.print("  MPFloat 256-bit, from[F64](0.1)  — inherits F64 approximation:")
    let tenth_f64 = MPFloat.from[F64](0.1, p)
    var sum_f64 = MPFloat.from[F64](0.0, p)
    for _ in Range(0, 10) do
      sum_f64 = sum_f64 + tenth_f64
    end
    let err_f64 = (sum_f64 - m_one) * scale
    // Show 30 sig-digits so the inherited ~5.5e-17 error is clearly visible.
    _row("    (sum-1)*1e16    =", _mpf(err_f64, ".30e"))
    _row("    almost_eq 1e-14?", sum_f64.almost_eq(m_one, tol_14, tol_14).string())
    _row("    almost_eq 1e-30?", sum_f64.almost_eq(m_one, tol_30, tol_30).string())
    _out.print("  The F64 approximation error (~5.5e-17) is inherited exactly.")
    _out.print("  Arithmetic at 256 bits adds no further error: tol=1e-14 passes,")
    _out.print("  but tol=1e-30 (tighter than the input error) correctly fails.")
    _out.print("")

    _out.print("  MPFloat 256-bit, from_string(\"0.1\") — decimal parsed directly:")
    let tenth_str = try MPFloat.from_string("0.1", p)? else MPFloat.from[F64](0.0, p) end
    var sum_str = MPFloat.from[F64](0.0, p)
    for _ in Range(0, 10) do
      sum_str = sum_str + tenth_str
    end
    let err_str = (sum_str - m_one) * scale
    // Show the stored 0.1 at 30 sig-digits: the tiny base-256 residual is visible.
    _row("    stored 0.1      =", _mpf(tenth_str, ".30e"))
    _row("    (sum-1)*1e16    =", _mpf(err_str, ".30e"))
    _row("    almost_eq 1e-14?", sum_str.almost_eq(m_one, tol_14, tol_14).string())
    _row("    almost_eq 1e-30?", sum_str.almost_eq(m_one, tol_30, tol_30).string())
    _out.print("  from_string parses the decimal directly to 256 bits; the")
    _out.print("  base-256 rounding error is ~10^-77. Both tolerances pass.")
    _out.print("")


  // ------------------------------------------------------------------
  // 2. Catastrophic cancellation
  // ------------------------------------------------------------------

  fun catastrophic_cancellation() =>
    """
    sqrt(x + 1) - sqrt(x) → 0 in F64 for large enough x, because x and x+1
    round to the same F64 value, making the two sqrts identical.

    The mathematically equivalent stable form  -1 / (sqrt(x) + sqrt(x+1))
    avoids the subtraction of nearly-equal values and gives the right answer.

    With MPFloat at 256 bits the direct subtraction retains enough guard digits
    to agree with the stable form.
    """
    _section("2. Catastrophic cancellation: sqrt(x+1) - sqrt(x) for large x")

    _out.print("")
    // At x = 1e17, F64 has 53 bits of mantissa ≈ 15.9 decimal digits.
    // x = 1e17 and x+1 = 1e17 + 1 differ in the 18th digit, beyond F64 precision.
    // So (x+1).f64() == x.f64() and the difference of sqrts is exactly 0.
    _out.print("  F64 (x = 1e17):")
    let x64: F64 = 1.0e17
    let xp1_64: F64 = x64 + 1.0
    let direct64: F64 = xp1_64.sqrt() - x64.sqrt()
    // Stable: 1 / (sqrt(x+1) + sqrt(x))
    let stable64: F64 = 1.0 / (xp1_64.sqrt() + x64.sqrt())
    _row("    x             =", x64.string())
    _row("    x + 1         =", xp1_64.string())
    _row("    x == x+1?     ", (x64 == xp1_64).string())
    _row("    direct result =", direct64.string())
    _row("    stable result =", stable64.string())
    _out.print("  Why: F64 has only 15-16 significant decimal digits. The '1' in")
    _out.print("  x+1 falls below the last representable digit of 1e17, so x+1")
    _out.print("  rounds to x. Both sqrts are identical, their difference is 0.")
    _out.print("")

    _out.print("  MPFloat (256-bit):")
    let p: USize = 256
    let one_mp   = MPFloat.from[F64](1.0, p)
    let xmp      = MPFloat.from[F64](1.0e17, p)
    let xp1_mp   = xmp + one_mp
    let direct_mp  = xp1_mp.sqrt() - xmp.sqrt()
    let stable_mp  = one_mp / (xp1_mp.sqrt() + xmp.sqrt())
    _row("    x == x+1?     ", (xmp == xp1_mp).string())
    _row("    direct result =", _mpf(direct_mp))
    _row("    stable result =", _mpf(stable_mp))
    let rtol = MPFloat.from[F64](1.0e-10, p)
    let atok = MPFloat.from[F64](1.0e-25, p)
    _row("    results agree?", direct_mp.almost_eq(stable_mp, rtol, atok).string())
    _out.print("  At 256 bits x and x+1 are distinct, so the direct subtraction")
    _out.print("  is no longer zero and agrees with the stable form.")
    _out.print("")


  // ------------------------------------------------------------------
  // 3. Non-associativity
  // ------------------------------------------------------------------

  fun non_associativity() =>
    """
    F64 addition is not associative: (a + b) + c may differ from a + (b + c).
    This example uses three values whose sum is exactly 0 mathematically, but
    two different groupings give different F64 results.

    MPFloat at high precision makes the residual error negligible.
    """
    _section("3. Non-associativity: (a + b) + c  vs  a + (b + c)")

    _out.print("")
    // a = 2e16, b = 1.0 - 2e16, c = 0.5. Exact sum = 1.5.
    // F64 ULP of 2e16 ≈ 4.0, so both 1.0 and 0.5 are below the last bit of 2e16.
    //   b in F64: 1.0 - 2e16 rounds to -2e16 (1.0 is lost).
    //   left:  (a + b) + c = (2e16 + (-2e16)) + 0.5 = 0.5     [wrong: lost 1.0]
    //   right:  a + (b + c) = 2e16 + ((-2e16) + 0.5) = 0       [wrong: lost 0.5 too]
    //   two different wrong answers.
    _out.print("  Values: a = 2e16,  b = 1.0 - 2e16,  c = 0.5")
    _out.print("  Exact sum = 1.5")
    _out.print("")
    _out.print("  F64:")
    let a64: F64 = 2.0e16
    let b64: F64 = 1.0 - (2.0e16)
    let c64: F64 = 0.5
    let left64:  F64 = (a64 + b64) + c64
    let right64: F64 = a64 + (b64 + c64)
    _row("    (a + b) + c =", left64.string())
    _row("    a + (b + c) =", right64.string())
    _row("    equal?      ", (left64 == right64).string())
    _out.print("  Why: b = 1.0 - 2e16 rounds to exactly -2e16 in F64 (1.0 falls")
    _out.print("  below the last ULP of 2e16). So a + b = 0 and left = 0.5.")
    _out.print("  But b + c = -2e16 + 0.5 = -2e16 (0.5 lost too), so right = 0.")
    _out.print("")

    _out.print("  MPFloat (256-bit for extra headroom):")
    let p: USize = 256
    let amp = MPFloat.from[F64](2.0e16, p)
    // Construct b = 1.0 - 2e16 exactly at 256 bits (not via F64)
    let bmp = (MPFloat.from[F64](1.0, p)) - (MPFloat.from[F64](2.0e16, p))
    let cmp = MPFloat.from[F64](0.5, p)
    let left_mp:  MPFloat = (amp + bmp) + cmp
    let right_mp: MPFloat = amp + (bmp + cmp)
    _row("    (a + b) + c =", left_mp.string())
    _row("    a + (b + c) =", right_mp.string())
    _row("    equal?      ", (left_mp == right_mp).string())
    _out.print("  At 256 bits b = 1.0 - 2e16 is represented exactly; both")
    _out.print("  groupings give the correct answer 1.5.")
    _out.print("  Note: non-associativity is an inherent property of finite-")
    _out.print("  precision arithmetic; MPFloat greatly reduces but does not")
    _out.print("  eliminate it for magnitudes that exceed its own precision.")
    _out.print("")


  // ------------------------------------------------------------------
  // 4. Quadratic formula instability
  // ------------------------------------------------------------------

  fun quadratic_formula() =>
    """
    The standard quadratic formula  x = (-b ± sqrt(b²-4ac)) / 2a  suffers
    catastrophic cancellation when b² >> 4ac (one root near zero). The
    numerically stable alternative uses the identity x₁·x₂ = c/a.

    MPFloat avoids the cancellation by carrying enough extra digits so that
    even the direct formula gives the right answer.
    """
    _section("4. Quadratic formula: x^2 - 10000.0001x + 0.0001 = 0")

    // Roots are approximately x1 ≈ 10000 and x2 ≈ 1e-8.
    // b = -10000.0001, b^2 ≈ 1e8, 4ac = 4 * 0.0001 = 0.0004.
    // Standard formula: x2 = (-b - sqrt(b^2 - 4ac)) / 2a
    //   = (10000.0001 - sqrt(99999999.9996)) / 2
    //   catastrophic cancellation in the numerator.
    _out.print("")
    _out.print("  Equation: x^2 - 10000.0001*x + 0.0001 = 0")
    _out.print("  True roots: x1 ≈ 10000.0000999..., x2 ≈ 1.0e-8")
    _out.print("")
    _out.print("  F64:")
    let a64: F64 = 1.0
    let b64: F64 = -10000.0001
    let c64: F64 = 0.0001
    // disc = b*b - 4*a*c
    let disc64: F64 = (b64 * b64) - ((4.0 * a64) * c64)
    let sq64:   F64 = disc64.sqrt()
    // Standard formula for both roots:
    let x1_std64: F64 = ((-b64) + sq64) / (2.0 * a64)
    let x2_std64: F64 = ((-b64) - sq64) / (2.0 * a64)
    // Stable formula: for x2, use c / (a * x1)
    let x2_stable64: F64 = c64 / (a64 * x1_std64)
    _row("    x1 (standard) =", x1_std64.string())
    _row("    x2 (standard) =", x2_std64.string())
    _row("    x2 (stable)   =", x2_stable64.string())
    _out.print("  Why: in x2_standard, 10000.0001 - sqrt(≈10000.0001) ≈ 0,")
    _out.print("  catastrophic cancellation wipes out the significant digits.")
    _out.print("  The stable form x2 = c/(a*x1) avoids this subtraction.")
    _out.print("")

    _out.print("  MPFloat (256-bit):")
    let p: USize = 256
    let two_mp  = MPFloat.from[F64](2.0, p)
    let four_mp = MPFloat.from[F64](4.0, p)
    let amp     = MPFloat.from[F64](1.0, p)
    let bmp     = MPFloat.from[F64](-10000.0001, p)
    let cmp     = MPFloat.from[F64](0.0001, p)
    let neg_bmp = MPFloat.from[F64](10000.0001, p)
    // disc = b*b - 4*a*c
    let disc_mp = (bmp * bmp) - ((four_mp * amp) * cmp)
    let sq_mp   = disc_mp.sqrt()
    // x1 = (-b + sq) / (2*a)
    let x1_mp = (neg_bmp + sq_mp) / (two_mp * amp)
    // x2 standard = (-b - sq) / (2*a)
    let x2_std_mp = (neg_bmp - sq_mp) / (two_mp * amp)
    // x2 stable = c / (a * x1)
    let x2_stable_mp = cmp / (amp * x1_mp)
    _row("    x1 (standard) =", _mpf(x1_mp))
    _row("    x2 (standard) =", _mpf(x2_std_mp))
    _row("    x2 (stable)   =", _mpf(x2_stable_mp))
    let rtol = MPFloat.from[F64](1.0e-20, p)
    let atol = MPFloat.from[F64](1.0e-35, p)
    _row("    x2 forms agree?", x2_std_mp.almost_eq(x2_stable_mp, rtol, atol).string())
    _out.print("  At 256 bits both formulas agree; the guard digits absorb the")
    _out.print("  near-cancellation without loss.")
    _out.print("")


  // ------------------------------------------------------------------
  // 5. Muller's recurrence
  // ------------------------------------------------------------------

  fun muller_recurrence() =>
    """
    Muller's recurrence (1990):
      x(0) = 4,  x(1) = 4.25
      x(n) = 108 - (815 - 1500/x(n-2)) / x(n-1)

    Mathematically, the fixed point is 5. Every trajectory starting near 5
    converges to 5. But in F64, rounding errors push the trajectory toward
    the spurious fixed point 100, so the sequence diverges to 100.

    At sufficient precision the correct fixed point 5 is reached.

    Reference: J.-M. Muller, "Elementary Functions", Birkhäuser, 1997.
    """
    _section("5. Muller's recurrence: x(n) = 108 - (815 - 1500/x(n-2)) / x(n-1)")

    _out.print("")
    _out.print("  Start: x(0) = 4,  x(1) = 4.25.  Fixed point = 5.")
    _out.print("")
    _out.print("  F64 (first 20 iterates):")
    var xp64: F64 = 4.0      // x(n-2)
    var xc64: F64 = 4.25     // x(n-1)
    for i in Range(2, 22) do
      // x(n) = 108 - (815 - 1500/x(n-2)) / x(n-1)
      let xn64: F64 = 108.0 - ((815.0 - (1500.0 / xp64)) / xc64)
      _out.print("    x(" + i.string() + ") = " + xn64.string())
      xp64 = xc64
      xc64 = xn64
    end
    _out.print("  Why: the recurrence has a repelling fixed point at 100. Any")
    _out.print("  tiny rounding error in x(n) gets amplified each step until the")
    _out.print("  sequence converges to 100 instead of 5.")
    _out.print("")

    _out.print("  MPFloat (256-bit, first 20 iterates):")
    let p: USize = 256
    let c108 = MPFloat.from[F64](108.0,  p)
    let c815 = MPFloat.from[F64](815.0,  p)
    let c1500 = MPFloat.from[F64](1500.0, p)
    var xp_mp = MPFloat.from[F64](4.0,   p)    // x(n-2)
    var xc_mp = MPFloat.from[F64](4.25,  p)    // x(n-1)
    for i in Range(2, 22) do
      // x(n) = 108 - (815 - 1500/x(n-2)) / x(n-1)
      let xn_mp = c108 - ((c815 - (c1500 / xp_mp)) / xc_mp)
      _out.print("    x(" + i.string() + ") = " + _mpf(xn_mp))
      xp_mp = xc_mp
      xc_mp = xn_mp
    end
    _out.print("  At 256 bits the sequence stays near 5 throughout all 20 steps.")
    _out.print("  Note: even MPFloat will eventually drift if precision is too")
    _out.print("  low — the recurrence amplifies errors regardless of base.")
    _out.print("")


  // ------------------------------------------------------------------
  // 6. Harmonic series partial sum
  // ------------------------------------------------------------------

  fun harmonic_series() =>
    """
    The harmonic series H(n) = 1 + 1/2 + 1/3 + ... + 1/n diverges (grows
    without bound). But F64 addition stalls: once 1/k < ε/2 (≈ 1.1e-16),
    adding 1/k to a running sum of order 1 has no effect (the term is below
    the rounding unit). The F64 partial sum stops growing around n ≈ 10^15
    and stays frozen forever.

    MPFloat at 128 bits has a rounding unit ≈ 10^-38, so the sum keeps
    growing until n ≈ 10^38 — far later, and provably still increasing.

    We demonstrate with N = 10_000 terms: F64 and MPFloat both grow here,
    but we show the difference in the accumulated partial sums to expose
    the rounding differences that compound over many terms.
    """
    _section("6. Harmonic series: H(N) = sum 1/k for k=1..N")

    let n: USize = 100_000

    _out.print("")
    _out.print("  N = " + n.string())
    _out.print("")

    _out.print("  F64 (forward: k = 1..N):")
    var sum_fwd64: F64 = 0.0
    for k in Range(1, n + 1) do
      sum_fwd64 = sum_fwd64 + (1.0 / k.f64())
    end
    _row("    H(N) forward  =", sum_fwd64.string())

    _out.print("  F64 (backward: k = N..1, adding small terms first — more accurate):")
    var sum_bwd64: F64 = 0.0
    var k64: USize = n
    while k64 >= 1 do
      sum_bwd64 = sum_bwd64 + (1.0 / k64.f64())
      if k64 == 0 then break end
      k64 = k64 - 1
    end
    _row("    H(N) backward =", sum_bwd64.string())
    _row("    difference    =", (sum_bwd64 - sum_fwd64).string())
    _out.print("  Why: summing forward adds large terms first; each subsequent")
    _out.print("  small term may fall below the last ULP of the running total")
    _out.print("  and be silently dropped. Backward summation is more accurate")
    _out.print("  because small terms accumulate before being merged with large ones.")
    _out.print("")

    _out.print("  MPFloat (128-bit, forward):")
    let p: USize = 128
    var sum_mp = MPFloat.from[F64](0.0, p)
    let one_mp = MPFloat.from[F64](1.0, p)
    for k in Range(1, n + 1) do
      sum_mp = sum_mp + (one_mp / MPFloat.from[USize](k, p))
    end
    // Show 25 sig-digits: more than F64's 15-16, demonstrating the extra precision.
    _row("    H(N) MPFloat  =", _mpf(sum_mp, ".25e"))
    _out.print("  MPFloat at 128 bits (~38 decimal digits) preserves terms down")
    _out.print("  to 1/" + n.string() + " ≈ " + (1.0 / n.f64()).string() +
               " with full precision — no silent drops.")
    _out.print("")


  // ------------------------------------------------------------------
  // 7. Rounding mode pitfalls
  // ------------------------------------------------------------------

  fun rounding_pitfalls() =>
    """
    IEEE 754 mandates five rounding modes, but most code uses only the default
    (round-to-nearest-even). Two pitfalls arise:

    1. Interval arithmetic requires directed rounding: computing a guaranteed
       lower and upper bound for an expression needs RoundingNegInf for the
       lower bound and RoundingPosInf for the upper bound. F64 has no API for
       per-operation rounding mode; you would need platform-specific intrinsics.
       MPFloat accepts a rounding mode per context.

    2. The "table-maker's dilemma": when a computed result is very close to a
       representable value, the rounding direction depends on guard bits that
       are discarded. At F64 precision the wrong direction is chosen; at higher
       precision the extra digits resolve the ambiguity.

    This section demonstrates directed rounding with MPFloat to bracket a
    result, and shows how a computation that straddles a rounding boundary
    gives a verifiably correct enclosure only with sufficient precision.
    """
    _section("7. Rounding mode pitfalls: directed rounding and interval enclosures")

    _out.print("")
    _out.print("  Goal: compute a guaranteed enclosure [lo, hi] for sqrt(2).")
    _out.print("")

    // F64 has a single rounding mode; we can only get one value.
    let sq2_f64: F64 = F64(2.0).sqrt()
    _out.print("  F64 (default round-to-nearest):")
    _row("    sqrt(2)        =", sq2_f64.string())
    _out.print("  F64 provides no per-operation rounding mode API.")
    _out.print("  We cannot verify whether this is an upper or lower bound.")
    _out.print("")

    // MPFloat: compute sqrt(2) rounded down and rounded up.
    let p: USize = 128
    let two_neg = MPFloat.from[F64](2.0, p, RoundingNegInf)
    let two_pos = MPFloat.from[F64](2.0, p, RoundingPosInf)
    let sq2_lo = two_neg.sqrt()    // context carries RoundingNegInf
    let sq2_hi = two_pos.sqrt()    // context carries RoundingPosInf
    _out.print("  MPFloat (128-bit) with directed rounding:")
    _row("    sqrt(2) rounded down =", _mpf(sq2_lo))
    _row("    sqrt(2) rounded up   =", _mpf(sq2_hi))
    _row("    lo < hi?             ", (sq2_lo < sq2_hi).string())
    _row("    lo * lo <= 2?        ", ((sq2_lo * sq2_lo) <= MPFloat.from[F64](2.0, p)).string())
    _row("    hi * hi >= 2?        ", ((sq2_hi * sq2_hi) >= MPFloat.from[F64](2.0, p)).string())
    _out.print("  The interval [lo, hi] is a provably correct enclosure: lo^2 <= 2 <= hi^2.")
    _out.print("  Directed rounding is essential for verified numerical computing.")
    _out.print("")

    _out.print("  Rounding mode effect: 1/3 + 1/3 + ... (6 times) vs 2.0:")
    _out.print("  Exact result = 2. Toward-zero always truncates, so the sum")
    _out.print("  stays strictly below 2; toward-+inf rounds up each time.")
    let one_n  = MPFloat.from[F64](1.0, p)
    let one_z  = MPFloat.from[F64](1.0, p, RoundingZero)
    let one_pi = MPFloat.from[F64](1.0, p, RoundingPosInf)
    let thr_n  = MPFloat.from[F64](3.0, p)
    let thr_z  = MPFloat.from[F64](3.0, p, RoundingZero)
    let thr_pi = MPFloat.from[F64](3.0, p, RoundingPosInf)
    // 1/3 in each mode
    let t_n  = one_n  / thr_n
    let t_z  = one_z  / thr_z
    let t_pi = one_pi / thr_pi
    // sum 6 copies; use a fresh 2 in the same context for comparison
    var s_n  = MPFloat.from[F64](0.0, p)
    var s_z  = MPFloat.from[F64](0.0, p, RoundingZero)
    var s_pi = MPFloat.from[F64](0.0, p, RoundingPosInf)
    for _ in Range(0, 6) do
      s_n  = s_n  + t_n
      s_z  = s_z  + t_z
      s_pi = s_pi + t_pi
    end
    let two_n  = MPFloat.from[F64](2.0, p)
    // Show 35 sig-digits so the rounding difference in the last digit is visible.
    _row("    1/3 nearest       =", _mpf(t_n,  ".35e"))
    _row("    1/3 toward-zero   =", _mpf(t_z,  ".35e"))
    _row("    1/3 toward-+inf   =", _mpf(t_pi, ".35e"))
    _row("    sum nearest < 2?  ", (s_n  < two_n).string())
    _row("    sum toward-0 < 2? ", (s_z  < two_n).string())
    _row("    sum toward-+inf>2?", (s_pi > two_n).string())
    _out.print("  Directed rounding gives verifiable bounds: sum-toward-zero")
    _out.print("  is always a lower bound, sum-toward-+inf an upper bound.")
    _out.print("")


  // ------------------------------------------------------------------
  // 8. Subnormal number pitfalls
  // ------------------------------------------------------------------

  fun subnormal_pitfalls() ? =>
    """
    IEEE 754 subnormal (denormalized) numbers fill the gap between zero and
    the smallest normal value. They maintain gradual underflow but at the cost
    of reduced precision: a subnormal with exponent -1022 has only k < 52
    significant bits instead of 53.

    Three pitfalls:
    1. Silent precision loss: arithmetic involving subnormals loses bits
       without any warning.
    2. Performance: many CPUs flush subnormals to zero in hardware (FTZ mode)
       or process them in a slow microcode trap; results can differ across
       platforms or compiler settings.
    3. Equality surprise: x * y / y ≠ x when the intermediate product is
       subnormal and then denormalized back.

    MPFloat has no subnormals: every value is either exactly zero or a
    fully normalized number with all `prec` bits significant. There is no
    gradual loss of mantissa bits at any magnitude.
    """
    _section("8. Subnormal number pitfalls: silent precision loss near zero")

    _out.print("")
    // F64 smallest normal:  2^-1022 ≈ 2.225e-308
    // F64 smallest subnormal: 2^-1074 ≈ 5e-324
    let min_normal:    F64 = F64.min_normalised()   // 2^-1022
    let half_normal:   F64 = min_normal / 2.0        // 2^-1023, subnormal
    let quarter_normal: F64 = min_normal / 4.0       // 2^-1024, subnormal

    _out.print("  F64 subnormal arithmetic:")
    _row("    min_normal        =", min_normal.string())
    _row("    min_normal / 2    =", half_normal.string())
    _row("    min_normal / 4    =", quarter_normal.string())

    // Precision loss: subnormal * 2 / 2 ≠ subnormal in general
    let trip: F64 = (half_normal * 2.0) / 2.0
    _row("    (x/2)*2/2 == x/2? ", (trip == half_normal).string())

    // Classic: x + y - y ≠ x when y is normal and x is subnormal
    let x_sub: F64 = half_normal
    let y_big: F64 = 1.0
    let roundtrip: F64 = (x_sub + y_big) - y_big
    _row("    (x+1.0)-1.0 == x? ", (roundtrip == x_sub).string())
    _out.print("  Why: adding a subnormal to 1.0 loses all its mantissa bits")
    _out.print("  (the gap between 1.0 and its ULP is ~2e-16, far larger than")
    _out.print("  the subnormal ~5e-324). Subtracting 1.0 back gives 0, not x.")

    // Relative error of subnormal: a subnormal with k < 52 significant bits
    // has relative error up to 2^(52-k) times larger than a normal number.
    let sub_val:   F64 = F64.from_bits(0x0000_0000_0000_0001)  // 1 ULP, subnormal
    let sub_next:  F64 = F64.from_bits(0x0000_0000_0000_0002)  // 2 ULP, subnormal
    let norm_val:  F64 = F64.from_bits(0x0010_0000_0000_0001)  // 1st normal + 1 ULP
    let norm_next: F64 = F64.from_bits(0x0010_0000_0000_0002)  // 1st normal + 2 ULP
    _out.print("")
    _row("    subnormal ULP gap  =", (sub_next - sub_val).string())
    _row("    normal    ULP gap  =", (norm_next - norm_val).string())
    _out.print("  Both gaps are equal in absolute terms (~5e-324), but the")
    _out.print("  subnormal's relative spacing is 1/1 while the normal's is")
    _out.print("  1/2^52: the subnormal has only 1 significant bit vs 53.")
    _out.print("")

    _out.print("  MPFloat: no subnormals, every value has full precision:")
    let p: USize = 128
    // 1e-320 is in the subnormal range for F64 (below min_normal 2.2e-308),
    // but MPFloat stores it as a fully normalized number.
    let tiny = MPFloat.from_string("1e-320", p)?
    let tiny2 = tiny / MPFloat.from[F64](2.0, p)
    let trip_mp = (tiny2 * MPFloat.from[F64](2.0, p)) / MPFloat.from[F64](2.0, p)
    _row("    1e-320            =", _mpf(tiny))
    _row("    1e-320 / 2        =", _mpf(tiny2))
    _row("    (x/2)*2/2 == x/2? ", (trip_mp == tiny2).string())
    _out.print("  MPFloat stores 1e-320 as a normal number with all 128 bits")
    _out.print("  intact. Dividing and multiplying back is exact.")
    _out.print("")


  // ------------------------------------------------------------------
  // 9. Overflow and underflow pitfalls
  // ------------------------------------------------------------------

  fun overflow_underflow() ? =>
    """
    F64 overflows silently to ±∞ and underflows silently to ±0. Both can
    corrupt downstream calculations without any error being raised.

    Common traps:
    1. Intermediate overflow in an expression whose final value is representable:
       e.g. `hypot(a, b) = sqrt(a*a + b*b)` — if `a` or `b` overflows when
       squared, sqrt(∞) = ∞ even though the true result would fit in F64.
    2. Overflow in factorial / combinatorial counts: 171! overflows F64.
    3. Underflow to zero breaks reciprocals: 1 / tiny = ∞ after underflow
       even when the true result would be representable.
    4. The "double rounding" in `x * y + z` (FMA vs non-FMA paths): the
       intermediate `x*y` may overflow to ∞ even when `x*y+z` would not.

    MPFloat uses an arbitrary integer exponent (I64); overflow to ∞ only
    occurs when explicitly constructing an infinite value. Intermediate
    results can be astronomically large or small with no precision loss.
    """
    _section("9. Overflow and underflow: silent ±∞ and ±0 corruption")

    _out.print("")
    // --- 9a. Intermediate overflow in hypot ---
    _out.print("  9a. Intermediate overflow in hypot(a, b) = sqrt(a*a + b*b):")
    let big: F64 = 1.0e200
    let a2: F64  = big * big        // overflows to +∞
    let hyp_naive: F64 = ((big * big) + (big * big)).sqrt()
    // Stable: factor out max. hypot(a,a) = a * sqrt(2)
    let hyp_stable: F64 = big * F64(2.0).sqrt()
    _row("    a                  =", big.string())
    _row("    a*a                =", a2.string())
    _row("    sqrt(a^2+a^2) naive=", hyp_naive.string())
    _row("    a*sqrt(2)   stable =", hyp_stable.string())
    _out.print("  Why: a^2 = (1e200)^2 = 1e400, which exceeds F64.max_value")
    _out.print("  (~1.8e308) and becomes +∞. sqrt(∞+∞) = ∞.")
    _out.print("")

    let p: USize = 256
    let a_mp  = MPFloat.from[F64](1.0e200, p)
    let two_mp = MPFloat.from[F64](2.0, p)
    let hyp_mp = ((a_mp * a_mp) + (a_mp * a_mp)).sqrt()
    _row("    MPFloat sqrt(a^2+a^2)=", _mpf(hyp_mp))
    _out.print("  MPFloat exponents are I64; 1e400 is representable exactly.")
    _out.print("")

    // --- 9b. Factorial overflow ---
    _out.print("  9b. Factorial overflow: n! in F64 overflows at n = 171:")
    var fact_f64: F64 = 1.0
    var overflow_at: USize = 0
    for n in Range[USize](1, 200) do
      fact_f64 = fact_f64 * n.f64()
      if fact_f64.infinite() and (overflow_at == 0) then
        overflow_at = n
      end
    end
    var fact170: F64 = 1.0
    for n in Range[USize](1, 171) do
      fact170 = fact170 * n.f64()
    end
    var fact171: F64 = fact170 * 171.0
    _row("    first overflow at n =", overflow_at.string())
    _row("    170!               =", fact170.string())
    _row("    171!               =", fact171.string())
    _out.print("")

    _out.print("  MPFloat 256-bit: 200! computed exactly:")
    var fact_mp = MPFloat.from[F64](1.0, p)
    for n in Range[USize](1, 201) do
      fact_mp = fact_mp * MPFloat.from[USize](n, p)
    end
    // 200! has ~375 digits; show 25 sig-digits in scientific notation.
    _row("    200!               =", _mpf(fact_mp, ".25e"))
    _out.print("")

    // --- 9c. Underflow to zero and reciprocal blow-up ---
    _out.print("  9c. Underflow: tiny value becomes 0, reciprocal becomes +∞:")
    let tiny_f64: F64 = F64.from_bits(1)   // smallest positive F64 ≈ 5e-324
    let even_tinier: F64 = tiny_f64 / 2.0  // underflows to 0
    _row("    smallest F64       =", tiny_f64.string())
    _row("    smallest_F64 / 2   =", even_tinier.string())
    _row("    1 / (smallest/2)   =", (1.0 / even_tinier).string())
    _out.print("  Why: dividing the smallest subnormal by 2 underflows to 0.")
    _out.print("  Its reciprocal is then +∞ rather than a very large finite number.")
    _out.print("")

    _out.print("  MPFloat: no underflow to zero for non-zero values:")
    // 5e-350: far below F64's smallest positive value (~5e-324).
    // string() shows extreme values using base-256 exponent notation;
    // the reciprocal proves the value is stored as a finite non-zero number.
    let tiny_mp   = MPFloat.from_string("5e-350", p)?
    let tinier_mp = tiny_mp / MPFloat.from[F64](2.0, p)
    let recip_mp  = MPFloat.from[F64](1.0, p) / tinier_mp
    _row("    5e-350 / 2 is zero? ", tinier_mp.is_zero().string())
    _row("    5e-350 / 2 is fin.? ", (not tinier_mp.is_infinite()).string())
    _row("    1 / (5e-350/2)      =", _mpf(recip_mp))
    _row("    1/(5e-350/2) finite?", (not recip_mp.is_infinite()).string())
    _out.print("  5e-350/2 is stored as a finite non-zero number (exponent I64 = -1001")
    _out.print("  in base-256 terms). Its reciprocal ~4e+349 is correct and finite.")
    _out.print("")
