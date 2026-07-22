// BigRealFactory interface and generic test suites for BigReal[T].
// Shared between tests_bignum/ (pure-Pony backend) and tests_bignum_gmp/ (GMP backend).
// Each suite is a primitive with a `run[T]` method that accepts a factory.

use "../bignum"
use "../pony_testx"


// ── Factory interface ───────────────────────────────────────────────────────

interface val BigRealFactory[T: BigReal[T] val]
  """
  Companion factory for a BigReal[T] implementation.
  Supplies construction and approximate equality, which are outside BigReal[T].
  """
  fun from_f64(v: F64): T
  fun nan(): T
  fun pos_inf(): T
  fun neg_inf(): T
  fun pi(): T
  fun zero(): T
  fun ae(h: TestHelper, got: T, expected: T, msg: String)


// ── Test suites (generic primitives) ────────────────────────────────────────

primitive _SuitePredicates
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let zero  = f.zero()
    let one   = f.from_f64(1.0)
    let neg   = f.from_f64(-2.5)
    let nan   = f.nan()
    let pinf  = f.pos_inf()
    let ninf  = f.neg_inf()

    h.assert_true(nan.is_nan(),          "nan.is_nan()")
    h.assert_false(one.is_nan(),         "1.is_nan() = false")

    h.assert_true(pinf.is_infinite(),    "+inf.is_infinite()")
    h.assert_true(ninf.is_infinite(),    "-inf.is_infinite()")
    h.assert_false(one.is_infinite(),    "1.is_infinite() = false")

    h.assert_true(one.is_finite(),       "1.is_finite()")
    h.assert_false(pinf.is_finite(),     "+inf.is_finite() = false")
    h.assert_false(nan.is_finite(),      "nan.is_finite() = false")

    h.assert_true(zero.is_zero(),        "0.is_zero()")
    h.assert_false(one.is_zero(),        "1.is_zero() = false")

    h.assert_true(neg.is_negative(),     "-2.5.is_negative()")
    h.assert_false(one.is_negative(),    "1.is_negative() = false")
    h.assert_false(zero.is_negative(),   "0.is_negative() = false")

    h.assert_true(f.from_f64(3.0).is_integer(),  "3.is_integer()")
    h.assert_false(f.from_f64(3.5).is_integer(), "3.5.is_integer() = false")
    h.assert_false(nan.is_integer(),     "nan.is_integer() = false")
    h.assert_false(pinf.is_integer(),    "+inf.is_integer() = false")


primitive _SuiteComparisons
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let one  = f.from_f64(1.0)
    let two  = f.from_f64(2.0)
    let neg  = f.from_f64(-1.0)
    let pinf = f.pos_inf()
    let ninf = f.neg_inf()

    h.assert_true(one.lt(two),   "1 < 2")
    h.assert_true(two.gt(one),   "2 > 1")
    h.assert_true(one.le(one),   "1 <= 1")
    h.assert_true(one.ge(one),   "1 >= 1")
    h.assert_true(one.eq(one),   "1 == 1")
    h.assert_true(one.ne(two),   "1 != 2")
    h.assert_true(neg.lt(one),   "-1 < 1")
    h.assert_true(ninf.lt(pinf), "-inf < +inf")
    h.assert_true(ninf.lt(one),  "-inf < 1")
    h.assert_true(one.lt(pinf),  "1 < +inf")

    h.assert_true(one.compare(two) is Less,    "compare 1 < 2")
    h.assert_true(two.compare(one) is Greater, "compare 2 > 1")
    h.assert_true(one.compare(one) is Equal,   "compare 1 == 1")

    h.assert_true(one.sign() is Greater,      "sign(1) = Greater")
    h.assert_true(neg.sign() is Less,         "sign(-1) = Less")
    h.assert_true(f.zero().sign() is Equal,   "sign(0) = Equal")


primitive _SuiteArithmetic
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let ae = {(got: T, exp: T, msg: String) => f.ae(h, got, exp, msg)}

    let zero  = f.zero()
    let one   = f.from_f64(1.0)
    let two   = f.from_f64(2.0)
    let three = f.from_f64(3.0)
    let four  = f.from_f64(4.0)
    let half  = f.from_f64(0.5)
    let neg2  = f.from_f64(-2.0)
    let seven = f.from_f64(7.0)

    ae(one.add(one),   two,   "1 + 1 = 2")
    ae(three.sub(one), two,   "3 - 1 = 2")
    ae(two.mul(two),   four,  "2 * 2 = 4")
    ae(one.div(two),   half,  "1 / 2 = 0.5")
    ae(two.neg(),      neg2,  "neg(2) = -2")
    ae(neg2.abs(),     two,   "abs(-2) = 2")
    ae(two.inv(),      half,  "inv(2) = 0.5")

    (let q, let r) = seven.divrem(three)
    ae(q, two, "7 divrem 3 → quotient 2")
    ae(r, one, "7 divrem 3 → remainder 1")

    h.assert_true(f.nan().add(one).is_nan(), "NaN + 1 = NaN")
    h.assert_true(one.mul(f.nan()).is_nan(),  "1 * NaN = NaN")
    h.assert_true(f.pos_inf().add(one).is_infinite(), "+inf + 1 = +inf")
    h.assert_true(f.neg_inf().neg().is_infinite(),    "neg(-inf) is_infinite")
    h.assert_false(f.neg_inf().neg().is_negative(),   "neg(-inf) is positive")


primitive _SuiteRoots
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let ae = {(got: T, exp: T, msg: String) => f.ae(h, got, exp, msg)}

    let one  = f.from_f64(1.0)
    let two  = f.from_f64(2.0)
    let four = f.from_f64(4.0)
    let eight = f.from_f64(8.0)

    ae(four.sqrt(),              two,   "sqrt(4) = 2")
    ae(eight.cbrt(),             two,   "cbrt(8) = 2")
    ae(f.from_f64(16.0).rootn(4), two,  "rootn(16,4) = 2")
    ae(one.sqrt(),               one,   "sqrt(1) = 1")
    ae(one.cbrt(),               one,   "cbrt(1) = 1")

    ae(two.powi(10), f.from_f64(1024.0), "2^10 = 1024")
    ae(two.powi(-1), f.from_f64(0.5),    "2^(-1) = 0.5")
    ae(two.powi(0),  one,                "2^0 = 1")
    ae(two.pow(f.from_f64(3.0)), eight,  "2^3 = 8")

    h.assert_true(f.from_f64(-1.0).sqrt().is_nan(),  "sqrt(-1) = NaN")
    h.assert_true(f.pos_inf().sqrt().is_infinite(),   "sqrt(+inf) = +inf")


primitive _SuiteMinMax
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let ae = {(got: T, exp: T, msg: String) => f.ae(h, got, exp, msg)}

    let one = f.from_f64(1.0)
    let two = f.from_f64(2.0)
    let neg = f.from_f64(-3.0)

    ae(one.min(two), one, "min(1,2) = 1")
    ae(one.max(two), two, "max(1,2) = 2")
    ae(neg.min(one), neg, "min(-3,1) = -3")
    ae(neg.max(one), one, "max(-3,1) = 1")

    let above = one.next_above()
    let below = one.next_below()
    h.assert_true(above.gt(one), "next_above(1) > 1")
    h.assert_true(below.lt(one), "next_below(1) < 1")
    h.assert_true(above.lt(two), "next_above(1) < 2")


primitive _SuiteRounding
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let ae = {(got: T, exp: T, msg: String) => f.ae(h, got, exp, msg)}

    let two   = f.from_f64(2.0)
    let three = f.from_f64(3.0)
    let neg2  = f.from_f64(-2.0)
    let neg3  = f.from_f64(-3.0)

    ae(f.from_f64(2.7).trunc(),   two,  "trunc(2.7) = 2")
    ae(f.from_f64(-2.7).trunc(),  neg2, "trunc(-2.7) = -2")
    ae(f.from_f64(2.7).floor(),   two,  "floor(2.7) = 2")
    ae(f.from_f64(-2.3).floor(),  neg3, "floor(-2.3) = -3")
    ae(f.from_f64(2.3).ceil(),    three, "ceil(2.3) = 3")
    ae(f.from_f64(-2.7).ceil(),   neg2,  "ceil(-2.7) = -2")
    ae(f.from_f64(2.5).round(),   three, "round(2.5) = 3")
    ae(f.from_f64(2.4).round(),   two,   "round(2.4) = 2")


primitive _SuiteLogExp
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let ae = {(got: T, exp: T, msg: String) => f.ae(h, got, exp, msg)}

    let zero = f.zero()
    let one  = f.from_f64(1.0)
    let two  = f.from_f64(2.0)

    ae(one.ln(),               zero, "ln(1) = 0")
    ae(one.log(),              zero, "log(1) = 0")
    ae(two.log2(),             one,  "log2(2) = 1")
    ae(one.log2(),             zero, "log2(1) = 0")
    ae(f.from_f64(10.0).log10(), one, "log10(10) = 1")
    ae(zero.exp(),             one,  "exp(0) = 1")
    ae(zero.exp2(),            one,  "exp2(0) = 1")
    ae(one.exp2(),             two,  "exp2(1) = 2")

    // ln(e) = 1: compare via f64() since e_v is only F64-accurate
    let e_v = f.from_f64(F64.e())
    let ln_e_err = (e_v.ln().f64() - 1.0).abs()
    h.assert_true(ln_e_err < 1e-14, "ln(e) ≈ 1 (err=" + ln_e_err.string() + ")")
    let exp1_err = (one.exp().f64() - F64.e()).abs()
    h.assert_true(exp1_err < 1e-14, "exp(1) ≈ e (err=" + exp1_err.string() + ")")

    // round-trip: ln(exp(x)) = x
    ae(f.from_f64(1.5).exp().ln(), f.from_f64(1.5), "ln(exp(1.5)) = 1.5")

    h.assert_true(f.from_f64(-1.0).ln().is_nan(),   "ln(-1) = NaN")
    h.assert_true(zero.ln().is_infinite(),           "ln(0) = -inf")
    h.assert_true(f.pos_inf().exp().is_infinite(),   "exp(+inf) = +inf")


primitive _SuiteTrig
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let ae = {(got: T, exp: T, msg: String) => f.ae(h, got, exp, msg)}

    let zero = f.zero()
    let one  = f.from_f64(1.0)
    let pi   = f.pi()
    let pi2  = pi.div(f.from_f64(2.0))

    ae(zero.sin(), zero, "sin(0) = 0")
    ae(pi2.sin(),  one,  "sin(π/2) = 1")
    ae(zero.cos(), one,  "cos(0) = 1")
    ae(pi2.cos(),  zero, "cos(π/2) ≈ 0")
    ae(zero.tan(), zero, "tan(0) = 0")

    let x = f.from_f64(0.7)
    let s = x.sin()
    let c = x.cos()
    ae(s.mul(s).add(c.mul(c)), one, "sin²(0.7)+cos²(0.7) = 1")

    h.assert_true(f.nan().sin().is_nan(), "sin(NaN) = NaN")
    h.assert_true(f.nan().cos().is_nan(), "cos(NaN) = NaN")


primitive _SuiteInverseTrig
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let ae = {(got: T, exp: T, msg: String) => f.ae(h, got, exp, msg)}

    let zero    = f.zero()
    let one     = f.from_f64(1.0)
    let neg_one = f.from_f64(-1.0)
    let pi      = f.pi()
    let pi2     = pi.div(f.from_f64(2.0))
    let pi4     = pi.div(f.from_f64(4.0))

    ae(zero.asin(),    zero,      "asin(0) = 0")
    ae(one.asin(),     pi2,       "asin(1) = π/2")
    ae(neg_one.asin(), pi2.neg(), "asin(-1) = -π/2")
    ae(one.acos(),     zero,      "acos(1) = 0")
    ae(neg_one.acos(), pi,        "acos(-1) = π")
    ae(zero.acos(),    pi2,       "acos(0) = π/2")
    ae(zero.atan(),    zero,      "atan(0) = 0")
    ae(one.atan(),     pi4,       "atan(1) = π/4")
    ae(neg_one.atan(), pi4.neg(), "atan(-1) = -π/4")

    let x = f.from_f64(0.5)
    ae(x.sin().asin(), x, "asin(sin(0.5)) = 0.5")
    ae(x.tan().atan(), x, "atan(tan(0.5)) = 0.5")

    h.assert_true(f.from_f64(1.5).asin().is_nan(), "asin(1.5) = NaN")
    h.assert_true(f.from_f64(1.5).acos().is_nan(), "acos(1.5) = NaN")


primitive _SuiteHyp
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    let ae = {(got: T, exp: T, msg: String) => f.ae(h, got, exp, msg)}

    let zero = f.zero()
    let one  = f.from_f64(1.0)

    ae(zero.sinh(), zero, "sinh(0) = 0")
    ae(zero.cosh(), one,  "cosh(0) = 1")
    ae(zero.tanh(), zero, "tanh(0) = 0")

    let x  = f.from_f64(1.2)
    let ch = x.cosh()
    let sh = x.sinh()
    ae(ch.mul(ch).sub(sh.mul(sh)), one, "cosh²(1.2) - sinh²(1.2) = 1")

    ae(x.sinh().asinh(), x, "asinh(sinh(1.2)) = 1.2")
    ae(x.cosh().acosh(), x, "acosh(cosh(1.2)) = 1.2")
    ae(f.from_f64(0.5).tanh().atanh(), f.from_f64(0.5), "atanh(tanh(0.5)) = 0.5")


primitive _SuiteConversions
  fun run[T: BigReal[T] val](h: TestHelper, f: BigRealFactory[T]) =>
    // f64 round-trip
    let vals: Array[F64] = [0.0; 1.0; -1.0; 0.5; -1e-5]
    for v in vals.values() do
      let got = f.from_f64(v).f64()
      let delta = (got - v).abs()
      let scale = (v.abs() * 1e-12)
      h.assert_true(delta <= scale.max(1e-30),
        "f64 round-trip: " + v.string())
    end

    h.assert_false(f.from_f64(1.0).f32().nan(),      "f32(1.0) not NaN")
    h.assert_false(f.from_f64(1.0).f32().infinite(), "f32(1.0) not inf")
    h.assert_true(f.nan().f64().nan(),                "f64(NaN) = NaN")
    h.assert_true(f.pos_inf().f64().infinite(),       "f64(+inf) = inf")
    h.assert_true(f.from_f64(1.0).get_base() > 0,    "get_base() > 0")
    h.assert_true(f.from_f64(1.0).get_precision() >= 64, "get_precision() >= 64")
