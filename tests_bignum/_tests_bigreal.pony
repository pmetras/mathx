// Concrete BigReal[T] test wrappers for the pure-Pony MPFloat backend.
// The BigRealFactory interface and all _Suite* primitives are in
// _tests_bigreal_suites.pony (same package).

use "../bignum"
use "../pony_testx"


// ── MPFloat factory ─────────────────────────────────────────────────────────

primitive MPFloatFactory is BigRealFactory[MPFloat]
  fun _p(): USize => 192

  fun from_f64(v: F64): MPFloat => MPFloat.from[F64](v, _p())
  fun nan(): MPFloat            => MPFloat.nan_val()
  fun pos_inf(): MPFloat        => MPFloat.inf_val(true)
  fun neg_inf(): MPFloat        => MPFloat.inf_val(false)
  fun pi(): MPFloat             => MPFloat.pi(_p())
  fun zero(): MPFloat           => MPFloat.create(_p())

  fun ae(h: TestHelper, got: MPFloat, expected: MPFloat, msg: String) =>
    let tol = MPFloat.epsilon(_p()).sqrt()
    h.assert_true(got.almost_eq(expected, tol, tol), msg)


// ── Concrete UnitTest wrappers: MPFloat ──────────────────────────────────────

class iso _TestBigRealMPFloatPredicates is UnitTest
  fun name(): String => "BigReal[MPFloat]/predicates"
  fun apply(h: TestHelper) => _SuitePredicates.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatComparisons is UnitTest
  fun name(): String => "BigReal[MPFloat]/comparisons"
  fun apply(h: TestHelper) => _SuiteComparisons.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatArithmetic is UnitTest
  fun name(): String => "BigReal[MPFloat]/arithmetic"
  fun apply(h: TestHelper) => _SuiteArithmetic.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatRoots is UnitTest
  fun name(): String => "BigReal[MPFloat]/roots"
  fun apply(h: TestHelper) => _SuiteRoots.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatMinMax is UnitTest
  fun name(): String => "BigReal[MPFloat]/min_max"
  fun apply(h: TestHelper) => _SuiteMinMax.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatRounding is UnitTest
  fun name(): String => "BigReal[MPFloat]/rounding"
  fun apply(h: TestHelper) => _SuiteRounding.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatLogExp is UnitTest
  fun name(): String => "BigReal[MPFloat]/log_exp"
  fun apply(h: TestHelper) => _SuiteLogExp.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatTrig is UnitTest
  fun name(): String => "BigReal[MPFloat]/trig"
  fun apply(h: TestHelper) => _SuiteTrig.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatInverseTrig is UnitTest
  fun name(): String => "BigReal[MPFloat]/inverse_trig"
  fun apply(h: TestHelper) => _SuiteInverseTrig.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatHyp is UnitTest
  fun name(): String => "BigReal[MPFloat]/hyp"
  fun apply(h: TestHelper) => _SuiteHyp.run[MPFloat](h, MPFloatFactory)

class iso _TestBigRealMPFloatConversions is UnitTest
  fun name(): String => "BigReal[MPFloat]/conversions"
  fun apply(h: TestHelper) => _SuiteConversions.run[MPFloat](h, MPFloatFactory)
