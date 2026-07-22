// Concrete BigReal[T] test wrappers for the GMP/MPFR MPFloat backend.
// The BigRealFactory interface and all _Suite* primitives are in
// _tests_bigreal_suites.pony (symlinked from tests_bignum/).

use "../bignum/gmp"
use "../pony_testx"


// ── GMP MPFloat factory ──────────────────────────────────────────────────────

primitive GMPFloatFactory is BigRealFactory[MPFloat]
  fun _p(): USize => 192

  fun from_f64(v: F64): MPFloat => MPFloat.from_f64(v, _p())
  fun nan(): MPFloat            => MPFloat.nan_val(_p())
  fun pos_inf(): MPFloat        => MPFloat.inf_val(true, _p())
  fun neg_inf(): MPFloat        => MPFloat.inf_val(false, _p())
  fun pi(): MPFloat             => MPFloat.pi(_p())
  fun zero(): MPFloat           => MPFloat.from_f64(0.0, _p())

  fun ae(h: TestHelper, got: MPFloat, expected: MPFloat, msg: String) =>
    let tol = F64.epsilon().sqrt()
    h.assert_true(got.almost_eq(expected, tol, tol), msg)


// ── Concrete UnitTest wrappers: GMP MPFloat ──────────────────────────────────

class iso _TestBigRealGMPFloatPredicates is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/predicates"
  fun apply(h: TestHelper) => _SuitePredicates.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatComparisons is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/comparisons"
  fun apply(h: TestHelper) => _SuiteComparisons.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatArithmetic is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/arithmetic"
  fun apply(h: TestHelper) => _SuiteArithmetic.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatRoots is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/roots"
  fun apply(h: TestHelper) => _SuiteRoots.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatMinMax is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/min_max"
  fun apply(h: TestHelper) => _SuiteMinMax.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatRounding is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/rounding"
  fun apply(h: TestHelper) => _SuiteRounding.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatLogExp is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/log_exp"
  fun apply(h: TestHelper) => _SuiteLogExp.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatTrig is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/trig"
  fun apply(h: TestHelper) => _SuiteTrig.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatInverseTrig is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/inverse_trig"
  fun apply(h: TestHelper) => _SuiteInverseTrig.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatHyp is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/hyp"
  fun apply(h: TestHelper) => _SuiteHyp.run[MPFloat](h, GMPFloatFactory)

class iso _TestBigRealGMPFloatConversions is UnitTest
  fun name(): String => "BigReal[gmp.MPFloat]/conversions"
  fun apply(h: TestHelper) => _SuiteConversions.run[MPFloat](h, GMPFloatFactory)
