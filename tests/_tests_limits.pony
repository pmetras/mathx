// Tests for FLimits

use "../mathx"
use "../pony_testx"


class iso _TestFLimits[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Property-based tests for FLimits[F]: verify definitional guarantees that hold
  on any IEEE 754 platform, independent of the specific values.
  """
  fun name(): String =>
    "FLimits[F]/properties"

  fun apply(h: TestHelper) =>
    let lim = FLimits[F]
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)

    // IEEE 754: binary, round-to-nearest, no guard digits needed
    h.assert_eq[ISize](2, lim.radix())
    h.assert_true(lim.round_style() >= 2)
    h.assert_eq[ISize](0, lim.guard_digits())

    // Structural sanity
    h.assert_true(lim.digit() > 0)
    h.assert_true(lim.exponent() > 0)
    h.assert_true(lim.machep() < 0)
    h.assert_true(lim.negeps() < 0)
    h.assert_true(lim.min_exponent() < 0)
    h.assert_true(lim.max_exponent() > 0)

    // epsilon: the definining property is 1 + eps != 1 but 1 + eps/2 == 1
    let eps = lim.epsilon()
    h.assert_true(eps > zero)
    h.assert_true((one + eps) != one)
    h.assert_true((one + (eps / two)) == one)
    // must match the stdlib constant
    h.assert_true(eps == F.epsilon())

    // negative_epsilon: 1 - neg_eps != 1
    let neg_eps = lim.negative_epsilon()
    h.assert_true(neg_eps > zero)
    h.assert_true((one - neg_eps) != one)

    // min_value: positive, finite, and smaller than 1
    let xmin = lim.min_value()
    h.assert_true(xmin > zero)
    h.assert_true(xmin.finite())
    h.assert_true(xmin < one)

    // max_value: finite, but doubling overflows to infinity
    let xmax = lim.max_value()
    h.assert_true(xmax.finite())
    h.assert_false((xmax + xmax).finite())
    h.assert_true(xmax > one)


class iso _TestFLimitsF32Specific is UnitTest
  """
  Test FLimits[F32] against known IEEE 754 single-precision constants.

  IEEE 754 binary32: 1 sign bit, 8 exponent bits, 23 stored mantissa bits
  (24 significant bits counting the implicit leading 1).
  """
  fun name(): String =>
    "FLimits[F32]/IEEE754-constants"

  fun apply(h: TestHelper) =>
    let lim = FLimits[F32]

    h.assert_eq[ISize](24, lim.digit())     // 23 stored + 1 implicit leading bit
    h.assert_eq[ISize](8, lim.exponent())   // 8-bit biased exponent field
    h.assert_eq[ISize](-23, lim.machep())   // eps = 2^{-23} = ULP of 1.0
    h.assert_eq[ISize](-24, lim.negeps())   // 1 - 2^{-24} is next float below 1.0


class iso _TestFLimitsF64Specific is UnitTest
  """
  Test FLimits[F64] against known IEEE 754 double-precision constants.

  IEEE 754 binary64: 1 sign bit, 11 exponent bits, 52 stored mantissa bits
  (53 significant bits counting the implicit leading 1).
  """
  fun name(): String =>
    "FLimits[F64]/IEEE754-constants"

  fun apply(h: TestHelper) =>
    let lim = FLimits[F64]

    h.assert_eq[ISize](53, lim.digit())     // 52 stored + 1 implicit leading bit
    h.assert_eq[ISize](11, lim.exponent())  // 11-bit biased exponent field
    h.assert_eq[ISize](-52, lim.machep())   // eps = 2^{-52} = ULP of 1.0
    h.assert_eq[ISize](-53, lim.negeps())   // 1 - 2^{-53} is next float below 1.0
