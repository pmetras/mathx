// Tests for MPFloat using the GMP/MPFR binding (mathx/gmp/gmpfloat.pony).
// Tests for the pure Pony implementation are in _tests_mpfloat_pony.pony.
use "../mathx/gmp"

use "../pony_testx"

use "collections"
use "random"


class iso _TestGMPFloatFromF64 is UnitTest
  """
  Tests on MPFloat from_f64.
  """

  fun name(): String =>
    "GMPFloat/from_f64"

  fun apply(h: TestHelper) =>
    let zeron = MPFloat.from_f64(-0.0, 10, RoundingAwayZ)
    let zerop = MPFloat.from_f64(0.0, 10, RoundingAwayZ)
    h.assert_true(zeron == zerop, "+0 and -0 are equal")
    h.assert_true(zeron.zero(), "-0 is zero")
    //TODO How to read sign of -0.0? Pony bug?
    //h.assert_true(zeron.sign() == Less, "-0 is negative")
    h.assert_true(zerop.zero(), "+0 is zero")
    //TODO How to generate +0.0
    //h.assert_true(zerop.sign() == Greater, "+0 is positive")

    let pi1 = MPFloat.from_f64(F64.pi())
    let pi2 = MPFloat.from_f64(F64.pi(), 112, RoundingNearest)
    h.assert_true(pi1 == pi2, "Default parameters of MPFloat.create")
    h.assert_true((pi1 - pi2).zero())
    h.assert_true((pi1 + pi2) == MPFloat.from_f64(F64.pi() + F64.pi()))
    h.assert_true(pi1.finite(), "Pi is a number")

    let nann = MPFloat.from_f64(0.0 / 0.0)
    h.assert_true(nann.nan(), "Can create NaN")

    let infp = MPFloat.from_f64(1.0 / 0.0)
    let infn = MPFloat.from_f64(-1.0 / 0.0)
    h.assert_true(infp.infinite(), "Can create +inf")
    h.assert_true(infp.sign() == Greater, "+inf is positive")
    h.assert_true(infn.infinite(), "Can create -inf")
    h.assert_true(infn.sign() == Less, "-inf is negative")
    h.assert_true((infp + infn).nan(), "+inf + -inf is NaN")
    h.assert_true((infp - infn).infinite(), "+inf - -inf is infinite")
    h.assert_true((infp + nann).nan(), "NaN is absorbant")
    h.assert_true((infp + pi1).infinite(), "inf is absorbant")


class iso _TestGMPFloatCreate is UnitTest
  """
  Tests on MPFloat.create
  """

  fun name(): String =>
    "GMPFloat/create"

  fun apply(h: TestHelper) =>
    try
      // This number has an exact binary representation and can be created with
      // from_f64 constructor.
      let a = MPFloat.from_f64(0.3125)
      let b = MPFloat.from_string("0.3125")?
      h.assert_true(a == b, "Can create from string")
      let c = MPFloat.from_string("-0.3125")?
      h.assert_true(c == -b, "Can create negative number")
    else
      h.fail("Can't create from string")
    end

    try
      // These numbers have more digits that the default precision of MPFloat (112 bits)
      let d = MPFloat.from_string("1.2345678901234567890123456789012345678901234567890e123456789")?
      let e = MPFloat.from_string("-1.2345678901234567890123456789012345678901234567890e123456789")?
      h.assert_true(d == -e, "Opposite numbers")
      h.assert_true((d + e).zero())
      h.assert_true((d * d) == MPFloat.from_string("1.524157875323883675049535156256668194500838287337570492365005335e+246913578")?)
    else
      h.fail("Opposite numbers from string")
    end

    try
      let f = MPFloat.from_string("1")?
      let g = MPFloat.from_string("-7")?
      h.assert_true((f / g) == MPFloat.from_string("-0.142857142857142857142857142857142857")?, "Infinite fraction")
    else
      h.fail("Recurring fraction")
    end

    try
      h.assert_true(MPFloat.from_string("0.3125e12")? == MPFloat.from_f64(0.3125e12), "Valid number with exponent")
      // IEEE-754 F64 has 53 bits precision. We have to limit MPFloat to that
      // precision to test that the from_string value is correct.
      h.assert_true(MPFloat.from_string("-3.8125e-74", 53)? == MPFloat.from_f64(-3.8125e-74, 53), "Valid number with exponent")
      h.assert_true(MPFloat.from_string("-11.5625e+77", 53)? == MPFloat.from_f64(-11.5625e+77, 53), "Valid number with exponent")
    else
      h.fail("Create with exponents")
    end

    try
      h.assert_no_error({() ? => MPFloat.from_string("1234567890abcdef", 112, 16)?}, "Valid number in base 16")
      h.assert_no_error({() ? => MPFloat.from_string("1234567890ABCDEF", 112, 16)?}, "Valid number in base 16")
      h.assert_true(MPFloat.from_string("1234567890ABCDEF", 112, 16)? == MPFloat.from_string("1311768467294899695")?)
    else
      h.fail("Number in base 16")
    end

    try
      // The exponent must use character @ when the base is larger than 10
      h.assert_no_error({() ? => MPFloat.from_string("+1234567890ABCDEF@-1234", 112, 16)?}, "Valid number in base 16")
      h.assert_no_error({() ? => MPFloat.from_string("+1234567890abcdef@-1234", 112, 16)?}, "Valid number in base 16")
      h.assert_true(MPFloat.from_string("+1234567890ABCDEF@-1234", 112, 16)? == MPFloat.from_string("1.713162618988664507353932983315274e-1468")?)
    else
      h.fail("Number in base 16 with negative exponent")
    end

    try
      h.assert_no_error({() ? => MPFloat.from_string("-1234567890ABCDEF@+56", 112, 16)?}, "Valid number in base 16")
      h.assert_no_error({() ? => MPFloat.from_string("-1234567890abcdef@+56", 112, 16)?}, "Valid number in base 16")
      h.assert_true(MPFloat.from_string("-1234567890ABCDEF@+56", 112, 16)? == MPFloat.from_string("-3.536520791792043407087853535066778e+85")?)
    else
      h.fail("Negative number in base 16 with exponent")
    end

    try
      h.assert_no_error({() ? => MPFloat.from_string("1234567890ABCDEF@7", 112, 16)?}, "Valid number in base 16")
      h.assert_no_error({() ? => MPFloat.from_string("1234567890abcdef@7", 112, 16)?}, "Valid number in base 16")
      h.assert_true(MPFloat.from_string("1234567890ABCDEF@7", 112, 16)? == MPFloat.from_string("352125166684727486101585920")?)
    else
      h.fail("Number in base 16 with exponent")
    end

    try
      h.assert_no_error({() ? => MPFloat.from_string("0.111001010010010010010101000", 112, 2)?}, "Valid binary float")
      h.assert_true(MPFloat.from_string("0.111001010010010010010101000", 112, 2)? == MPFloat.from_string("0.895089447498321533203125")?)
    else
      h.fail("Binary float")
    end

    try
      h.assert_no_error({() ? => MPFloat.from_string("0.111001010010010010010101001e24", 112, 2)?}, "Valid binary float")
      h.assert_true(MPFloat.from_string("0.111001010010010010010101001e24", 112, 2)? == MPFloat.from_string("15017109.125")?)
    else
      h.fail("Binary float with exponent")
    end

    try
      h.assert_no_error({() ? => MPFloat.from_string("-0.111001010010010010010101001E-98", 112, 2)?}, "Valid binary float")
      h.assert_true(MPFloat.from_string("-0.111001010010010010010101001E-98", 112, 2)? == MPFloat.from_string("-2.824404310739091998286531786854574e-30")?)
    else
      h.log("Negative binary float with exponent")
    end

    try
      h.assert_true(MPFloat.from_string("inf")?.infinite(), "Can create inf")
      h.assert_true(MPFloat.from_string("+inf")?.infinite(), "Can create +inf")
      h.assert_true(MPFloat.from_string("-inf")?.infinite(), "Can create -inf")
    else
      h.fail("Can't create inf from_string")
    end

    try
      h.assert_true(MPFloat.from_string("-0")?.zero(), "Can create -0")
      h.assert_true(MPFloat.from_string("+0")?.zero(), "Can create +0")
      h.assert_true(MPFloat.from_string("0")?.zero(), "Can create 0")
      h.assert_true(MPFloat.from_string("0.0")?.zero(), "Can create 0.0")
      h.assert_true(MPFloat.from_string("0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000")?.zero(), "Can create 0.000...")
    else
      h.fail("Can't create 0 from_string")
    end

    try
      h.assert_true(MPFloat.from_string("nan")?.nan(), "Can create NaN")
      h.assert_true(MPFloat.from_string("NaN")?.nan(), "Can create NaN")
    else
      h.fail("Can't create NaN from_string")
    end

    h.assert_error({() ? => MPFloat.from_string("abcd")?}, "Invalid decimal number")
    h.assert_no_error({() ? => MPFloat.from_string("    3.14159")?}, "Leading spaces are allowed")
    h.assert_error({() ? => MPFloat.from_string("-3.14159    ")?}, "Trailing spaces not allowed")
    h.assert_error({() ? => MPFloat.from_string("3.141_59")?}, "Digit separator _ not allowed")
    h.assert_error({() ? => MPFloat.from_string("1234567890abcdef", 112, 12)?}, "Invalid characters in base 12")


class iso _TestGMPFloatFrom is UnitTest
  """
  Tests on MPFloat.from with various types
  """

  fun name(): String =>
    "GMPFloat/from"

  fun apply(h: TestHelper) =>
    let one = try MPFloat.from_string("1")? else MPFloat.from_f64(1.0) end

    h.assert_true(MPFloat.from[I8](I8(1)) == one, "Create from I8")
    h.assert_true(MPFloat.from[I16](I16(1)) == one, "Create from I16")
    h.assert_true(MPFloat.from[I32](I32(1)) == one, "Create from I32")
    h.assert_true(MPFloat.from[I64](I64(1)) == one, "Create from I64")
    h.assert_true(MPFloat.from[I128](I128(1)) == one, "Create from I128")
    h.assert_true(MPFloat.from[ISize](ISize(1)) == one, "Create from ISize")
    h.assert_true(MPFloat.from[ILong](ILong(1)) == one, "Create from ILong")

    h.assert_true(MPFloat.from[U8](U8(1)) == one, "Create from U8")
    h.assert_true(MPFloat.from[U16](U16(1)) == one, "Create from U16")
    h.assert_true(MPFloat.from[U32](U32(1)) == one, "Create from U32")
    h.assert_true(MPFloat.from[U64](U64(1)) == one, "Create from U64")
    h.assert_true(MPFloat.from[U128](U128(1)) == one, "Create from U128")
    h.assert_true(MPFloat.from[USize](USize(1)) == one, "Create from USize")
    h.assert_true(MPFloat.from[ULong](ULong(1)) == one, "Create from ULong")

    h.assert_true(MPFloat.from[F32](F32(1)) == one, "Create from F32")
    h.assert_true(MPFloat.from[F64](F64(1)) == one, "Create from F64")


class iso _TestGMPFloatFromMPFloat is UnitTest
  """
  Tests on MPFloat.from_mpfloat: assignment with copy
  """

  fun name(): String =>
    "GMPFloat/from_mpfloat"

  fun apply(h: TestHelper) =>
    let one = try MPFloat.from_string("1")? else MPFloat.from_f64(1.0) end
    let another = MPFloat.from_mpfloat(one)

    h.assert_true(another == one, "Both copies are equal")
    h.assert_false(another is one, "Different instances")

    let big_one = MPFloat.from_mpfloat(one, 1024)
    h.assert_true(big_one.string() == one.string(), "Compare with different precisions")
    h.assert_true(big_one == one, "Equal with different precisions")

    let small_one = MPFloat.from_mpfloat(one, 10)
    h.assert_true(small_one.string() == one.string(), "Compare with different precisions")
    h.assert_true(small_one == one, "Equal with different precisions")

    let inf = MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)
    let nan = MPFloat.from_f64(0.0) / MPFloat.from_f64(0.0)
    let inf' = MPFloat.from_mpfloat(inf, 1024)
    let nan' = MPFloat.from_mpfloat(nan, 1024)
    h.assert_true(inf'.infinite(), "Set infinite")
    h.assert_true(nan'.nan(), "Set NaN")


class iso _TestGMPFloatConversionF64 is UnitTest
  """
  Tests on MPFloat conversion to F64
  """

  fun name(): String =>
    "GMPFloat/f64"

  fun apply(h: TestHelper) =>
    h.assert_true(MPFloat.from_f64(1.0).f64() == 1.0, "Conversion to F64")
    h.assert_true(MPFloat.pi(1024).f64().finite(), "Conversion to number")
    h.assert_true((MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)).f64() == (-1.0 / 0.0), "Conversion to -inf")
    h.assert_true((MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)).f64().infinite(), "Conversion to -inf")
    h.assert_true((MPFloat.from_f64(53.0) / MPFloat.from_f64(0.0)).f64() == (136.12 / 0.0), "Conversion to +inf")
    h.assert_true((MPFloat.from_f64(53.0) / MPFloat.from_f64(0.0)).f64().infinite(), "Conversion to +inf")
    h.assert_true(((MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)) - (MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0))).f64().nan(), "Conversion to NaN")


class iso _TestGMPFloatConversionF32 is UnitTest
  """
  Tests on MPFloat conversion to F32
  """

  fun name(): String =>
    "GMPFloat/f32"

  fun apply(h: TestHelper) =>
    h.assert_true(MPFloat.from_f64(1.0).f32() == 1.0, "Conversion to F32")
    h.assert_true(MPFloat.pi(1024).f32().finite(), "Conversion to number")
    h.assert_true((MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)).f32() == (-1.0 / 0.0), "Conversion to -inf")
    h.assert_true((MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)).f32().infinite(), "Conversion to -inf")
    h.assert_true((MPFloat.from_f64(53.0) / MPFloat.from_f64(0.0)).f32() == (136.12 / 0.0), "Conversion to +inf")
    h.assert_true((MPFloat.from_f64(53.0) / MPFloat.from_f64(0.0)).f32().infinite(), "Conversion to +inf")
    h.assert_true(((MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)) - (MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0))).f32().nan(), "Conversion to NaN")


class iso _TestGMPFloatConversionISize is UnitTest
  """
  Tests on MPFloat conversion to ISize
  """

  fun name(): String =>
    "GMPFloat/isize"

  fun apply(h: TestHelper) =>
    h.assert_no_error({() ? => MPFloat.from_f64(1.0).isize()? == 1}, "Conversion to ISize")
    h.assert_no_error({() ? => MPFloat.from_f64(-42.8).isize()? == -43}, "Conversion negative number to ISize")
    h.assert_no_error({() ? => MPFloat.pi(1024).isize()? == 3}, "Pi == 3!")
    h.assert_error({() ? => (MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)).isize()? == 0}, "Conversion to -inf")
    h.assert_error({() ? => (MPFloat.from_f64(53.0) / MPFloat.from_f64(0.0)).isize()? == 0}, "Conversion to +inf")
    h.assert_error({() ? => ((MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)) - (MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0))).isize()? == 0}, "Conversion to NaN")


class iso _TestGMPFloatConversionUSize is UnitTest
  """
  Tests on MPFloat conversion to USize
  """

  fun name(): String =>
    "GMPFloat/usize"

  fun apply(h: TestHelper) =>
    h.assert_no_error({() ? => MPFloat.from_f64(1.0).usize()? == 1}, "Conversion to USize")
    h.assert_error({() ? => MPFloat.from_f64(-42.8).usize()? == -43}, "Can't convert negative number to USize")
    h.assert_no_error({() ? => MPFloat.pi(1024).usize()? == 3}, "Pi == 3!")
    h.assert_error({() ? => (MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)).usize()? == 0}, "Conversion to -inf")
    h.assert_error({() ? => (MPFloat.from_f64(53.0) / MPFloat.from_f64(0.0)).usize()? == 0}, "Conversion to +inf")
    h.assert_error({() ? => ((MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)) - (MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0))).usize()? == 0}, "Conversion to NaN")


class iso _TestGMPFloatConversionILong is UnitTest
  """
  Tests on MPFloat conversion to ILong
  """

  fun name(): String =>
    "GMPFloat/ilong"

  fun apply(h: TestHelper) =>
    h.assert_no_error({() ? => MPFloat.from_f64(1.0).ilong()? == 1}, "Conversion to ILong")
    h.assert_no_error({() ? => MPFloat.from_f64(-42.8).ilong()? == -43}, "Conversion negative number to ILong")
    h.assert_no_error({() ? => MPFloat.pi(1024).ilong()? == 3}, "Pi == 3!")
    h.assert_error({() ? => (MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)).ilong()? == 0}, "Conversion to -inf")
    h.assert_error({() ? => (MPFloat.from_f64(53.0) / MPFloat.from_f64(0.0)).ilong()? == 0}, "Conversion to +inf")
    h.assert_error({() ? => ((MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)) - (MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0))).ilong()? == 0}, "Conversion to NaN")


class iso _TestGMPFloatConversionULong is UnitTest
  """
  Tests on MPFloat conversion to ULong
  """

  fun name(): String =>
    "GMPFloat/ulong"

  fun apply(h: TestHelper) =>
    h.assert_no_error({() ? => MPFloat.from_f64(1.0).ulong()? == 1}, "Conversion to ULong")
    h.assert_error({() ? => MPFloat.from_f64(-42.8).ulong()? == -43}, "Can't convert negative number to ULong")
    h.assert_no_error({() ? => MPFloat.pi(1024).ulong()? == 3}, "Pi == 3!")
    h.assert_error({() ? => (MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)).ulong()? == 0}, "Conversion to -inf")
    h.assert_error({() ? => (MPFloat.from_f64(53.0) / MPFloat.from_f64(0.0)).ulong()? == 0}, "Conversion to +inf")
    h.assert_error({() ? => ((MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)) - (MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0))).ulong()? == 0}, "Conversion to NaN")


class iso _TestGMPFloatConversionString is UnitTest
  """
  Tests on string representation of MPFloat
  """

  fun name(): String =>
    "GMPFloat/string"

  fun apply(h: TestHelper) =>
    h.assert_true(MPFloat.from_f64(1.0).string() == "1", "Conversion to String")
    h.assert_true(MPFloat.from_f64(-42.8, 53).string() == "-42.8", "Conversion to String")
    try
      h.assert_true(MPFloat.from_string("-1234567890")?.string() == "-1234567890", "Conversion to String")
    else
      h.fail("Can't create from string")
    end

    try
      h.assert_true(MPFloat.from_string("-12345.67890")?.string() == "-12345.6789", "Conversion to String")
    else
      h.fail("Can't create from string")
    end

    try
      h.assert_true(MPFloat.from_string("-.1234567890")?.string() == "-0.123456789", "Conversion to String")
    else
      h.fail("Can't create from string")
    end

    try
      h.assert_true(MPFloat.from_string("-0.000001234567890")?.string() == "-1.23456789e-06", "Conversion to String")
    else
      h.fail("Can't create from string")
    end

    try
      h.assert_true(MPFloat.from_string("-12345.67890e1234")?.string() == "-1.23456789e+1238", "Conversion to String")
    else
      h.fail("Can't create from string")
    end

    h.assert_true(MPFloat.pi(120).string() == "3.141592653589793238462643383279502885", "Pi")
    h.assert_true((MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)).string() == "-inf", "Conversion to -inf")
    h.assert_true((MPFloat.from_f64(53.0) / MPFloat.from_f64(0.0)).string() == "inf", "Conversion to +inf")
    h.assert_true(((MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)) - (MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0))).string() == "nan", "Conversion to NaN")


class iso _TestGMPFloatConversionExactString is UnitTest
  """
  Tests on exact string representation of MPFloat
  """

  fun name(): String =>
    "GMPFloat/exact_string"

  fun apply(h: TestHelper) =>
    try
      let a = MPFloat.from_string("47752.9875836e749")?.exact_string()
      h.assert_true(a._1 == "47752987583599999999999999999999999")
      h.assert_true(a._2 == 754)
      h.assert_true(a._3)
    else
      h.fail("Can't create from string")
    end

    let p = MPFloat.from_f64(F64.pi(), 20).exact_string()
    h.assert_true(p._1 == "31415939")
    h.assert_true(p._2 == 1)
    h.assert_true(p._3)

    try
      let a = (MPFloat.from_string("1")? / MPFloat.from_string("0")?).exact_string()
      h.assert_true(a._1 == "inf")
    else
      h.fail("Can't create from string")
    end

    try
      let a = ((MPFloat.from_string("-1")? / MPFloat.from_string("0")?) + (MPFloat.from_string("1")? / MPFloat.from_string("0")?)).exact_string()
      h.assert_true(a._1 == "nan")
    else
      h.fail("Can't create from string")
    end

    try
      let a = MPFloat.from_string("BadBeef@12", 112, 16)?.exact_string(16)
      h.assert_true(a._1 == "badbeef0000000000000000000000")
      h.assert_true(a._2 == 19)
      h.assert_true(a._3)
    else
      h.log("Can't create hexadecimal number")
    end

    try
      let a = MPFloat.from_string("1101011.111001001100@-23", 112, 2)?.exact_string(2)
      h.assert_true(a._1 == "1101011111001001100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000")
      h.assert_true(a._2 == -16)
      h.assert_true(a._3)
    else
      h.log("Can't create hexadecimal number")
    end


class iso _TestGMPFloatConstants is UnitTest
  """
  Tests on calculated constants
  """

  fun name(): String =>
    "GMPFloat/constants"

  fun apply(h: TestHelper) =>
    try
      h.assert_true(MPFloat.pi(1000) == MPFloat.from_string("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737", 1000)?)
    else
      h.fail("Can't calculate Pi")
    end


class iso _TestGMPFloatComparisons is UnitTest
  """
  Tests on comparisons and boolean functions
  """

  fun name(): String =>
    "GMPFloat/comparisons"

  fun apply(h: TestHelper) =>
    h.assert_true(MPFloat.from_f64(0.0).zero(), "Test 0")
    h.assert_true(MPFloat.from_f64(1.0).finite(), "Test finite")
    h.assert_true(MPFloat.from_f64(2.0).number(), "Test number != inf, NaN or 0")
    h.assert_true((MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)).infinite(), "Test infinite")
    h.assert_true((MPFloat.from_f64(0.0) / MPFloat.from_f64(0.0)).nan(), "Test NaN")

    try
      h.assert_true(MPFloat.from_string("-0")? == MPFloat.from_string("+0")?, "Equality of 0")
    else
      h.log("Can't create from string")
    end

    h.assert_true(MPFloat.from_f64(-1234.0) == -MPFloat.from_f64(1234.0), "Signed equality")

    try
      h.assert_true((MPFloat.from_string("9")? / MPFloat.from_string("0")?) == (MPFloat.from_string("6")? / MPFloat.from_string("0")?), "Equality of inf")
    else
      h.log("Can't create from string")
    end

    // Let's try full range comparisons...
    try
      let a = MPFloat.from_string("-123.456e789", 1024)?
      let b = MPFloat.from_string("-830.252e-323", 1024)?
      let c = MPFloat.from_string("-0", 1024)?
      let d = MPFloat.from_string("+0", 1024)?
      let e = MPFloat.from_string("999.999e999", 1024)?
      let infm = e / c // Sign given by 0
      let infp = e / d // Sign given by 0
      let nan = c / c

      h.assert_true(a.finite(), "Finite")
      h.assert_true(b.finite(), "Finite")
      h.assert_true(c.finite(), "Finite")
      h.assert_true(d.finite(), "Finite")
      h.assert_true(e.finite(), "Finite")
      h.assert_false(infm.finite(), "Finite")
      h.assert_false(infp.finite(), "Finite")
      h.assert_false(nan.finite(), "Finite")

      h.assert_false(a.zero(), "Zero")
      h.assert_false(b.zero(), "Zero")
      h.assert_true(c.zero(), "Zero")
      h.assert_true(d.zero(), "Zero")
      h.assert_false(e.zero(), "Zero")
      h.assert_false(infm.zero(), "Zero")
      h.assert_false(infp.zero(), "Zero")
      h.assert_false(nan.zero(), "Zero")

      h.assert_true(a.number(), "Number")
      h.assert_true(b.number(), "Number")
      h.assert_false(c.number(), "Number")
      h.assert_false(d.number(), "Number")
      h.assert_true(e.number(), "Number")
      h.assert_false(infm.number(), "Number")
      h.assert_false(infp.number(), "Number")
      h.assert_false(nan.number(), "Number")

      h.assert_false(a.infinite(), "Infinite")
      h.assert_false(b.infinite(), "Infinite")
      h.assert_false(c.infinite(), "Infinite")
      h.assert_false(d.infinite(), "Infinite")
      h.assert_false(e.infinite(), "Infinite")
      h.assert_true(infm.infinite(), "Infinite")
      h.assert_true(infp.infinite(), "Infinite")
      h.assert_false(nan.infinite(), "Infinite")

      h.assert_false(a.nan(), "NaN")
      h.assert_false(b.nan(), "NaN")
      h.assert_false(c.nan(), "NaN")
      h.assert_false(d.nan(), "NaN")
      h.assert_false(e.nan(), "NaN")
      h.assert_false(infm.nan(), "NaN")
      h.assert_false(infp.nan(), "NaN")
      h.assert_true(nan.nan(), "NaN")

      h.assert_true(a == a, "==")
      h.assert_true(a < b, "<")
      h.assert_true(a <= b, "<=")
      h.assert_true(a != b, "!=")
      h.assert_false(a > b, ">")
      h.assert_false(a >= b, ">=")
      h.assert_false(a == b, "==")
      h.assert_true(a < c, "<")
      h.assert_true(a <= c, "<=")
      h.assert_true(a != c, "!=")
      h.assert_false(a > c, ">")
      h.assert_false(a >= c, ">=")
      h.assert_false(a == c, "==")
      h.assert_true(a < d, "<")
      h.assert_true(a <= d, "<=")
      h.assert_true(a != d, "!=")
      h.assert_false(a > d, ">")
      h.assert_false(a >= d, ">=")
      h.assert_false(a == d, "==")
      h.assert_true(a < e, "<")
      h.assert_true(a <= e, "<=")
      h.assert_true(a != e, "!=")
      h.assert_false(a > e, ">")
      h.assert_false(a >= e, ">=")
      h.assert_false(a == e, "==")
      h.assert_true(b < c, "<")
      h.assert_true(b <= c, "<=")
      h.assert_true(b != c, "!=")
      h.assert_false(b > c, ">")
      h.assert_false(b >= c, ">=")
      h.assert_false(b == c, "==")
      h.assert_true(b < d, "<")
      h.assert_true(b <= d, "<=")
      h.assert_true(b != d, "!=")
      h.assert_false(b > d, ">")
      h.assert_false(b >= d, ">=")
      h.assert_false(b == d, "==")
      h.assert_true(b < e, "<")
      h.assert_true(b <= e, "<=")
      h.assert_true(b != e, "!=")
      h.assert_false(b > e, ">")
      h.assert_false(b >= e, ">=")
      h.assert_false(b == e, "==")
      h.assert_false(c < d, "<")   // not -0 < +0
      h.assert_true(c <= d, "<=")
      h.assert_false(c != d, "!=") // not -0 != +0
      h.assert_false(c > d, ">")
      h.assert_true(c >= d, ">=")  // -0 >= +0
      h.assert_true(c == d, "==")  // -0 == +0
      h.assert_true(c < e, "<")
      h.assert_true(c <= e, "<=")
      h.assert_true(c != e, "!=")
      h.assert_false(c > e, ">")
      h.assert_false(c >= e, ">=")
      h.assert_false(c == e, "==")
      h.assert_true(d < e, "<")
      h.assert_true(d <= e, "<=")
      h.assert_true(d != e, "!=")
      h.assert_false(d > e, ">")
      h.assert_false(d >= e, ">=")
      h.assert_false(d == e, "==")

      // Infinity
      h.assert_true(infm < infp, "Infinite order")
      h.assert_true(infm <= infp, "Infinite order")
      h.assert_true(infp > infm, "Infinite order")
      h.assert_true(infp >= infm, "Infinite order")
      h.assert_true(infm == infm, "Infinite equality")
      h.assert_true(infp == infp, "Infinite equality")
      h.assert_true(infm != infp, "Infinite equality")

      // NaN is a special beast...
      h.assert_false(nan == nan, "NaN equality")
      h.assert_true(nan != nan, "NaN equality")
      h.assert_false(nan < nan, "NaN equality")
      h.assert_false(nan <= nan, "NaN equality")
      h.assert_false(nan > nan, "NaN equality")
      h.assert_false(nan >= nan, "NaN equality")
      h.assert_true(nan != a, "NaN equality")
      h.assert_true(a != nan, "NaN equality")
      h.assert_true(nan != infm, "NaN equality")
      h.assert_true(infm != nan, "NaN equality")
      h.assert_true(nan != infp, "NaN equality")
      h.assert_true(infp != nan, "NaN equality")
    else
      h.fail("Can't create from string")
    end


class iso _TestGMPFloatMisc is UnitTest
  """
  Tests on miscellaneous functions
  """

  fun name(): String =>
    "GMPFloat/miscellaneous"

  fun apply(h: TestHelper) =>
    let e = MPFloat.from_f64(1.0)
    h.assert_true(e.get_precision() == 112, "Precision in bits")
    h.assert_true(e.sign() == Greater, "Positive sign")
    h.assert_true(e.get_base() == 2, "Base 2")
    try
      h.assert_true(e.get_exponent()? == 1, "Exponent in base")
      h.assert_true(e.get_exponent10()? == 0, "Exponent in base 10")
    else
      h.fail("Exponent of non-regular number")
    end

    let p = MPFloat.pi(1000)
    h.assert_true(p.get_precision() == 1000, "Precision in bits")
    h.assert_true(p.sign() == Greater, "Positive sign")
    h.assert_true(p.get_base() == 2, "Base 2")
    try
      h.assert_true(p.get_exponent()? == 2, "Exponent in base")
      h.assert_true(p.get_exponent10()? == 1, "Exponent in base 10")
    else
      h.fail("Exponent of non-regular number")
    end

    let q = p.clone()
    let r = p
    h.assert_true(p == q, "Clone")
    h.assert_false(p is q, "Clones are not the same")
    h.assert_true(p is r, "Assignments are the same")

    let infm = MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)
    h.assert_true(infm.get_precision() == 112, "Default precision")
    h.assert_true(infm.sign() == Less, "Negative sign")
    h.assert_true(infm.get_base() == 2, "Base 2")
    h.assert_error({() ? => infm.get_exponent()? == 0}, "Exponent in base of -inf")
    h.assert_error({() ? => infm.get_exponent10()? == 0}, "Exponent in base 10 of -inf")

    let nan = MPFloat.from_f64(0.0) / MPFloat.from_f64(0.0)
    h.assert_true(nan.get_precision() == 112, "Default precision")
    h.assert_true(nan.sign() == Equal, "NaN sign is Equal")
    h.assert_true(nan.get_base() == 2, "Base 2")
    h.assert_error({() ? => nan.get_exponent()? == 0}, "Exponent in base of NaN")
    h.assert_error({() ? => nan.get_exponent10()? == 0}, "Exponent in base 10 of NaN")

    let one = try MPFloat.from_string("1")? else MPFloat.from_f64(1.0) end
    let zero = try MPFloat.from_string("0")? else MPFloat.from_f64(0.0) end
    let one_epsilon = one.next_toward(zero)
    let one_epsilon' = one.next_below()
    let one' = one_epsilon.next_above()
    let one_epsilon'' = one.next_toward(infm)
    let nan' = one.next_toward(nan)
    let nan'' = nan.next_toward(infm)
    let infm' = infm.next_above()
    let zero' = zero.next_below()
    let zero'' = zero.next_above()

    h.assert_true(one_epsilon < one, "Next toward")
    h.assert_true(one_epsilon' < one, "Next below")
    h.assert_true(one_epsilon == one_epsilon', "Next below")
    h.assert_true(one == one', "Next above")
    h.assert_true(one.min(one') == one, "Min")
    h.assert_true(one'.min(one) == one, "Min")
    h.assert_true(one.max(one') == one, "Min")
    h.assert_true(one'.max(one) == one, "Min")
    h.assert_true(one.max(one_epsilon) == one, "Max")
    h.assert_true(one_epsilon.max(one) == one, "Max")
    h.assert_true(one_epsilon'' == one_epsilon, "Next toward -inf")
    h.assert_true(nan'.nan(), "NaN doesn't move")
    h.assert_true(nan''.nan(), "NaN doesn't move")
    h.assert_false(infm'.infinite(), "Move from -inf")
    h.assert_true(infm < infm', "Move from -inf")
    h.assert_true(zero.sign() == Equal, "True 0")
    h.assert_true(zero' < zero, "Below 0")
    h.assert_true(zero'.sign() == Less, "Below 0")
    h.assert_true(zero'' > zero, "Above 0")
    h.assert_true(zero''.sign() == Greater, "Above 0")
    h.assert_true(one.next_toward(one) == one, "Don't move")


class iso _TestGMPFloatArithmetic is UnitTest
  """
  Tests on arithmetic operations
  """

  fun name(): String =>
    "GMPFloat/arithmetic"

  fun apply(h: TestHelper) =>
    let zero = try MPFloat.from_string()? else MPFloat.from_f64(0.0) end
    let nan = MPFloat.from_f64(0.0) / MPFloat.from_f64(0.0)
    let infm = MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)
    let infp = MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)
    let x = try MPFloat.from_string("13.78e22")? else MPFloat.from_f64(13.78e22) end
    let y = try MPFloat.from_string("864.993e21")? else MPFloat.from_f64(864.993e21) end

    h.assert_true((x + zero) == x, "Zero addition")
    h.assert_true((x - zero) == x, "Zero substraction")
    h.assert_true((x * zero) == zero, "Zero multiplication")
    h.assert_true((zero / x) == zero, "Zero division")

    try
      h.assert_true((x + y) == MPFloat.from_string("1.002793e24")?, "Addition")
      h.assert_true((x - y) == MPFloat.from_string("-7.27193e23")?, "Substraction")
      h.assert_true((x * y) == MPFloat.from_string("1.191960354e+47")?, "Multiplication")
      h.assert_true((x / y) == MPFloat.from_string("0.1593076475763387680593947003039331")?, "Division")
    else
      h.fail("Can't create from string")
    end

    h.assert_true((infm + nan).nan(), "NaN absorbs for addition")
    h.assert_true((infp + nan).nan(), "NaN absorbs for addition")
    h.assert_true((zero + nan).nan(), "NaN absorbs for addition")
    h.assert_true((x + nan).nan(), "NaN absorbs for addition")
    h.assert_true((nan + infm).nan(), "NaN absorbs for addition")
    h.assert_true((nan + infp).nan(), "NaN absorbs for addition")
    h.assert_true((nan + zero).nan(), "NaN absorbs for addition")
    h.assert_true((nan + x).nan(), "NaN absorbs for addition")
    h.assert_true((infm - nan).nan(), "NaN absorbs for substraction")
    h.assert_true((infp - nan).nan(), "NaN absorbs for substraction")
    h.assert_true((zero - nan).nan(), "NaN absorbs for substraction")
    h.assert_true((x - nan).nan(), "NaN absorbs for substraction")
    h.assert_true((nan - infm).nan(), "NaN absorbs for substraction")
    h.assert_true((nan - infp).nan(), "NaN absorbs for substraction")
    h.assert_true((nan - zero).nan(), "NaN absorbs for substraction")
    h.assert_true((nan - x).nan(), "NaN absorbs for substraction")
    h.assert_true((infm * nan).nan(), "NaN absorbs for multiplication")
    h.assert_true((infp * nan).nan(), "NaN absorbs for multiplication")
    h.assert_true((zero * nan).nan(), "NaN absorbs for multiplication")
    h.assert_true((x * nan).nan(), "NaN absorbs for multiplication")
    h.assert_true((nan * infm).nan(), "NaN absorbs for multiplication")
    h.assert_true((nan * infp).nan(), "NaN absorbs for multiplication")
    h.assert_true((nan * zero).nan(), "NaN absorbs for multiplication")
    h.assert_true((nan * x).nan(), "NaN absorbs for multiplication")
    h.assert_true((infm / nan).nan(), "NaN absorbs for division")
    h.assert_true((infp / nan).nan(), "NaN absorbs for division")
    h.assert_true((zero / nan).nan(), "NaN absorbs for division")
    h.assert_true((x / nan).nan(), "NaN absorbs for division")
    h.assert_true((nan / infm).nan(), "NaN absorbs for division")
    h.assert_true((nan / infp).nan(), "NaN absorbs for division")
    h.assert_true((nan / zero).nan(), "NaN absorbs for division")
    h.assert_true((nan / x).nan(), "NaN absorbs for division")

    h.assert_true((infm + x).infinite(), "Infinite is far away")
    h.assert_true((infp + x).infinite(), "Infinite is far away")
    h.assert_true((x + infm).infinite(), "Infinite is far away")
    h.assert_true((x + infp).infinite(), "Infinite is far away")
    h.assert_true((infm - x).infinite(), "Infinite is far away")
    h.assert_true((infp - x).infinite(), "Infinite is far away")
    h.assert_true((x - infm).infinite(), "Infinite is far away")
    h.assert_true((x - infp).infinite(), "Infinite is far away")
    h.assert_true((infm * x).infinite(), "Infinite is far away")
    h.assert_true((infp * x).infinite(), "Infinite is far away")
    h.assert_true((x * infm).infinite(), "Infinite is far away")
    h.assert_true((x * infp).infinite(), "Infinite is far away")
    h.assert_true((infm / x).infinite(), "Infinite is far away")
    h.assert_true((infp / x).infinite(), "Infinite is far away")
    h.assert_true((x / infm).zero(), "Infinite is huge")
    h.assert_true((x / infp).zero(), "Infinite is huge")


class iso _TestGMPFloatRoots is UnitTest
  """
  Tests on square and cubic roots
  """

  fun name(): String =>
    "GMPFloat/roots"

  fun apply(h: TestHelper) =>
    let infm = MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)
    let infp = MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)
    let nan = infm / infm

    h.assert_true(MPFloat.from_f64(64.0).sqrt() == MPFloat.from_f64(8.0))
    try
      h.assert_true(MPFloat.from_string("975461057789971041")?.sqrt() == MPFloat.from_f64(987654321.0))
    else
      h.fail("Can't create from string")
    end
    h.assert_true(MPFloat.from_f64(-2.0).sqrt().nan(), "Square root of negative number")
    h.assert_true(infm.sqrt().nan(), "Square root of -inf")
    h.assert_true(infp.sqrt().infinite(), "Square root of inf")
    h.assert_true(nan.sqrt().nan(), "Square root of nan")

    h.assert_true(MPFloat.from_f64(512.0).cbrt() == MPFloat.from_f64(8.0))
    h.assert_true(MPFloat.from_f64(-512.0).cbrt() == MPFloat.from_f64(-8.0))
    h.assert_true(MPFloat.from_f64(-512.0).cbrt() == -MPFloat.from_f64(8.0))
    try
      h.assert_true(MPFloat.from_string("963418328693495609108518161")?.cbrt() == MPFloat.from_string("987654321")?)
    else
      h.fail("Can't create from string")
    end
    h.assert_true(infm.cbrt().infinite(), "Cubic root of -inf")
    h.assert_true(infp.cbrt().infinite(), "Cubic root of inf")
    h.assert_true(nan.cbrt().nan(), "Cubic root of nan")

    h.assert_true(MPFloat.from_f64(16777216.0).rootn(3) == MPFloat.from_f64(256.0))
    h.assert_true(MPFloat.from_f64(16777216.0).rootn(2) == MPFloat.from_f64(4096.0))
    h.assert_true(MPFloat.from_f64(16777216.0).rootn(1) == MPFloat.from_f64(16777216.0))
    h.assert_true(MPFloat.from_f64(16777216.0).rootn(0).nan())
    h.assert_true(MPFloat.from_f64(-512.0).rootn(3) == MPFloat.from_f64(-8.0))
    h.assert_true(MPFloat.from_f64(-512.0).rootn(3) == -MPFloat.from_f64(8.0))
    try
      h.assert_true(MPFloat.from_string("-80397083515082019093")?.rootn(5) == -MPFloat.from_string("9573")?)
    else
      h.fail("Can't create from string")
    end
    h.assert_true(infm.rootn(11).infinite(), "Odd root of -inf")
    h.assert_true(infm.rootn(12).nan(), "Even root of -inf")
    h.assert_true(infp.rootn(13).infinite(), "Root of inf")
    h.assert_true(nan.rootn(14).nan(), "Root of nan")


class iso _TestGMPFloatAbs is UnitTest
  """
  Tests on absolute value and sign
  """

  fun name(): String =>
    "GMPFloat/abs"

  fun apply(h: TestHelper) =>
    let infm = MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)
    let infp = MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)
    let zerom = try MPFloat.from_string("-0")? else MPFloat.from_f64(0.0) end
    let zerop = try MPFloat.from_string("+0")? else MPFloat.from_f64(0.0) end
    let nan = infp / infm  // Negative NaN

    h.assert_true(infm.abs() == infp, "Absolute value of -inf")
    h.assert_true(infp.abs() == infp, "Absolute value of +inf")
    h.assert_true(zerom.abs() == zerop, "Absolute value of -0")
    h.assert_true(zerop.abs() == zerop, "Absolute value of +0")
    h.assert_true(nan.abs().nan(), "Absolute value of NaN")

    h.assert_true(MPFloat.from_f64(-84922.844772).abs() == MPFloat.from_f64(84922.844772))
    h.assert_true(MPFloat.from_f64(8448.99472e24).abs() == MPFloat.from_f64(8448.99472e24))
    h.assert_true(MPFloat.from_f64(-3348.229801e-38).abs().sign() == Greater)

    h.assert_true(infp.copysign(zerom) == infm, "Signed zero")
    h.assert_true(MPFloat.from_f64(98440.23884e19).copysign(infm).sign() == Less, "Sign from -inf")
    h.assert_true(MPFloat.from_f64(98440.23884e19).copysign(zerom).sign() == Less, "Sign from -0")
    h.assert_true(MPFloat.from_f64(98440.23884e19).copysign(MPFloat.from_f64(-3.14159)) == MPFloat.from_f64(-98440.23884e19), "Copysign")
    h.assert_true(MPFloat.from_f64(98440.23884e19).copysign(nan) == MPFloat.from_f64(-98440.23884e19), "Copysign of NaN")
    h.assert_true(nan.copysign(MPFloat.from_f64(-1.5)).nan(), "Copysign of NaN")


class iso _TestGMPFloatExponential is UnitTest
  """
  Tests on exponential functions
  """

  fun name(): String =>
    "GMPFloat/exponential"

  fun apply(h: TestHelper) =>
    let infm = MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)
    let infp = MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)
    let zerom = try MPFloat.from_string("-0")? else MPFloat.from_f64(0.0) end
    let zerop = try MPFloat.from_string("+0")? else MPFloat.from_f64(0.0) end
    let nan = infp / infm
    let one = MPFloat.from_f64(1.0)

    h.assert_true(MPFloat.e().log() == one, "log(e) == 1")
    h.assert_true((MPFloat.e() * MPFloat.e()).log() == (one + one), "log(e) == 1")
    h.assert_true(MPFloat.from_f64(2.0).log2() == one, "log2(2) == 1")
    h.assert_true(MPFloat.const_log2() == MPFloat.from_f64(2.0).log(), "log(2) constant")
    h.assert_true(MPFloat.from_f64(10.0).log10() == one, "log10(10) == 1")
    h.assert_true((MPFloat.from_f64(3749.884) * MPFloat.from_f64(83001.394)).log() == (MPFloat.from_f64(3749.884).log() + MPFloat.from_f64(83001.394).log()), "Log of products = sum of logs")
    h.assert_true((MPFloat.from_f64(3749.884) * MPFloat.from_f64(83001.394)).log2() == (MPFloat.from_f64(3749.884).log2() + MPFloat.from_f64(83001.394).log2()), "Log2 of products = sum of log2s")
    h.assert_true((MPFloat.from_f64(3749.884) * MPFloat.from_f64(83001.394)).log10() == (MPFloat.from_f64(3749.884).log10() + MPFloat.from_f64(83001.394).log10()), "Log10 of products = sum of log10s")

    h.assert_true(zerom.log() == infm, "log(-0) == -inf")
    h.assert_true(zerop.log() == infm, "log(+0) == -inf")
    h.assert_true(one.log() == zerom, "log(1) == 0")
    h.assert_true(one.log() == zerop, "log(1) == 0")
    h.assert_true((-one).log().nan(), "Negative numbers don't have log")
    h.assert_true(nan.log().nan(), "Log of NaN")
    h.assert_true(zerom.log2() == infm, "log2(-0) == -inf")
    h.assert_true(zerop.log2() == infm, "log2(+0) == -inf")
    h.assert_true(one.log2() == zerom, "log2(1) == 0")
    h.assert_true(one.log2() == zerop, "log2(1) == 0")
    h.assert_true((-one).log2().nan(), "Negative numbers don't have log")
    h.assert_true(nan.log2().nan(), "Log2 of NaN")
    h.assert_true(zerom.log10() == infm, "log10(-0) == -inf")
    h.assert_true(zerop.log10() == infm, "log10(+0) == -inf")
    h.assert_true(one.log10() == zerom, "log10(1) == 0")
    h.assert_true(one.log10() == zerop, "log10(1) == 0")
    h.assert_true((-one).log10().nan(), "Negative numbers don't have log")
    h.assert_true(nan.log10().nan(), "Log10 of NaN")

    h.assert_true((MPFloat.from_f64(-3722.48).exp() * MPFloat.from_f64(7983.88402).exp()) == (MPFloat.from_f64(-3722.48) + MPFloat.from_f64(7983.88402)).exp(), "exp(a) * exp(b) == exp(a + b)")
    h.assert_true((MPFloat.from_f64(-3722.48).exp2() * MPFloat.from_f64(7983.88402).exp2()) == (MPFloat.from_f64(-3722.48) + MPFloat.from_f64(7983.88402)).exp2(), "exp2(a) * exp2(b) == exp2(a + b)")
    h.assert_true(one.exp() == MPFloat.e(), "exp == e")
    h.assert_true(zerom.exp() == one, "exp(-0) == 1")
    h.assert_true(zerop.exp() == one, "exp(+0) == 1")
    h.assert_true(infm.exp() == zerop, "exp(-inf) == 0")
    h.assert_true(nan.exp().nan(), "Exponential of NaN")
    h.assert_true(one.exp2() == MPFloat.from_f64(2.0), "exp2(1) == 2")
    h.assert_true(zerom.exp2() == one, "exp2(-0) == 1")
    h.assert_true(zerop.exp2() == one, "exp2(+0) == 1")
    h.assert_true(infm.exp2() == zerop, "exp2(-inf) == 0")
    h.assert_true(nan.exp2().nan(), "Exponential of NaN")

    h.assert_true(MPFloat.from_f64(3.0).pow(MPFloat.from_f64(56.0)) == (MPFloat.from_f64(3.0).pow(MPFloat.from_f64(49.0)) * MPFloat.from_f64(3.0).pow(MPFloat.from_f64(7.0))), "a.pow(b + c) == a.pow(b) * a.pow(c)")
    h.assert_true(MPFloat.from_f64(3.0).powi(56) == (MPFloat.from_f64(3.0).powi(49) * MPFloat.from_f64(3.0).powi(7)), "a.powi(b + c) == a.powi(b) * a.powi(c)")
    h.assert_true(one.pow(MPFloat.from_f64(1737393.3e23)) == one, "1^x = 1")
    h.assert_true(one.powi(1737393) == one, "1^x = 1")
    h.assert_true(zerom.pow(MPFloat.from_f64(940028.4e-56)) == zerom, "0^x = 0")
    h.assert_true(zerop.pow(MPFloat.from_f64(940028.4e-56)) == zerom, "0^x = 0")
    h.assert_true(zerom.powi(940028) == zerom, "0^x = 0")
    h.assert_true(zerop.powi(940028) == zerom, "0^x = 0")
    h.assert_true(MPFloat.pi().pow(infm) == zerom, "N^-inf == 0")
    h.assert_true(MPFloat.pi().pow(infp) == infp, "N^+inf == +inf")
    h.assert_true(MPFloat.pi().pow(nan).nan(), "NaN power")
    h.assert_true(nan.pow(MPFloat.pi()).nan(), "Power of NaN")


class iso _TestGMPFloatTrigonometric is UnitTest
  """
  Tests on trigonometric functions
  """

  fun name(): String =>
    "GMPFloat/trigonometric"

  fun apply(h: TestHelper) =>
    let infm = MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)
    let infp = MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)
    let zerom = try MPFloat.from_string("-0")? else MPFloat.from_f64(0.0) end
    let zerop = try MPFloat.from_string("+0")? else MPFloat.from_f64(0.0) end
    let nan = infp / infm
    let one = MPFloat.from_f64(1.0)
    let two = MPFloat.from_f64(2.0)
    let pi = MPFloat.pi()
    let epsilon = MPFloat.from_f64(1.0e-30) // Accepted error

    h.assert_true(zerom.cos() == one, "cos(0) == 1")
    h.assert_true(zerop.cos() == one, "cos(0) == 1")
    h.assert_true(pi.cos() == -one, "cos(pi) == -1")
    h.assert_true((pi / two).cos().abs() <= epsilon, "cos(pi/2) == 0")
    h.assert_true((pi / (two * two)).cos() == (two.sqrt() / two), "cos(pi/4) == sqrt(2)/2")
    h.assert_true(nan.cos().nan(), "cos(NaN) == NaN")

    h.assert_true(zerom.sin() == zerom, "sin(-0) == -0")
    h.assert_true(zerop.sin() == zerop, "sin(+0) == +0")
    h.assert_true(zerom.sin() == zerop, "sin(-0) == +0")
    h.assert_true(pi.sin().abs() <= epsilon, "sin(pi) == 0")
    h.assert_true((pi / two).sin() == one, "sin(pi/2) == 1")
    h.assert_true((-pi / two).sin() == -one, "sin(-pi/2) == -1")
    h.assert_true((pi / (two * two)).sin() == (two.sqrt() / two), "sin(pi/4) == sqrt(2)/2")
    h.assert_true(nan.sin().nan(), "sin(NaN) == NaN")

    for i in Range[F64](0.0, 10.0, 0.1) do
      let n = MPFloat.from_f64(i)
      h.assert_true(((n.cos() * n.cos()) + (n.sin() * n.sin())) >= (one - epsilon), "cos(n)^2 + sin(n)^2 == 1")
    end

    h.assert_true(zerom.tan() == zerom, "tan(-0) == -0")
    h.assert_true(zerop.tan() == zerop, "tan(+0) == +0")
    h.assert_true(zerom.tan() == zerop, "tan(-0) == +0")
    h.assert_true((pi / two).tan() >= (one / epsilon), "tan(pi/2) == +inf")
    h.assert_true((-pi / two).tan() <= (-one / epsilon), "tan(-pi/2) == -inf")
    h.assert_true((pi / (two * two)).tan() == one, "tan(pi/4) == 1")
    h.assert_true(nan.tan().nan(), "tan(NaN) == NaN")

    for i in Range[F64](0.0, 10.0, 0.1) do
      let n = MPFloat.from_f64(i)
      h.assert_true((n.tan() - (n.sin() / n.cos())).abs() <= epsilon, "tan(N) == sin(N) / cos(N)")
    end

    for i in Range[F64](0.0, 1.0, 0.01) do
      let n = MPFloat.from_f64(i)
      h.assert_true((n.cos().acos() - n).abs() <= epsilon, "acos(cos(N)) == N")
      h.assert_true((n.acos().cos() - n).abs() <= epsilon, "cos(acos(N)) == N")
      h.assert_true((n.sin().asin() - n).abs() <= epsilon, "asin(sin(N)) == N")
      h.assert_true((n.asin().sin() - n).abs() <= epsilon, "sin(asin(N)) == N")
      h.assert_true((n.tan().atan() - n).abs() <= epsilon, "atan(tan(N)) == N")
      h.assert_true((n.atan().tan() - n).abs() <= epsilon, "tan(atan(N)) == N")
    end

    h.assert_true(nan.acos().nan(), "acos(NaN) == NaN")
    h.assert_true(nan.asin().nan(), "asin(NaN) == NaN")
    h.assert_true(nan.atan().nan(), "atan(NaN) == NaN")


class iso _TestGMPFloatHyperbolic is UnitTest
  """
  Tests on hyperbolic functions
  """

  fun name(): String =>
    "GMPFloat/hyperbolic"

  fun apply(h: TestHelper) =>
    let infm = MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)
    let infp = MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)
    let zerom = try MPFloat.from_string("-0")? else MPFloat.from_f64(0.0) end
    let zerop = try MPFloat.from_string("+0")? else MPFloat.from_f64(0.0) end
    let nan = infp / infm
    let one = MPFloat.from_f64(1.0)
    let two = MPFloat.from_f64(2.0)
    let e = MPFloat.e()
    let epsilon = MPFloat.from_f64(1.0e-30) // Accepted error

    h.assert_true(zerom.cosh() == one, "cosh(0) == 1")
    h.assert_true(zerop.cosh() == one, "cosh(0) == 1")
    h.assert_true(one.cosh() == (((e * e) + one) / (two * e)), "cosh(1) == (e^2 + 1) / 2e")
    h.assert_true(infm.cosh() == infp, "cosh(-inf) == +inf")
    h.assert_true(infp.cosh() == infp, "cosh(+inf) == +inf")
    h.assert_true(nan.cosh().nan(), "cosh(NaN) == NaN")

    h.assert_true(zerom.sinh() == zerom, "sinh(0) == 0")
    h.assert_true(zerop.sinh() == zerop, "sinh(0) == 0")
    h.assert_true(one.sinh() == (((e * e) - one) / (two * e)), "sinh(1) == (e^2 - 1) / 2e")
    h.assert_true(infm.sinh() == infm, "sinh(-inf) == -inf")
    h.assert_true(infp.sinh() == infp, "sinh(+inf) == +inf")
    h.assert_true(nan.sinh().nan(), "sinh(NaN) == NaN")

    h.assert_true(zerom.tanh() == zerom, "tanh(0) == 0")
    h.assert_true(zerop.tanh() == zerop, "tanh(0) == 0")
    h.assert_true(one.tanh() == (((e * e) - one) / ((e * e) + one)), "tanh(1) == (e^2 - 1) / (e^2 + 1)")
    h.assert_true(infm.tanh() == -one, "tanh(-inf) == -1")
    h.assert_true(infp.tanh() == one, "tanh(+inf) == 1")
    h.assert_true(nan.tanh().nan(), "tanh(NaN) == NaN")

    for i in Range[F64](0.0, 1.0, 0.01) do
      let n = MPFloat.from_f64(i)
      h.assert_true((n.tanh() - (n.sinh() / n.cosh())).abs() <= epsilon, "tanh(N) == sinh(N) / cosh(N)")
      h.assert_true((n.cosh().acosh() - n).abs() <= epsilon, "acosh(cosh(N)) == N")
      h.assert_true(((n + one).acosh().cosh() - (n + one)).abs() <= epsilon, "cosh(acosh(N + 1)) == N + 1")
h.log(((n + one).acosh().cosh() - (n + one)).abs().string())
      h.assert_true((n.sinh().asinh() - n).abs() <= epsilon, "asinh(sinh(N)) == N")
      h.assert_true((n.asinh().sinh() - n).abs() <= epsilon, "sinh(asinh(N)) == N")
      h.assert_true((n.tanh().atanh() - n).abs() <= epsilon, "atanh(tanh(N)) == N")
      h.assert_true((n.atanh().tanh() - n).abs() <= epsilon, "tanh(atanh(N)) == N")
    end


class iso _TestGMPFloatInteger is UnitTest
  """
  Tests on integer and ceiling functions
  """

  fun name(): String =>
    "GMPFloat/integer"

  fun apply(h: TestHelper) =>
    let infm = MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)
    let infp = MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)
    let zerom = try MPFloat.from_string("-0")? else MPFloat.from_f64(0.0) end
    let zerop = try MPFloat.from_string("+0")? else MPFloat.from_f64(0.0) end
    let nan = infp / infm
    let one = MPFloat.from_f64(1.0)
    let two = MPFloat.from_f64(2.0)
    let epsilon = MPFloat.from_f64(1.0e-20) // Accepted error

    h.assert_true(zerom.ceil() == zerop, "ceil(0) == 0")
    h.assert_true(zerop.ceil() == zerop, "ceil(0) == 0")
    h.assert_true((zerom + epsilon).ceil() == one, "ceil(eps) == 1")
    h.assert_true(MPFloat.from_f64(274.8984).ceil() == MPFloat.from_f64(275.0))
    h.assert_true(MPFloat.from_f64(-274.8984).ceil() == MPFloat.from_f64(-274.0))
    h.assert_true(infm.ceil() == infm, "ceil(-inf) == -inf")
    h.assert_true(infp.ceil() == infp, "ceil(+inf) == +inf")
    h.assert_true(nan.ceil().nan(), "ceil(NaN) == NaN")

    h.assert_true(zerom.floor() == zerop, "floor(0) == 0")
    h.assert_true(zerop.floor() == zerop, "floor(0) == 0")
    h.assert_true((zerom + epsilon).floor() == zerop, "floor(eps) == 0")
    h.assert_true(MPFloat.from_f64(274.8984).floor() == MPFloat.from_f64(274.0))
    h.assert_true(MPFloat.from_f64(-274.8984).floor() == MPFloat.from_f64(-275.0))
    h.assert_true(infm.floor() == infm, "floor(-inf) == -inf")
    h.assert_true(infp.floor() == infp, "floor(+inf) == +inf")
    h.assert_true(nan.floor().nan(), "floor(NaN) == NaN")

    h.assert_true(zerom.round() == zerop, "round(0) == 0")
    h.assert_true(zerop.round() == zerop, "round(0) == 0")
    h.assert_true((zerom + epsilon).round() == zerop, "round(eps) == 0")
    h.assert_true((zerom - epsilon).round() == zerop, "round(-eps) == 0")
    h.assert_true(MPFloat.from_f64(274.8984).round() == MPFloat.from_f64(275.0))
    h.assert_true(MPFloat.from_f64(-274.8984).round() == MPFloat.from_f64(-275.0))
    h.assert_true(infm.round() == infm, "round(-inf) == -inf")
    h.assert_true(infp.round() == infp, "round(+inf) == +inf")
    h.assert_true(nan.round().nan(), "round(NaN) == NaN")

    h.assert_true(zerom.trunc() == zerop, "trunc(0) == 0")
    h.assert_true(zerop.trunc() == zerop, "trunc(0) == 0")
    h.assert_true((zerom + epsilon).trunc() == zerop, "trunc(eps) == 0")
    h.assert_true((zerom - epsilon).trunc() == zerom, "trunc(-eps) == 0")
    h.assert_true(MPFloat.from_f64(274.8984).trunc() == MPFloat.from_f64(274.0))
    h.assert_true(MPFloat.from_f64(-274.8984).trunc() == MPFloat.from_f64(-274.0))
    h.assert_true(infm.trunc() == infm, "trunc(-inf) == -inf")
    h.assert_true(infp.trunc() == infp, "trunc(+inf) == +inf")
    h.assert_true(nan.trunc().nan(), "trunc(NaN) == NaN")

    h.assert_true((zerom %% two) == zerop, "0 mod 2 == 0")
    h.assert_true((zerop %% two) == zerop, "0 mod 2 == 0")
    h.assert_true((epsilon %% two) == epsilon, "esp mod 2 == eps")
    h.assert_true((-epsilon %% two) == -epsilon, "-eps mod 2 == -eps")
    let a = try MPFloat.from_string("274.8984")? else MPFloat.from_f64(274.8984) end
    let b = try MPFloat.from_string("0.8984")? else MPFloat.from_f64(0.8984) end
    h.assert_true(((a %% two) - b).abs() <= epsilon, "274.8984 mod 2 == 0.8984")
    h.assert_true(((a %% -two) - b).abs() <= epsilon, "274.8984 mod -2 == 0.8984")
    h.assert_true(((-a %% two) - -b).abs() <= epsilon, "-274.8984 mod 2 == -0.8984")
    h.assert_true(((-a %% -two) - -b).abs() <= epsilon, "-274.8984 mod 2 == -0.8984")
    h.assert_true((a %% zerop).nan(), "2 mod 0 == NaN")
    h.assert_true((infm %% two).nan(), "-inf mod 2 == NaN")
    h.assert_true((infp %% two).nan(), "+inf mod 2 == NaN")
    h.assert_true((two %% infm) == two, "2 mod -inf == 2")
    h.assert_true((two %% infp) == two, "+inf mod 2 == 2")
    h.assert_true((two %% nan).nan(), "2 mod NaN == NaN")
    h.assert_true((nan %% two).nan(), "NaN mod 2 == NaN")

    h.assert_true(zerom.divrem(two)._1 == MPFloat.from_f64(0.0), "0 divrem 2 = (0, 0)")
    h.assert_true(zerom.divrem(two)._2 == MPFloat.from_f64(0.0), "0 divrem 2 = (0, 0)")

    let x = MPFloat.from_f64(74.80038)
    for i in Range[F64](-100.0, 100.0) do
      let n = MPFloat.from_f64(i)
      h.assert_true((x + n).ceil() == (x.ceil() + n), "ceil(x + N) == ceil(x) + N")
      h.assert_true((x + n).floor() == (x.floor() + n), "floor(x + N) == floor(x) + N")
      h.assert_true(((n / two).floor() + (n / two).ceil()) == n, "n = floor(n/2) + ceil(n/2)")
      if n > zerop then
        h.assert_true((x %% n) == (x - (n * (x / n).floor())), "x mod n = x - n * floor(x / n)")
      end
    end


class iso _TestGMPFloatPrecision is UnitTest
  """
  Tests on precision change
  """

  fun name(): String =>
    "GMPFloat/precision"

  fun apply(h: TestHelper) =>
    let infm = MPFloat.from_f64(-1.0) / MPFloat.from_f64(0.0)
    let infp = MPFloat.from_f64(1.0) / MPFloat.from_f64(0.0)
    let zerom = try MPFloat.from_string("-0")? else MPFloat.from_f64(0.0) end
    let zerop = try MPFloat.from_string("+0")? else MPFloat.from_f64(0.0) end
    let nan = infp / infm
    let test = try MPFloat.from_string("1984.69857721920012888340998141948e237", 256)? else MPFloat.from_f64(1984.69e23, 256) end

    h.assert_true(infm.get_precision() == 112, "Default precision")
    h.assert_true(infp.get_precision() == 112, "Default precision")
    h.assert_true(zerom.get_precision() == 112, "Default precision")
    h.assert_true(zerom.get_precision() == 112, "Default precision")
    h.assert_true(nan.get_precision() == 112, "Default precision")
    h.assert_true(test.get_precision() == 256, "Test has larger precision")

    infm.change_precision(1024)
    infp.change_precision(1024)
    zerom.change_precision(24)
    zerop.change_precision(24)
    nan.change_precision(1)
    test.change_precision(128)

    h.assert_true(infm.get_precision() == 1024, "Increase precision")
    h.assert_true(infp.get_precision() == 1024, "Increase precision")
    h.assert_true(zerom.get_precision() == 24, "Decrease precision")
    h.assert_true(zerom.get_precision() == 24, "Decrease precision")
    h.assert_true(nan.get_precision() == 1, "No precision")
    h.assert_true(test.get_precision() == 128, "Reduce precision")

    h.assert_true(infm.min_precision() == 0, "Special value")
    h.assert_true(infp.min_precision() == 0, "Special value")
    h.assert_true(zerom.min_precision() == 0, "Special value")
    h.assert_true(zerom.min_precision() == 0, "Special value")
    h.assert_true(nan.min_precision() == 0, "Special value")
    h.assert_true(test.min_precision() == 125, "Does not fit in F64")

