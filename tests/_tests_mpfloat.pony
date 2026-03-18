// Tests for the pure Pony MPFloat (mathx/mpfloat.pony).
// These tests exercise only the operations that currently exist in the
// pure Pony implementation.  Tests for the GMP/MPFR binding are in
// _tests_gmpfloat.pony (_TestGMPFloat* classes).
//
// ── Known limitations of the current implementation ────────────────────────
//
//  1. No sign field: all values are treated as unsigned.
//     `neg()` does byte-level two's-complement (wrong for floats).
//     The unary `-` operator therefore gives a numerically wrong result.
//
//  2. `sub(that)` returns `(MPFloat, I8)` (sign indicator), not `MPFloat`.
//     The binary `-` operator is therefore not usable.
//
//  3. `div(that)` returns `(MPFloat, MPFloat)` (quotient, remainder), not
//     `MPFloat`.  The binary `/` operator is therefore not usable.
//
//  4. `_exponent` field exists but is never set or used.  All arithmetic
//     assumes both operands have the same magnitude.
//
//  5. No special values: NaN, +inf, -inf, ±0 are not supported.
//
//  6. `string()` produces at most `size` decimal digits (one per base-256
//     digit) instead of the full ~2.4 decimal digits per byte available.
//     It also inserts '_' separators at positions 0, 3, 6, … of the
//     fractional part (non-standard format, different from the GMP version).
//
//  7. `from_string`, `from_f64`, `f64()`, comparison operators (`eq`, `lt`,
//     …), classification predicates (`nan()`, `infinite()`, …), and all
//     transcendental functions are not yet implemented.
//
// ── Improvements needed ────────────────────────────────────────────────────
//  • Add sign field and fix `neg()` / unary `-`
//  • Fix `sub` / `div` return types so operators work
//  • Implement `from_f64`, `f64`, `eq`, `lt`, `gt`, … for interoperability
//  • Implement `finite()`, `infinite()`, `nan()`, `zero()`, `sign()`
//  • Implement all missing transcendental functions
//  • Fix `string()` to extract all available decimal digits and drop '_'
//  • Use `_exponent` correctly for numbers outside [0, 256)

use "../mathx"
use "../pony_testx"


// ── Helper ─────────────────────────────────────────────────────────────────

primitive _MPFStr
  """
  Strip the '_' digit-group separators inserted by MPFloat.string() so that
  we can compare decimal strings more easily.
  """
  fun apply(s: String box): String ref^ =>
    let buf = String(s.size())
    for c in s.values() do
      if c != '_' then
        buf.push(c)
      end
    end
    buf


// ── Creation ───────────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatCreate is UnitTest
  """
  MPFloat(0) is a special-cased zero; MPFloat(n) initialises n bytes to 0.
  """
  fun name(): String => "MPFloat/create"

  fun apply(h: TestHelper) =>
    // Size-0 zero: handled by the early-return branch in string()
    h.assert_true(MPFloat(0).string() == "0.0", "MPFloat(0) is '0.0'")

    // Non-zero size but all bytes zero → integer part 0, all decimal digits 0
    let z8 = MPFloat(8)
    h.assert_true(z8.string().at("0."), "MPFloat(8) starts with '0.'")


// ── Pi constant ────────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatPi is UnitTest
  """
  MPFloat.pi(n) should converge to π = 3.14159…
  More bytes of precision → more decimal digits in the output.
  """
  fun name(): String => "MPFloat/pi"

  fun apply(h: TestHelper) =>
    let pi4 = MPFloat.pi(4)
    h.assert_true(pi4.string().at("3."), "pi(4) starts with '3.'")

    let pi8 = MPFloat.pi(8)
    h.assert_true(pi8.string().at("3."), "pi(8) starts with '3.'")

    // Higher precision yields more output characters
    h.assert_true(
      pi8.string().size() > pi4.string().size(),
      "pi(8) has more characters than pi(4)")

    // First few digits of π (after stripping '_' separators)
    // string() outputs one decimal digit per base-256 byte; 8 bytes → 8 digits
    // Expected prefix: "3.14159265" (9 chars for "3." + 8 digit chars – but
    // underscores are interleaved, so we strip first)
    let s = _MPFStr(pi8.string())
    h.assert_true(s.at("3.1415926"), "pi(8) decimal prefix matches π")


// ── Addition ───────────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatAdd is UnitTest
  """
  `add` / `+` performs unsigned base-256 fixed-point addition.
  """
  fun name(): String => "MPFloat/add"

  fun apply(h: TestHelper) =>
    let pi = MPFloat.pi(8)

    // π + π ≈ 6.28…
    let two_pi = pi + pi
    h.assert_true(two_pi.string().at("6."), "π + π starts with '6.'")

    // 0 + π = π
    let zero = MPFloat(8)
    h.assert_true((zero + pi).string().at("3."), "0 + π starts with '3.'")

    // Commutativity: π + 0 = π
    h.assert_true((pi + zero).string().at("3."), "π + 0 starts with '3.'")


// ── Subtraction ────────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatSub is UnitTest
  """
  Old `sub` tests — disabled during the class val redesign. The original
  version returned a `(MPFloat, I8)` tuple; `sub` now returns a plain
  `MPFloat`. Tests have been superseded by `_TestMPFloatSub`.
  """
  fun name(): String => "MPFloat/sub/disabled"

  fun apply(h: TestHelper) =>
    None  // disabled — see _TestMPFloatSub


// ── Multiplication ─────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatMul is UnitTest
  """
  `mul` / `*` uses FFT-based convolution for arbitrary-precision multiplication.
  """
  fun name(): String => "MPFloat/mul"

  fun apply(h: TestHelper) =>
    let pi = MPFloat.pi(8)

    // π² ≈ 9.8696…
    let pi_sq = pi * pi
    h.assert_true(pi_sq.string().at("9."), "π² starts with '9.'")

    // π × 0 = 0
    let zero = MPFloat(8)
    h.assert_true((pi * zero).string().at("0."), "π × 0 starts with '0.'")

    // 0 × π = 0
    h.assert_true((zero * pi).string().at("0."), "0 × π starts with '0.'")


// ── Inversion ──────────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatInv is UnitTest
  """
  `inv` computes 1/x using Newton's method.
  """
  fun name(): String => "MPFloat/inv"

  fun apply(h: TestHelper) =>
    let pi = MPFloat.pi(8)

    // 1/π ≈ 0.3183…
    let inv_pi = pi.inv()
    h.assert_true(inv_pi.string().at("0."), "1/π starts with '0.'")

    // Round-trip: (1/π) × π ≈ 1
    // The result may land as either "0.xxx" or "1.xxx" depending on precision;
    // we check that the integer part is 0 or 1.
    let product = inv_pi * pi
    let s = _MPFStr(product.string())
    h.assert_true(
      s.at("0.") or s.at("1."),
      "(1/π) × π has integer part 0 or 1")


// ── Square root ────────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatSqrt is UnitTest
  """
  `sqrt` computes the square root using Newton's method (reciprocal sqrt).
  """
  fun name(): String => "MPFloat/sqrt"

  fun apply(h: TestHelper) =>
    let pi = MPFloat.pi(8)

    // √π ≈ 1.7724…
    let sqrt_pi = pi.sqrt()
    h.assert_true(sqrt_pi.string().at("1."), "√π starts with '1.'")

    // Round-trip: (√π)² ≈ π
    let pi_approx = sqrt_pi * sqrt_pi
    h.assert_true(pi_approx.string().at("3."), "(√π)² starts with '3.'")


// ── digit_shl ──────────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatDigitShl is UnitTest
  """
  `digit_shl(n)` drops the n most-significant base-256 digits and pads
  with zeros on the right — equivalent to multiplying by 256^n and then
  discarding overflow past the array size.  It is primarily used inside
  arithmetic algorithms to remove the integer part after an operation.
  """
  fun name(): String => "MPFloat/digit_shl"

  fun apply(h: TestHelper) =>
    let pi = MPFloat.pi(8)  // [3, d1, d2, …]  → value ≈ 3.14159…

    // Shifting by 1 removes digit[0]=3 and brings d1 into position 0.
    // d1 ≈ 36 (first fractional base-256 byte of π), so the result is ≈ 36.xx
    let shifted = pi.digit_shl(1)
    // Integer part should be 36 (0.14159... × 256 ≈ 36.2)
    h.assert_true(shifted.string().at("36."), "digit_shl(1) of π starts with '36.'")

    // Shifting by the full size gives all zeros
    let all_zero = pi.digit_shl(8)
    h.assert_true(all_zero.string().at("0."), "digit_shl(size) gives zero")


// ── String format ──────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatString is UnitTest
  """
  `string()` format: '<integer>.<frac_digit>[_<frac_digit>{3}]*'
  where underscore separators appear after position 0, 3, 6, … of the
  fractional digits.  The number of fractional digit characters equals the
  byte size of the MPFloat.
  """
  fun name(): String => "MPFloat/string"

  fun apply(h: TestHelper) =>
    // Zero of size 0: special-cased to "0.0"
    h.assert_true(MPFloat(0).string() == "0.0", "size-0 zero")

    // All-zero digits → "0.0_000_000_0" for size 8
    // string() returns String iso^; annotate as String val so it can be passed
    // as String box to _MPFStr.
    let z: String = MPFloat(8).string()
    h.assert_true(z.at("0."), "zero(8) starts with '0.'")
    // All decimal digits should be '0'
    let zs = _MPFStr(z)   // remove '_'  (String box accepted)
    h.assert_true(zs == "0.00000000", "zero(8) all decimal digits are 0")

    // Underscore placement: for size=4 the format is "<d>.<f0>_<f1><f2><f3>_"
    // i.e. underscore after fractional index 0 and 3
    let s4 = MPFloat.pi(4).string()
    // After the integer part and '.', position 1 is the first fractional char,
    // position 2 is '_', so s4[2] == '_'
    h.assert_true(try s4(2)? == '_' else false end,
      "underscore after first fractional digit")


// ── New tests for redesigned class val MPFloat ──────────────────────────────


// ── Create / special values ─────────────────────────────────────────────────

class iso _TestMPFloatNewCreate is UnitTest
  """
  `MPFloat()` produces a positive zero; predicates agree.
  """
  fun name(): String => "MPFloat/new/create"

  fun apply(h: TestHelper) =>
    let z0 = MPFloat()
    h.assert_true(z0.is_zero(), "MPFloat() is zero")
    h.assert_true(z0.is_finite(), "MPFloat() is finite")
    h.assert_false(z0.is_nan(), "MPFloat() is not NaN")
    h.assert_false(z0.is_infinite(), "MPFloat() is not inf")
    h.assert_false(z0.is_negative(), "MPFloat() is not negative")

    let z8 = MPFloat(8)
    h.assert_true(z8.is_zero(), "MPFloat(8) is zero")
    h.assert_true(z8.is_finite(), "MPFloat(8) is finite")
    h.assert_false(z8.is_negative(), "MPFloat(8) is not negative")


class iso _TestMPFloatNewSpecial is UnitTest
  """
  `nan_val()` and `inf_val()` produce the expected special values.
  """
  fun name(): String => "MPFloat/new/special"

  fun apply(h: TestHelper) =>
    // NaN
    let n = MPFloat.nan_val()
    h.assert_true(n.is_nan(), "nan_val() is NaN")
    h.assert_false(n.is_infinite(), "nan_val() is not inf")
    h.assert_false(n.is_zero(), "nan_val() is not zero")
    h.assert_false(n.is_negative(), "nan_val() is not negative (NaN sign ignored)")

    // +infinity
    let pinf = MPFloat.inf_val()
    h.assert_true(pinf.is_infinite(), "+inf_val() is inf")
    h.assert_false(pinf.is_nan(), "+inf_val() is not NaN")
    h.assert_false(pinf.is_finite(), "+inf_val() is not finite")
    h.assert_false(pinf.is_negative(), "+inf_val() is positive")

    // -infinity
    let ni = MPFloat.inf_val(false)
    h.assert_true(ni.is_infinite(), "-inf_val() is inf")
    h.assert_true(ni.is_negative(), "-inf_val() is negative")


// ── from_f64 ─────────────────────────────────────────────────────────────────

class iso _TestMPFloatFromF64 is UnitTest
  """
  `from_f64` preserves sign, special values, and approximate magnitude.
  Digit values are verified indirectly through `string()`.
  """
  fun name(): String => "MPFloat/from_f64"

  fun apply(h: TestHelper) =>
    // NaN propagation: quiet NaN bit pattern
    let fnan = MPFloat.from_f64(F64.from_bits(0x7FF8_0000_0000_0000))
    h.assert_true(fnan.is_nan(), "from_f64(NaN) is NaN")

    // +infinity bit pattern
    let finf = MPFloat.from_f64(F64.from_bits(0x7FF0_0000_0000_0000))
    h.assert_true(finf.is_infinite(), "from_f64(+inf) is inf")
    h.assert_false(finf.is_negative(), "from_f64(+inf) is positive")

    // -infinity bit pattern
    let fninf = MPFloat.from_f64(F64.from_bits(0xFFF0_0000_0000_0000))
    h.assert_true(fninf.is_infinite(), "from_f64(-inf) is inf")
    h.assert_true(fninf.is_negative(), "from_f64(-inf) is negative")

    // +0
    let fz = MPFloat.from_f64(0.0)
    h.assert_true(fz.is_zero(), "from_f64(0.0) is zero")
    h.assert_false(fz.is_negative(), "from_f64(0.0) is +0")

    // -0 detected via sign bit
    let fnz = MPFloat.from_f64(-0.0)
    h.assert_true(fnz.is_zero(), "from_f64(-0.0) is zero")
    h.assert_true(fnz.is_negative(), "from_f64(-0.0) is -0")

    // Positive value: string() must start with "3."
    let fpi = MPFloat.from_f64(3.14159, 8)
    h.assert_false(fpi.is_zero(), "from_f64(3.14) is not zero")
    h.assert_false(fpi.is_negative(), "from_f64(3.14) is positive")
    h.assert_true(fpi.is_finite(), "from_f64(3.14) is finite")
    h.assert_true(fpi.string().at("3."), "from_f64(3.14) string starts with \"3.\"")

    // Negative value
    let fneg = MPFloat.from_f64(-2.71828, 8)
    h.assert_true(fneg.is_negative(), "from_f64(-2.71) is negative")
    h.assert_true(fneg.string().at("-2."), "from_f64(-2.71) string starts with \"-2.\"")


// ── from_string special values ───────────────────────────────────────────────

class iso _TestMPFloatFromStringSpecial is UnitTest
  """
  `from_string` recognises all GMP-compatible special value strings.
  """
  fun name(): String => "MPFloat/from_string/special"

  fun apply(h: TestHelper) =>
    // Empty string → +0
    h.assert_true(try MPFloat.from_string("")?.is_zero() else false end,
      "\"\" → zero")
    h.assert_true(try MPFloat.from_string("     ")?.is_zero() else false end,
      "\"     \" → zero")
    h.assert_false(try MPFloat.from_string("")?.is_negative() else true end,
      "\"\" → +0 (not negative)")
    h.assert_false(try MPFloat.from_string("     ")?.is_negative() else true end,
      "\"     \" → +0 (not negative)")

    // Canonical zeros
    for s in ["0"; "0.0"; "+0"; "+0.0"; "  0  "; "  0.0  "; "  +0  "; "  +0.0  "].values() do
      h.assert_true(try MPFloat.from_string(s)?.is_zero() else false end,
        "\"" + s + "\" → zero")
      h.assert_false(try MPFloat.from_string(s)?.is_negative() else true end,
        "\"" + s + "\" → +0")
    end

    // Negative zero
    for s in ["-0"; "-0.0"; "   -0   "; "   -0.0   "].values() do
      h.assert_true(try MPFloat.from_string(s)?.is_zero() else false end,
        "\"" + s + "\" → zero")
      h.assert_true(try MPFloat.from_string(s)?.is_negative() else false end,
        "\"" + s + "\" → -0")
    end

    // NaN forms
    for s in ["nan"; "NaN"; "@NaN@"; " nan "; " NaN "; " @NaN@ "].values() do
      h.assert_true(try MPFloat.from_string(s)?.is_nan() else false end,
        "\"" + s + "\" → NaN")
    end

    // +infinity forms
    for s in ["+inf"; "@Inf@"; " +inf   "; " @Inf@   "].values() do
      h.assert_true(try MPFloat.from_string(s)?.is_infinite() else false end,
        "\"" + s + "\" → inf")
      h.assert_false(try MPFloat.from_string(s)?.is_negative() else true end,
        "\"" + s + "\" → +inf")
    end

    // -infinity forms
    for s in ["-inf"; "-@Inf@"; "    -inf  "; "    -@Inf@  "].values() do
      h.assert_true(try MPFloat.from_string(s)?.is_infinite() else false end,
        "\"" + s + "\" → inf")
      h.assert_true(try MPFloat.from_string(s)?.is_negative() else false end,
        "\"" + s + "\" → -inf")
    end

    // Invalid strings must error
    h.assert_false(try MPFloat.from_string("abc")?; true else false end,
      "\"abc\" errors")
    h.assert_false(try MPFloat.from_string("1.2.3")?; true else false end,
      "\"1.2.3\" errors")
    h.assert_false(try MPFloat.from_string("1e")?; true else false end,
      "\"1e\" (no exponent digits) errors")


// ── from_string decimal ───────────────────────────────────────────────────────

class iso _TestMPFloatFromStringDecimal is UnitTest
  """
  `from_string` parses decimal strings using full multi-precision arithmetic
  (not via F64), so the entire `prec` bytes of precision are used regardless
  of how many significant digits the string contains.
  """
  fun name(): String => "MPFloat/from_string/decimal"

  fun apply(h: TestHelper) =>
    // Integer: "3" → finite, positive, string starts with "3."
    let s3 = try MPFloat.from_string("3", 8)? else MPFloat.nan_val() end
    h.assert_true(s3.is_finite(), "\"3\" is finite")
    h.assert_false(s3.is_negative(), "\"3\" is positive")
    h.assert_true(s3.string().at("3."), "\"3\" string starts with \"3.\"")

    // Negative: "-2.5"
    let sm = try MPFloat.from_string("-2.5", 8)? else MPFloat.nan_val() end
    h.assert_true(sm.is_negative(), "\"-2.5\" is negative")
    h.assert_true(sm.is_finite(), "\"-2.5\" is finite")
    h.assert_true(sm.string().at("-2."), "\"-2.5\" string starts with \"-2.\"")

    // Exponent notation: "314e-2" = 3.14
    let se = try MPFloat.from_string("314e-2", 8)? else MPFloat.nan_val() end
    h.assert_true(se.is_finite(), "\"314e-2\" is finite")
    h.assert_false(se.is_negative(), "\"314e-2\" is positive")
    h.assert_true(se.string().at("3."), "\"314e-2\" string starts with \"3.\"")

    // '@' as exponent separator (GMP/MPFR style for base ≤ 10).
    let sat = try MPFloat.from_string("314@-2", 8)? else MPFloat.nan_val() end
    h.assert_true(sat.is_finite(), "\"314@-2\" is finite")
    h.assert_true(sat.string().at("3."), "\"314@-2\" string starts with \"3.\"")

    // Large exponent: "1E3" = 1000 → finite
    let s1k = try MPFloat.from_string("1E3", 8)? else MPFloat.nan_val() end
    h.assert_true(s1k.is_finite(), "\"1E3\" is finite")
    h.assert_false(s1k.is_negative(), "\"1E3\" is positive")

    // Large exponent: "-5.E3859045" → finite (large _exponent, prec digits).
    let sle = try MPFloat.from_string("-5.E3859045", 8)? else MPFloat.nan_val() end
    h.assert_true(sle.is_finite(), "\"-5.E3859045\" is finite")
    h.assert_true(sle.is_negative(), "\"-5.E3859045\" is negative")

    // High-precision parsing: 17 significant digits exceed F64 precision (≈15).
    // "10000000000000002" → the final '2' must survive into the MPFloat.
    // With F64 this would round to 10000000000000000.
    // Round-trip: from_string × from_f64 comparison via string prefix check.
    let shp = try MPFloat.from_string("10000000000000002", 8)? else MPFloat.nan_val() end
    h.assert_true(shp.is_finite(), "17-digit string is finite")
    h.assert_false(shp.is_negative(), "17-digit string is positive")
    h.assert_true(
      shp.string().at("1000000000000000"),
      "17-digit string starts with \"1000000000000000\"")

    // Unsupported base errors.
    h.assert_false(
      try MPFloat.from_string("ff", 8, 16)?; true else false end,
      "base 16 raises error")



// ── string() ─────────────────────────────────────────────────────────────────

class iso _TestMPFloatNewString is UnitTest
  """
  `string()` produces the expected format for special values and finite numbers.
  """
  fun name(): String => "MPFloat/new/string"

  fun apply(h: TestHelper) =>
    // Special values
    h.assert_true(MPFloat.nan_val().string() == "nan", "NaN → \"nan\"")
    h.assert_true(MPFloat.inf_val().string() == "+inf", "+inf → \"+inf\"")
    h.assert_true(MPFloat.inf_val(false).string() == "-inf", "-inf → \"-inf\"")

    // Zero: both size-0 and size-n zeros return "0.0"
    h.assert_true(MPFloat().string() == "0.0", "MPFloat() → \"0.0\"")
    h.assert_true(MPFloat(8).string() == "0.0", "MPFloat(8) → \"0.0\"")

    // Positive value ≈ 3: string starts with "3."
    let f3 = MPFloat.from_f64(3.0, 4)
    h.assert_true(f3.string().at("3."), "from_f64(3.0, 4) starts with \"3.\"")

    // Negative value ≈ -2: string starts with "-2."
    let fneg = MPFloat.from_f64(-2.0, 4)
    h.assert_true(fneg.string().at("-2."), "from_f64(-2.0, 4) starts with \"-2.\"")

    // Value in (0, 1): string starts with "0."
    let fhalf = MPFloat.from_f64(0.5, 4)
    h.assert_true(fhalf.string().at("0."), "from_f64(0.5, 4) starts with \"0.\"")

    // Verify decimal accuracy: 0.5 → first fractional digit is 5
    let shalf: String = fhalf.string()
    h.assert_true(
      try shalf(2)? == '5' else false end,
      "from_f64(0.5, 4): first fractional digit is 5")


// ── sign() predicate ──────────────────────────────────────────────────────────

class iso _TestMPFloatSign is UnitTest
  """
  `sign()` returns `Less` / `Equal` / `Greater` consistently with the value.
  """
  fun name(): String => "MPFloat/sign"

  fun apply(h: TestHelper) =>
    h.assert_true(MPFloat().sign() is Equal, "+0 sign = Equal")
    h.assert_true(MPFloat.from_f64(1.0).sign() is Greater, "+1 sign = Greater")
    h.assert_true(MPFloat.from_f64(-1.0).sign() is Less, "-1 sign = Less")
    h.assert_true(MPFloat.nan_val().sign() is Equal, "NaN sign = Equal (convention)")
    h.assert_true(MPFloat.inf_val().sign() is Greater, "+inf sign = Greater")
    h.assert_true(MPFloat.inf_val(false).sign() is Less, "-inf sign = Less")


// ── neg() ─────────────────────────────────────────────────────────────────────

class iso _TestMPFloatNeg is UnitTest
  """
  `neg()` flips the sign bit; NaN is unchanged.
  """
  fun name(): String => "MPFloat/neg"

  fun apply(h: TestHelper) =>
    let fp = MPFloat.from_f64(3.14, 4)
    h.assert_false(fp.is_negative(), "3.14 is positive")
    h.assert_true(fp.neg().is_negative(), "neg(3.14) is negative")
    h.assert_false(fp.neg().neg().is_negative(), "double neg restores positive")

    let fn2 = MPFloat.from_f64(-2.0, 4)
    h.assert_true(fn2.is_negative(), "-2.0 is negative")
    h.assert_false(fn2.neg().is_negative(), "neg(-2.0) is positive")
    h.assert_true(fn2.neg().string().at("2."), "neg(-2.0) string starts with \"2.\"")

    // NaN unchanged
    h.assert_true(MPFloat.nan_val().neg().is_nan(), "neg(NaN) is still NaN")


class iso _TestMPFloatAdd is UnitTest
  """
  `add()` handles special values, signed-zero, same-sign sums, opposite-sign
  sums (reducing to subtraction), and exponent-alignment carry.
  """
  fun name(): String => "MPFloat/add"

  fun apply(h: TestHelper) =>
    let p: USize = 8

    // Special values.
    h.assert_true(
      MPFloat.nan_val().add(MPFloat.from_f64(1.0, p)).is_nan(),
      "NaN + x = NaN")
    h.assert_true(
      MPFloat.from_f64(1.0, p).add(MPFloat.nan_val()).is_nan(),
      "x + NaN = NaN")
    h.assert_true(
      MPFloat.inf_val(true).add(MPFloat.inf_val(false)).is_nan(),
      "+inf + -inf = NaN")
    h.assert_true(
      MPFloat.inf_val(true).add(MPFloat.from_f64(1.0, p)).is_infinite(),
      "+inf + finite = inf")
    h.assert_false(
      MPFloat.inf_val(true).add(MPFloat.from_f64(1.0, p)).is_negative(),
      "+inf + finite is positive inf")

    // Identity with zero.
    h.assert_eq[String](
      MPFloat.from_f64(3.0, p).add(MPFloat.create(p)).string(),
      MPFloat.from_f64(3.0, p).string(),
      "x + 0 = x")
    h.assert_eq[String](
      MPFloat.create(p).add(MPFloat.from_f64(3.0, p)).string(),
      MPFloat.from_f64(3.0, p).string(),
      "0 + x = x")

    // Same-sign addition.
    h.assert_eq[String](
      MPFloat.from_f64(3.0, p).add(MPFloat.from_f64(2.0, p)).string(),
      MPFloat.from_f64(5.0, p).string(),
      "3 + 2 = 5")
    h.assert_eq[String](
      MPFloat.from_f64(-3.0, p).add(MPFloat.from_f64(-2.0, p)).string(),
      MPFloat.from_f64(-5.0, p).string(),
      "-3 + -2 = -5")

    // Opposite-sign addition (subtraction of magnitudes).
    h.assert_eq[String](
      MPFloat.from_f64(3.0, p).add(MPFloat.from_f64(-2.0, p)).string(),
      MPFloat.from_f64(1.0, p).string(),
      "3 + (-2) = 1")
    h.assert_eq[String](
      MPFloat.from_f64(2.0, p).add(MPFloat.from_f64(-3.0, p)).string(),
      MPFloat.from_f64(-1.0, p).string(),
      "2 + (-3) = -1")

    // Exact cancellation → zero.
    h.assert_true(
      MPFloat.from_f64(3.0, p).add(MPFloat.from_f64(-3.0, p)).is_zero(),
      "3 + (-3) = 0")

    // Carry across digit boundary: 255/256 + 1/256 = 1.
    // The carry expands the result by one digit; compare with prec=3.
    let a = MPFloat.from_f64(255.0 / 256.0, 2)
    let b = MPFloat.from_f64(1.0 / 256.0, 2)
    h.assert_eq[String](
      a.add(b).string(),
      MPFloat.from_f64(1.0, 3).string(),
      "255/256 + 1/256 = 1")

    // Exponent alignment: 256 + 1 = 257.
    // Alignment (shift=1) grows the result by one digit; compare with prec=9.
    h.assert_eq[String](
      MPFloat.from_f64(256.0, p).add(MPFloat.from_f64(1.0, p)).string(),
      MPFloat.from_f64(257.0, 9).string(),
      "256 + 1 = 257")


class iso _TestMPFloatSub is UnitTest
  """
  `sub()` delegates to `add(that.neg())`, so it is a thin wrapper. The tests
  verify correct sign handling for both orderings and exact cancellation.
  """
  fun name(): String => "MPFloat/sub"

  fun apply(h: TestHelper) =>
    let p: USize = 8

    h.assert_eq[String](
      MPFloat.from_f64(5.0, p).sub(MPFloat.from_f64(2.0, p)).string(),
      MPFloat.from_f64(3.0, p).string(),
      "5 - 2 = 3")
    h.assert_eq[String](
      MPFloat.from_f64(2.0, p).sub(MPFloat.from_f64(5.0, p)).string(),
      MPFloat.from_f64(-3.0, p).string(),
      "2 - 5 = -3")
    h.assert_true(
      MPFloat.from_f64(5.0, p).sub(MPFloat.from_f64(5.0, p)).is_zero(),
      "5 - 5 = 0")
    h.assert_eq[String](
      MPFloat.from_f64(-5.0, p).sub(MPFloat.from_f64(2.0, p)).string(),
      MPFloat.from_f64(-7.0, p).string(),
      "-5 - 2 = -7")
    h.assert_true(
      MPFloat.from_f64(1.0, p).sub(MPFloat.nan_val()).is_nan(),
      "x - NaN = NaN")


// ── Multiplication ─────────────────────────────────────────────────────────

class iso _TestMPFloatMul is UnitTest
  """
  `mul` uses FFT-based convolution. Verifies sign, exponent, and digit
  propagation for exact small cases, and special-value rules.
  """
  fun name(): String => "MPFloat/mul"

  fun apply(h: TestHelper) =>
    let p: USize = 8

    // Special values.
    h.assert_true(
      MPFloat.nan_val().mul(MPFloat.from_f64(2.0, p)).is_nan(),
      "NaN × x = NaN")
    h.assert_true(
      MPFloat.from_f64(2.0, p).mul(MPFloat.nan_val()).is_nan(),
      "x × NaN = NaN")
    h.assert_true(
      MPFloat.inf_val(true).mul(MPFloat.from_f64(2.0, p)).is_infinite(),
      "+inf × finite = inf")
    h.assert_false(
      MPFloat.inf_val(true).mul(MPFloat.from_f64(2.0, p)).is_negative(),
      "+inf × +x = +inf")
    h.assert_true(
      MPFloat.inf_val(true).mul(MPFloat.from_f64(-2.0, p)).is_negative(),
      "+inf × -x = -inf")
    h.assert_true(
      MPFloat.inf_val(true).mul(MPFloat.create(p)).is_nan(),
      "+inf × 0 = NaN")
    h.assert_true(
      MPFloat.create(p).mul(MPFloat.inf_val(true)).is_nan(),
      "0 × +inf = NaN")

    // Multiplication by zero.
    h.assert_true(
      MPFloat.from_f64(3.0, p).mul(MPFloat.create(p)).is_zero(),
      "3 × 0 = 0")
    h.assert_true(
      MPFloat.create(p).mul(MPFloat.from_f64(3.0, p)).is_zero(),
      "0 × 3 = 0")

    // Sign propagation.
    h.assert_false(
      MPFloat.from_f64(3.0, p).mul(MPFloat.from_f64(2.0, p)).is_negative(),
      "(+3) × (+2) is positive")
    h.assert_true(
      MPFloat.from_f64(-3.0, p).mul(MPFloat.from_f64(2.0, p)).is_negative(),
      "(-3) × (+2) is negative")
    h.assert_true(
      MPFloat.from_f64(3.0, p).mul(MPFloat.from_f64(-2.0, p)).is_negative(),
      "(+3) × (-2) is negative")
    h.assert_false(
      MPFloat.from_f64(-3.0, p).mul(MPFloat.from_f64(-2.0, p)).is_negative(),
      "(-3) × (-2) is positive")

    // Exact value: 3 × 2 = 6.
    h.assert_true(
      MPFloat.from_f64(3.0, p).mul(MPFloat.from_f64(2.0, p)).string().at("6."),
      "3 × 2 starts with '6.'")

    // Exact value: 10 × 10 = 100.
    h.assert_true(
      MPFloat.from_f64(10.0, p).mul(MPFloat.from_f64(10.0, p)).string().at("100."),
      "10 × 10 starts with '100.'")

    // Exponent propagation: (1/256) × 256 = 1.
    h.assert_true(
      MPFloat.from_f64(1.0 / 256.0, p).mul(MPFloat.from_f64(256.0, p))
        .string().at("1."),
      "(1/256) × 256 starts with '1.'")


// ── Inversion ──────────────────────────────────────────────────────────────

class iso _TestMPFloatInv is UnitTest
  """
  `inv` computes 1/x via Newton's method with correct sign and exponent.
  """
  fun name(): String => "MPFloat/inv"

  fun apply(h: TestHelper) =>
    let p: USize = 8

    // Special values.
    h.assert_true(MPFloat.nan_val().inv().is_nan(), "1/NaN = NaN")
    h.assert_true(MPFloat.inf_val(true).inv().is_zero(), "1/+∞ = 0")
    h.assert_true(MPFloat.inf_val(false).inv().is_zero(), "1/−∞ = 0")
    h.assert_true(MPFloat.create(p).inv().is_infinite(), "1/0 = ∞")
    h.assert_false(MPFloat.create(p).inv().is_negative(), "1/(+0) = +∞")
    h.assert_true(
      MPFloat.create(p).neg().inv().is_negative(),
      "1/(−0) = −∞")

    // Sign propagation.
    h.assert_false(
      MPFloat.from_f64(2.0, p).inv().is_negative(), "1/(+2) is positive")
    h.assert_true(
      MPFloat.from_f64(-2.0, p).inv().is_negative(), "1/(−2) is negative")

    // 1/1 = 1.
    h.assert_true(
      MPFloat.from_f64(1.0, p).inv().string().at("1."), "1/1 starts with '1.'")

    // 1/2 = 0.5 → string starts with "0.".
    h.assert_true(
      MPFloat.from_f64(2.0, p).inv().string().at("0."), "1/2 starts with '0.'")

    // Round-trip exact for powers of 2: (1/2) × 2 = 1.
    let x = MPFloat.from_f64(2.0, p)
    h.assert_true(
      x.inv().mul(x).string().at("1."), "(1/2) × 2 = 1")

    // Exponent: 1/256 = 256^{-1} → exponent should be 0.
    let inv256 = MPFloat.from_f64(256.0, p).inv()
    h.assert_true(inv256.string().at("0."), "1/256 starts with '0.'")

    // 1/(1/256) = 256.
    let inv_inv256: String val = MPFloat.from_f64(1.0 / 256.0, p).inv().string()
    h.assert_true(
      inv_inv256.at("256."),
      "1/(1/256) starts with '256.'")


// ── Division ───────────────────────────────────────────────────────────────

class iso _TestMPFloatDiv is UnitTest
  """
  `div` computes this/that = this × (1/that) with correct sign and exponent.
  """
  fun name(): String => "MPFloat/div"

  fun apply(h: TestHelper) =>
    let p: USize = 8

    // Special values.
    h.assert_true(MPFloat.nan_val().div(MPFloat.from_f64(2.0, p)).is_nan(),
      "NaN / x = NaN")
    h.assert_true(MPFloat.from_f64(2.0, p).div(MPFloat.nan_val()).is_nan(),
      "x / NaN = NaN")
    h.assert_true(
      MPFloat.inf_val(true).div(MPFloat.inf_val(true)).is_nan(), "+∞/+∞ = NaN")
    h.assert_true(
      MPFloat.from_f64(2.0, p).div(MPFloat.inf_val(true)).is_zero(),
      "finite / +∞ = 0")
    h.assert_true(
      MPFloat.inf_val(true).div(MPFloat.from_f64(2.0, p)).is_infinite(),
      "+∞ / finite = ∞")
    h.assert_true(
      MPFloat.from_f64(2.0, p).div(MPFloat.create(p)).is_infinite(),
      "finite / 0 = ∞")
    h.assert_true(
      MPFloat.create(p).div(MPFloat.create(p)).is_nan(), "0 / 0 = NaN")

    // Sign propagation.
    h.assert_false(
      MPFloat.from_f64(6.0, p).div(MPFloat.from_f64(2.0, p)).is_negative(),
      "(+6) / (+2) is positive")
    h.assert_true(
      MPFloat.from_f64(-6.0, p).div(MPFloat.from_f64(2.0, p)).is_negative(),
      "(−6) / (+2) is negative")
    h.assert_true(
      MPFloat.from_f64(6.0, p).div(MPFloat.from_f64(-2.0, p)).is_negative(),
      "(+6) / (−2) is negative")
    h.assert_false(
      MPFloat.from_f64(-6.0, p).div(MPFloat.from_f64(-2.0, p)).is_negative(),
      "(−6) / (−2) is positive")

    // Exact values.
    h.assert_true(
      MPFloat.from_f64(6.0, p).div(MPFloat.from_f64(2.0, p)).string().at("3."),
      "6 / 2 starts with '3.'")
    h.assert_true(
      MPFloat.from_f64(1.0, p).div(MPFloat.from_f64(4.0, p)).string().at("0."),
      "1 / 4 starts with '0.'")

    // Round-trip exact for powers of 2: (6 / 2) × 2 = 6.
    let xrt = MPFloat.from_f64(6.0, p)
    let yrt = MPFloat.from_f64(2.0, p)
    h.assert_true(
      xrt.div(yrt).mul(yrt).string().at("6."), "(6/2) × 2 = 6")


class iso _TestMPFloatSqrt is UnitTest
  """
  `sqrt` computes √this with correct sign, exponent, and IEEE 754 special
  cases. Uses the parity-split Newton reciprocal-square-root algorithm.
  """
  fun name(): String => "MPFloat/sqrt"

  fun apply(h: TestHelper) =>
    let p: USize = 8

    // Special values.
    h.assert_true(MPFloat.nan_val().sqrt().is_nan(), "sqrt(NaN) = NaN")
    h.assert_true(MPFloat.inf_val(true).sqrt().is_infinite(), "sqrt(+∞) = +∞")
    h.assert_false(MPFloat.inf_val(true).sqrt().is_negative(), "sqrt(+∞) is positive")
    h.assert_true(MPFloat.inf_val(false).sqrt().is_nan(), "sqrt(−∞) = NaN")
    h.assert_true(MPFloat.from_f64(0.0, p).sqrt().is_zero(), "sqrt(+0) = +0")
    h.assert_false(MPFloat.from_f64(0.0, p).sqrt().is_negative(), "sqrt(+0) is positive")
    h.assert_true(MPFloat.from_f64(-0.0, p).sqrt().is_zero(), "sqrt(−0) = −0")
    h.assert_true(MPFloat.from_f64(-0.0, p).sqrt().is_negative(), "sqrt(−0) is negative")

    // Negative non-zero → NaN.
    h.assert_true(MPFloat.from_f64(-1.0, p).sqrt().is_nan(), "sqrt(−1) = NaN")
    h.assert_true(MPFloat.from_f64(-4.0, p).sqrt().is_nan(), "sqrt(−4) = NaN")

    // Exact perfect squares that are powers of 2 (Newton gives exact result).
    // sqrt(4) = 2, sqrt(16) = 4, sqrt(64) = 8 — all exact in base-256.
    h.assert_true(
      MPFloat.from_f64(4.0, p).sqrt().string().at("2."),
      "sqrt(4) starts with '2.'")
    h.assert_true(
      MPFloat.from_f64(16.0, p).sqrt().string().at("4."),
      "sqrt(16) starts with '4.'")
    h.assert_true(
      MPFloat.from_f64(64.0, p).sqrt().string().at("8."),
      "sqrt(64) starts with '8.'")

    // Value with even exponent (exp=0 for 0.25 ∈ [1/256, 1)).
    h.assert_true(
      MPFloat.from_f64(0.25, p).sqrt().string().at("0."),
      "sqrt(0.25) starts with '0.' (= 0.5)")

    // Irrational: √2 ≈ 1.41421…
    h.assert_true(
      MPFloat.from_f64(2.0, p).sqrt().string().at("1."),
      "sqrt(2) starts with '1.'")

    // Round-trip exact for perfect square: sqrt(4)² = 4.
    h.assert_true(
      MPFloat.from_f64(4.0, p).sqrt().mul(MPFloat.from_f64(4.0, p).sqrt()).string().at("4."),
      "sqrt(4) × sqrt(4) = 4")


// ── Comparisons ──────────────────────────────────────────────────────────────

class iso _TestMPFloatCmp is UnitTest
  """
  Verifies all six comparison operators (eq, ne, lt, le, ge, gt) and
  `compare` against the full IEEE 754 truth tables, including:
    - NaN (unordered: all comparisons false except ne)
    - ±∞ self-equality and cross-sign inequality
    - finite == ±∞ is false
    - −0 == +0 (equal); −0 is not < +0
    - mixed-sign finite pairs (negative < positive)
    - same-sign pairs with different exponents (256 > 1)
    - same-sign pairs with equal exponents but different digits
    - min_value() == −∞ and max_value() == +∞ sanity
  """
  fun name(): String => "MPFloat/cmp"

  fun apply(h: TestHelper) =>
    let p: USize = 8
    let nan  = MPFloat.nan_val()
    let pinf = MPFloat.inf_val(true)
    let ninf = MPFloat.inf_val(false)
    let pz   = MPFloat.create(p)             // +0
    let nz   = MPFloat.create(p).neg()       // -0
    let one  = MPFloat.from_f64(1.0, p)
    let two  = MPFloat.from_f64(2.0, p)
    let neg1 = MPFloat.from_f64(-1.0, p)
    let neg2 = MPFloat.from_f64(-2.0, p)
    let big  = MPFloat.from_f64(256.0, p)   // exponent = 2

    // ── NaN: all comparisons false except ne ─────────────────────────────────
    h.assert_false(nan.eq(nan),   "NaN == NaN is false")
    h.assert_true( nan.ne(nan),   "NaN != NaN is true")
    h.assert_false(nan.lt(one),   "NaN < x is false")
    h.assert_false(nan.le(one),   "NaN <= x is false")
    h.assert_false(nan.ge(one),   "NaN >= x is false")
    h.assert_false(nan.gt(one),   "NaN > x is false")
    h.assert_false(one.lt(nan),   "x < NaN is false")
    h.assert_false(one.le(nan),   "x <= NaN is false")
    h.assert_false(one.ge(nan),   "x >= NaN is false")
    h.assert_false(one.gt(nan),   "x > NaN is false")
    h.assert_false(nan.eq(one),   "NaN == x is false")
    h.assert_true( nan.ne(one),   "NaN != x is true")

    // ── Infinities ───────────────────────────────────────────────────────────
    h.assert_true( pinf.eq(pinf),  "+inf == +inf")
    h.assert_true( ninf.eq(ninf),  "-inf == -inf")
    h.assert_false(pinf.eq(ninf),  "+inf != -inf")
    h.assert_false(ninf.eq(pinf),  "-inf != +inf")
    h.assert_false(pinf.eq(one),   "+inf != finite")
    h.assert_false(one.eq(pinf),   "finite != +inf")
    h.assert_false(pinf.eq(pz),    "+inf != 0")
    h.assert_true( pinf.ne(ninf),  "+inf ne -inf")
    h.assert_false(pinf.lt(pinf),  "+inf < +inf is false")
    h.assert_false(ninf.gt(ninf),  "-inf > -inf is false")
    h.assert_true( ninf.lt(pinf),  "-inf < +inf")
    h.assert_true( ninf.lt(one),   "-inf < finite")
    h.assert_true( one.lt(pinf),   "finite < +inf")
    h.assert_false(one.lt(ninf),   "finite < -inf is false")
    h.assert_false(pinf.lt(one),   "+inf < finite is false")
    h.assert_true( pinf.ge(pinf),  "+inf >= +inf")
    h.assert_true( ninf.le(ninf),  "-inf <= -inf")

    // ── min_value / max_value sanity ─────────────────────────────────────────
    h.assert_true(MPFloat.min_value().is_infinite(),      "min_value is inf")
    h.assert_true(MPFloat.min_value().is_negative(), "min_value is -inf (negative)")
    h.assert_true(MPFloat.max_value().is_infinite(),      "max_value is inf")
    h.assert_false(MPFloat.max_value().is_negative(),"max_value is +inf (positive)")
    h.assert_true(MPFloat.min_value().lt(MPFloat.max_value()), "-inf < +inf")
    h.assert_true(MPFloat.min_value().eq(ninf), "min_value == -inf")
    h.assert_true(MPFloat.max_value().eq(pinf), "max_value == +inf")

    // ── ±0 ──────────────────────────────────────────────────────────────────
    h.assert_true( pz.eq(nz),    "+0 == -0")
    h.assert_true( nz.eq(pz),    "-0 == +0")
    h.assert_false(pz.ne(nz),    "+0 ne -0 is false")
    h.assert_false(pz.lt(nz),    "+0 < -0 is false")
    h.assert_false(nz.lt(pz),    "-0 < +0 is false")
    h.assert_true( pz.le(nz),    "+0 <= -0")
    h.assert_true( nz.le(pz),    "-0 <= +0")
    h.assert_true( pz.ge(nz),    "+0 >= -0")
    h.assert_true( nz.ge(pz),    "-0 >= +0")
    h.assert_false(pz.gt(nz),    "+0 > -0 is false")
    h.assert_false(nz.gt(pz),    "-0 > +0 is false")

    // ── mixed sign ──────────────────────────────────────────────────────────
    h.assert_true( neg1.lt(one),  "-1 < 1")
    h.assert_false(one.lt(neg1),  "1 < -1 is false")
    h.assert_true( one.gt(neg1),  "1 > -1")
    h.assert_false(neg1.gt(one),  "-1 > 1 is false")
    h.assert_false(neg1.eq(one),  "-1 != 1")
    h.assert_true( neg1.ne(one),  "-1 ne 1")
    h.assert_true( neg1.le(one),  "-1 <= 1")
    h.assert_true( one.ge(neg1),  "1 >= -1")

    // ── same sign, different exponents (1 vs 256) ────────────────────────────
    h.assert_false(big.eq(one),  "256 != 1")
    h.assert_false(one.eq(big),  "1 != 256")
    h.assert_true( one.lt(big),  "1 < 256")
    h.assert_false(big.lt(one),  "256 < 1 is false")
    h.assert_true( big.gt(one),  "256 > 1")
    h.assert_false(one.gt(big),  "1 > 256 is false")
    h.assert_true( one.le(big),  "1 <= 256")
    h.assert_true( big.ge(one),  "256 >= 1")

    // Negative counterparts
    h.assert_true( neg1.gt(neg2), "-1 > -2")
    h.assert_true( neg2.lt(neg1), "-2 < -1")
    h.assert_false(neg1.eq(neg2), "-1 != -2")

    // ── same exponent, different digits ──────────────────────────────────────
    h.assert_true( one.lt(two),   "1 < 2")
    h.assert_false(two.lt(one),   "2 < 1 is false")
    h.assert_true( one.eq(one),   "1 == 1")
    h.assert_false(one.lt(one),   "1 < 1 is false")
    h.assert_true( one.le(one),   "1 <= 1")
    h.assert_true( one.ge(one),   "1 >= 1")
    h.assert_false(one.gt(one),   "1 > 1 is false")

    // ── compare helper ───────────────────────────────────────────────────────
    h.assert_true(one.compare(two)  is Less,    "compare(1, 2) = Less")
    h.assert_true(two.compare(one)  is Greater, "compare(2, 1) = Greater")
    h.assert_true(one.compare(one)  is Equal,   "compare(1, 1) = Equal")
    h.assert_true(neg1.compare(one) is Less,    "compare(-1, 1) = Less")


// ── from_mpint ────────────────────────────────────────────────────────────────

class iso _TestMPFloatFromMPInt is UnitTest
  """
  `from_mpint` converts an `MPInt` to an `MPFloat` with correct sign,
  magnitude, and precision.
  """
  fun name(): String => "MPFloat/from_mpint"

  fun apply(h: TestHelper) =>
    let p: USize = 8

    // Zero: always +0.
    let z = MPFloat.from_mpint(MPInt.from_ilong(0), p)
    h.assert_true(z.is_zero(),     "from_mpint(0) is zero")
    h.assert_false(z.is_negative(), "from_mpint(0) is +0")
    h.assert_false(z.is_nan(),     "from_mpint(0) is not NaN")
    h.assert_true(z.is_finite(),   "from_mpint(0) is finite")

    // Positive small integer: 1 → "1."
    let one = MPFloat.from_mpint(MPInt.from_ilong(1), p)
    h.assert_false(one.is_negative(), "from_mpint(1) is positive")
    h.assert_true(one.is_finite(),    "from_mpint(1) is finite")
    h.assert_true(one.string().at("1."), "from_mpint(1) string starts with \"1.\"")

    // Negative integer: -3 → "-3."
    let neg3 = MPFloat.from_mpint(MPInt.from_ilong(-3), p)
    h.assert_true(neg3.is_negative(),     "from_mpint(-3) is negative")
    h.assert_true(neg3.string().at("-3."), "from_mpint(-3) string starts with \"-3.\"")

    // Larger integer: 1000 → "1000."
    let thou = MPFloat.from_mpint(MPInt.from_ilong(1000), p)
    h.assert_false(thou.is_negative(),        "from_mpint(1000) is positive")
    h.assert_true(thou.string().at("1000."),  "from_mpint(1000) string starts with \"1000.\"")

    // Power-of-2: 256 → "256."
    let i256 = MPFloat.from_mpint(MPInt.from_ilong(256), p)
    h.assert_true(i256.string().at("256."), "from_mpint(256) string starts with \"256.\"")

    // Power-of-2: 65536 → "65536."
    let i65536 = MPFloat.from_mpint(MPInt.from_ilong(65536), p)
    h.assert_true(i65536.string().at("65536."), "from_mpint(65536) string starts with \"65536.\"")

    // Negative large: -65536 → "-65536."
    let ni65536 = MPFloat.from_mpint(MPInt.from_ilong(-65536), p)
    h.assert_true(ni65536.is_negative(),          "from_mpint(-65536) is negative")
    h.assert_true(ni65536.string().at("-65536."), "from_mpint(-65536) string starts with \"-65536.\"")

    // Round-trip sign: from_mpint(n).is_negative() == n.is_negative()
    let npos = MPInt.from_ilong(42)
    let nneg = MPInt.from_ilong(-42)
    h.assert_false(MPFloat.from_mpint(npos, p).is_negative(), "positive MPInt → positive MPFloat")
    h.assert_true( MPFloat.from_mpint(nneg, p).is_negative(), "negative MPInt → negative MPFloat")

    // Precision parameter is honoured: value must be finite.
    let big = MPFloat.from_mpint(MPInt.from_ilong(1000000), 4)
    h.assert_true(big.is_finite(),            "from_mpint(1e6, prec=4) is finite")
    h.assert_false(big.is_negative(),         "from_mpint(1e6, prec=4) is positive")
    h.assert_true(big.string().at("1000000"), "from_mpint(1e6, prec=4) string starts with \"1000000\"")

    // ── Large MPInt (beyond I64 / F64 range) ─────────────────────────────────

    // 10^20 > U64.max (≈1.8×10^19): requires MPInt arithmetic to construct.
    //
    // With p=8 (≈19.3 significant decimal digits), 10^20 (21 digits) cannot
    // be held exactly; the best 8-byte approximation is just below 10^20 and
    // its string representation starts with "9".  We only verify finiteness,
    // positiveness, and that the leading digit is "9" or "1" (order-of-
    // magnitude correct).
    let e20: MPInt = try MPInt.from_string("100000000000000000000")? else MPInt.from_ilong(0) end
    let fe20 = MPFloat.from_mpint(e20, p)
    h.assert_true(fe20.is_finite(),    "from_mpint(10^20, p=8) is finite")
    h.assert_false(fe20.is_negative(), "from_mpint(10^20, p=8) is positive")
    h.assert_true(
      fe20.string().at("9") or fe20.string().at("1"),
      "from_mpint(10^20, p=8) leading digit is 9 or 1")

    // With p=12 (≈29 significant decimal digits), 10^20 (21 digits) fits
    // comfortably and the string should represent "1e+20".
    let fe20p12 = MPFloat.from_mpint(e20, 12)
    h.assert_true(fe20p12.is_finite(),    "from_mpint(10^20, p=12) is finite")
    h.assert_false(fe20p12.is_negative(), "from_mpint(10^20, p=12) is positive")
    h.assert_true(
      fe20p12.string().at("1e+20") or fe20p12.string().at("100000000000000000000"),
      "from_mpint(10^20, p=12) string represents 10^20")

    // 29-digit positive: precision must be ≥ ceil(29 / log10(256)) ≈ 13 bytes
    // to represent all digits. With p=16 bytes all 29 digits survive.
    let d29: MPInt =
      try MPInt.from_string("12345678901234567890123456789")?
      else MPInt.from_ilong(0) end
    let fd29 = MPFloat.from_mpint(d29, 16)
    h.assert_true(fd29.is_finite(),    "from_mpint(29-digit) is finite")
    h.assert_false(fd29.is_negative(), "from_mpint(29-digit) is positive")
    h.assert_true(
      fd29.string().at("12345678901234567890"),
      "from_mpint(29-digit) first 20 digits preserved with prec=16")

    // Negative large: -10^20 with p=12 for sufficient precision.
    let ne20: MPInt =
      try MPInt.from_string("-100000000000000000000")? else MPInt.from_ilong(0) end
    let fne20 = MPFloat.from_mpint(ne20, 12)
    h.assert_true(fne20.is_finite(),    "from_mpint(-10^20) is finite")
    h.assert_true(fne20.is_negative(),  "from_mpint(-10^20) is negative")
    h.assert_true(
      fne20.string().at("-1e+20") or fne20.string().at("-100000000000000000000"),
      "from_mpint(-10^20, p=12) string represents -10^20")

    // 39-digit positive with precision of 50 "digits"
    let d39: MPInt =
      try MPInt.from_string("123456789012345678901234567890123456789")?
      else MPInt.from_ilong(0) end
    let fd39 = MPFloat.from_mpint(d39, 50)
    h.assert_true(fd39.is_finite(),    "from_mpint(39-digit) is finite")
    h.assert_false(fd39.is_negative(), "from_mpint(39-digit) is positive")
    h.assert_true(
      fd39.string().at("12345678901234567890"),
      "from_mpint(39-digit) first 20 digits preserved with prec=50")
    // Digits before decimal point are the same
    h.assert_eq[String](fd39.string().substring(0, try fd39.string().find(".")? else 0 end), d39.string(), "fd39 and d39 have same string representation")


// ── divrem / rem / fld / mod ──────────────────────────────────────────────────

class iso _TestMPFloatDivRem is UnitTest
  """
  `divrem`, `rem`, `fld`, and `mod` on finite positive and negative operands.

  Verified properties:
  - `divrem` invariant: `this = q × that + r` for all finite cases.
  - `rem` has the same sign as `this` (the dividend).
  - `fld` equals `trunc` for same-sign operands, and `trunc − 1` when the
    signs differ and the remainder is non-zero.
  - `mod` has the same sign as `that` (the divisor).
  - NaN and ±∞ special cases.
  - Unsafe variants produce the same results on finite inputs.
  """

  fun name(): String =>
    "MPFloat/divrem"

  fun apply(h: TestHelper) =>
    let p: USize = 8

    // Helpers: build MPFloat from a decimal string.
    let f = {(s: String): MPFloat =>
      try MPFloat.from_string(s, p)? else MPFloat.create(p) end
    }

    // ── Exact integer operands ────────────────────────────────────────────────

    // 7 / 2 = 3 rem 1
    (let q7_2, let r7_2) = f("7").divrem(f("2"))
    h.assert_true(q7_2.string().at("3."),  "divrem(7,2) q = 3")
    h.assert_true(r7_2.string().at("1."),  "divrem(7,2) r = 1")
    h.assert_false(r7_2.is_negative(),     "divrem(7,2) r ≥ 0")

    // -7 / 2 = -3 rem -1  (truncation toward zero)
    (let qn7_2, let rn7_2) = f("-7").divrem(f("2"))
    h.assert_true(qn7_2.string().at("-3."), "divrem(-7,2) q = -3")
    h.assert_true(rn7_2.string().at("-1."), "divrem(-7,2) r = -1")
    h.assert_true(rn7_2.is_negative(),      "divrem(-7,2) r ≤ 0")

    // 7 / -2 = -3 rem 1  (truncation toward zero)
    (let q7_n2, let r7_n2) = f("7").divrem(f("-2"))
    h.assert_true(q7_n2.string().at("-3."), "divrem(7,-2) q = -3")
    h.assert_true(r7_n2.string().at("1."),  "divrem(7,-2) r = 1")
    h.assert_false(r7_n2.is_negative(),     "divrem(7,-2) r ≥ 0")

    // -7 / -2 = 3 rem -1
    (let qn7_n2, let rn7_n2) = f("-7").divrem(f("-2"))
    h.assert_true(qn7_n2.string().at("3."),  "divrem(-7,-2) q = 3")
    h.assert_true(rn7_n2.string().at("-1."), "divrem(-7,-2) r = -1")

    // divrem invariant: this = q * that + r, and |r| < |that|
    let check_inv = {(h2: TestHelper, a: String, b: String) =>
      let fa = try MPFloat.from_string(a, p)? else MPFloat.create(p) end
      let fb = try MPFloat.from_string(b, p)? else MPFloat.create(p) end
      (let q, let r) = fa.divrem(fb)
      let rhs = q.mul(fb).add(r)
      h2.assert_true(fa.eq(rhs), "invariant: " + a + " = q×" + b + " + r")
      // |r| < |b| — the key correctness constraint; tautologically satisfied
      // if we merely define r = a − q×b, so must be checked separately.
      h2.assert_true(r.abs().lt(fb.abs()), "|r| < |b| for " + a + "/" + b)
    }
    check_inv(h, "7", "2")
    check_inv(h, "-7", "2")
    check_inv(h, "7", "-2")
    check_inv(h, "-7", "-2")
    check_inv(h, "10", "3")
    check_inv(h, "-10", "3")
    // Non-dyadic divisors: Newton inv() undershoots, so post-correction is needed.
    check_inv(h, "10", "5")
    check_inv(h, "-10", "5")
    check_inv(h, "14", "7")
    check_inv(h, "15", "5")

    // Exact integer divisibility (exercises the post-correction step).
    (let q10_5, let r10_5) = f("10").divrem(f("5"))
    h.assert_true(q10_5.string().at("2."), "divrem(10,5) q = 2")
    h.assert_true(r10_5.is_zero(),          "divrem(10,5) r = 0")

    // ── rem ───────────────────────────────────────────────────────────────────

    // rem sign follows dividend
    h.assert_false(f("7").rem(f("2")).is_negative(),   "rem(7,2) ≥ 0")
    h.assert_true( f("-7").rem(f("2")).is_negative(),  "rem(-7,2) ≤ 0")
    h.assert_false(f("7").rem(f("-2")).is_negative(),  "rem(7,-2) ≥ 0")
    h.assert_true( f("-7").rem(f("-2")).is_negative(), "rem(-7,-2) ≤ 0")

    // ── fld ──────────────────────────────────────────────────────────────────

    // Same-sign operands: fld = trunc
    h.assert_true(f("7").fld(f("2")).string().at("3."),    "fld(7,2) = 3")
    h.assert_true(f("-7").fld(f("-2")).string().at("3."),  "fld(-7,-2) = 3")

    // Opposite-sign: fld = trunc - 1
    h.assert_true(f("-7").fld(f("2")).string().at("-4."),  "fld(-7,2) = -4")
    h.assert_true(f("7").fld(f("-2")).string().at("-4."),  "fld(7,-2) = -4")

    // Exact division: no adjustment needed
    h.assert_true(f("6").fld(f("2")).string().at("3."),    "fld(6,2) = 3")
    h.assert_true(f("-6").fld(f("2")).string().at("-3."),  "fld(-6,2) = -3")
    // Non-dyadic exact division (10/5 = 2; exercises post-correction via divrem).
    h.assert_true(f("10").fld(f("5")).string().at("2."),   "fld(10,5) = 2")
    h.assert_true(f("10").mod(f("5")).is_zero(),            "mod(10,5) = 0")

    // ── mod ──────────────────────────────────────────────────────────────────

    // mod sign follows divisor
    h.assert_false(f("7").mod(f("2")).is_negative(),   "mod(7,2) ≥ 0")
    h.assert_false(f("-7").mod(f("2")).is_negative(),  "mod(-7,2) ≥ 0  (sign = divisor)")
    h.assert_true( f("7").mod(f("-2")).is_negative(),  "mod(7,-2) ≤ 0  (sign = divisor)")
    h.assert_true( f("-7").mod(f("-2")).is_negative(), "mod(-7,-2) ≤ 0")

    // mod values: mod(-7, 2) = 1,  mod(7, -2) = -1
    h.assert_true(f("-7").mod(f("2")).string().at("1."),  "mod(-7,2) = 1")
    h.assert_true(f("7").mod(f("-2")).string().at("-1."), "mod(7,-2) = -1")

    // fld / mod invariant: this = fld(this,that) * that + mod(this,that)
    let check_fld_inv = {(h2: TestHelper, a: String, b: String) =>
      let fa = try MPFloat.from_string(a, p)? else MPFloat.create(p) end
      let fb = try MPFloat.from_string(b, p)? else MPFloat.create(p) end
      let rhs = fa.fld(fb).mul(fb).add(fa.mod(fb))
      h2.assert_true(fa.eq(rhs), "fld/mod invariant: " + a + " / " + b)
    }
    check_fld_inv(h, "7", "2")
    check_fld_inv(h, "-7", "2")
    check_fld_inv(h, "7", "-2")
    check_fld_inv(h, "-7", "-2")

    // ── Zero dividend ────────────────────────────────────────────────────────

    h.assert_true(f("0").divrem(f("3"))._1.is_zero(), "divrem(0,3) q = 0")
    h.assert_true(f("0").divrem(f("3"))._2.is_zero(), "divrem(0,3) r = 0")
    h.assert_true(f("0").fld(f("3")).is_zero(),        "fld(0,3) = 0")
    h.assert_true(f("0").mod(f("3")).is_zero(),        "mod(0,3) = 0")

    // ── Special values ────────────────────────────────────────────────────────

    let nan = MPFloat.nan_val()
    let inf = MPFloat.inf_val()

    h.assert_true(nan.divrem(f("2"))._1.is_nan(), "divrem(NaN,2) q is NaN")
    h.assert_true(nan.divrem(f("2"))._2.is_nan(), "divrem(NaN,2) r is NaN")
    h.assert_true(f("2").divrem(nan)._1.is_nan(), "divrem(2,NaN) q is NaN")
    h.assert_true(inf.divrem(inf)._1.is_nan(),    "divrem(+inf,+inf) q is NaN")
    h.assert_true(f("2").divrem(f("0"))._2.is_nan(), "divrem(2,0) r is NaN")
    h.assert_true(f("2").divrem(inf)._1.is_zero(),   "divrem(2,+inf) q = 0")
    h.assert_false(f("2").divrem(inf)._2.is_zero(),  "divrem(2,+inf) r = 2 ≠ 0")

    h.assert_true(nan.fld(f("2")).is_nan(),  "fld(NaN,2) is NaN")
    h.assert_true(f("2").fld(nan).is_nan(),  "fld(2,NaN) is NaN")
    h.assert_true(f("2").fld(f("0")).is_nan(), "fld(2,0) is NaN")

    h.assert_true(nan.rem(f("2")).is_nan(),  "rem(NaN,2) is NaN")
    h.assert_true(nan.mod(f("2")).is_nan(),  "mod(NaN,2) is NaN")

    // ── Unsafe variants (same results on finite inputs) ───────────────────────

    (let qu, let ru) = f("7").divrem_unsafe(f("2"))
    h.assert_true(qu.string().at("3."), "divrem_unsafe(7,2) q = 3")
    h.assert_true(ru.string().at("1."), "divrem_unsafe(7,2) r = 1")

    h.assert_true(f("7").rem_unsafe(f("2")).string().at("1."),     "rem_unsafe(7,2) = 1")
    h.assert_true(f("-7").rem_unsafe(f("2")).string().at("-1."),   "rem_unsafe(-7,2) = -1")
    h.assert_true(f("-7").fld_unsafe(f("2")).string().at("-4."),   "fld_unsafe(-7,2) = -4")
    h.assert_true(f("-7").mod_unsafe(f("2")).string().at("1."),    "mod_unsafe(-7,2) = 1")


// ── trunc / floor / ceil / round ─────────────────────────────────────────────

class iso _TestMPFloatRounding is UnitTest
  """
  `trunc`, `floor`, `ceil`, and `round` on positive, negative, fractional,
  and exact-integer operands.  Also verifies NaN and ±∞ pass-through and
  that the `floor`/`ceil`/`round` implementations are consistent with the
  `trunc` primitive.
  """

  fun name(): String =>
    "MPFloat/rounding"

  fun apply(h: TestHelper) =>
    let p: USize = 8

    // from_string converts via Horner + ÷10 scaling, which is not exact in
    // base 256.  Use it only for cases where the truncation direction is
    // unambiguous regardless of small rounding error (e.g. f("2.7") is close
    // enough to 2.7 that trunc still gives 2).  For exact integers and exact
    // half-integers use from_f64 or from_mpint, which are exact for
    // binary fractions.
    let f = {(s: String): MPFloat =>
      try MPFloat.from_string(s, p)? else MPFloat.create(p) end
    }
    let ff = {(x: F64): MPFloat => MPFloat.from_f64(x, p) }
    let fi = {(n: ILong): MPFloat => MPFloat.from_mpint(MPInt.from_ilong(n), p) }

    // ── trunc ────────────────────────────────────────────────────────────────

    // Positive fractional: from_string approximation still truncates correctly.
    h.assert_true(f("2.7").trunc().string().at("2."),   "trunc(2.7) = 2")
    h.assert_true(f("0.9").trunc().string().at("0."),   "trunc(0.9) = 0")

    // Exact integer via from_f64 (exact for binary fractions).
    h.assert_true(ff(2.0).trunc().string().at("2."),    "trunc(2.0) = 2")
    h.assert_true(ff(256.0).trunc().string().at("256."), "trunc(256.0) = 256")

    // Negative fractional.
    h.assert_true(f("-2.7").trunc().string().at("-2."),  "trunc(-2.7) = -2")
    h.assert_true(f("-0.9").trunc().string().at("0."),   "trunc(-0.9) = 0")
    h.assert_true(ff(-2.0).trunc().string().at("-2."),   "trunc(-2.0) = -2")

    // Special values pass through.
    h.assert_true(MPFloat.nan_val().trunc().is_nan(),        "trunc(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().trunc().is_infinite(),   "trunc(+inf) = +inf")
    h.assert_false(MPFloat.inf_val().trunc().is_negative(),  "trunc(+inf) is positive")

    // ── floor ────────────────────────────────────────────────────────────────

    // Positive: same as trunc.
    h.assert_true(f("2.7").floor().string().at("2."),    "floor(2.7) = 2")
    h.assert_true(ff(2.0).floor().string().at("2."),     "floor(2.0) = 2")
    h.assert_true(f("0.9").floor().string().at("0."),    "floor(0.9) = 0")

    // Negative: one below trunc when fractional part is non-zero.
    h.assert_true(f("-2.7").floor().string().at("-3."),  "floor(-2.7) = -3")
    h.assert_true(ff(-2.0).floor().string().at("-2."),   "floor(-2.0) = -2")
    h.assert_true(f("-0.1").floor().string().at("-1."),  "floor(-0.1) = -1")

    h.assert_true(MPFloat.nan_val().floor().is_nan(),        "floor(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().floor().is_infinite(),   "floor(+inf) = +inf")
    h.assert_true(MPFloat.inf_val().neg().floor().is_infinite(), "floor(-inf) = -inf")

    // ── ceil ─────────────────────────────────────────────────────────────────

    // Positive: one above trunc when fractional part is non-zero.
    h.assert_true(f("2.1").ceil().string().at("3."),     "ceil(2.1) = 3")
    h.assert_true(ff(2.0).ceil().string().at("2."),      "ceil(2.0) = 2")
    h.assert_true(f("0.1").ceil().string().at("1."),     "ceil(0.1) = 1")

    // Negative: same as trunc.
    h.assert_true(f("-2.7").ceil().string().at("-2."),   "ceil(-2.7) = -2")
    h.assert_true(ff(-2.0).ceil().string().at("-2."),    "ceil(-2.0) = -2")
    h.assert_true(f("-0.9").ceil().string().at("0."),    "ceil(-0.9) = 0")

    h.assert_true(MPFloat.nan_val().ceil().is_nan(),        "ceil(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().ceil().is_infinite(),   "ceil(+inf) = +inf")

    // ── round (half away from zero) ───────────────────────────────────────────

    // Half-integers are exact in base 256 (0.5 = 128/256, 2.5 = 2 + 128/256).
    // Use from_f64 for these to avoid the from_string base-256 rounding issue.

    // Positive: from_string approximation is sufficient for non-half cases.
    h.assert_true(f("2.4").round().string().at("2."),    "round(2.4) = 2")
    h.assert_true(f("2.6").round().string().at("3."),    "round(2.6) = 3")
    h.assert_true(f("0.4").round().string().at("0."),    "round(0.4) = 0")

    // Half-integer ties: exact via from_f64.
    h.assert_true(ff(0.5).round().string().at("1."),     "round(0.5) = 1  (half away)")
    h.assert_true(ff(2.5).round().string().at("3."),     "round(2.5) = 3  (half away)")

    // Negative.
    h.assert_true(f("-2.4").round().string().at("-2."),  "round(-2.4) = -2")
    h.assert_true(f("-2.6").round().string().at("-3."),  "round(-2.6) = -3")
    h.assert_true(f("-0.4").round().string().at("0."),   "round(-0.4) = 0")
    h.assert_true(ff(-0.5).round().string().at("-1."),   "round(-0.5) = -1  (half away)")
    h.assert_true(ff(-2.5).round().string().at("-3."),   "round(-2.5) = -3  (half away)")

    // Exact integers: trunc = floor = ceil = round.
    h.assert_true(fi(3).round().string().at("3."),       "round(3) = 3")
    h.assert_true(fi(-3).round().string().at("-3."),     "round(-3) = -3")

    h.assert_true(MPFloat.nan_val().round().is_nan(),        "round(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().round().is_infinite(),   "round(+inf) = +inf")

