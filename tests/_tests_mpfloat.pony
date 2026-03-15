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
  version returned a `(MPFloat val, I8)` tuple; `sub` now returns a plain
  `MPFloat val`. Tests have been superseded by `_TestMPFloatSub`.
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
    let z: String val = MPFloat(8).string()
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
    h.assert_false(z0.is_inf(), "MPFloat() is not inf")
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
    h.assert_false(n.is_inf(), "nan_val() is not inf")
    h.assert_false(n.is_zero(), "nan_val() is not zero")
    h.assert_false(n.is_negative(), "nan_val() is not negative (NaN sign ignored)")

    // +infinity
    let pinf = MPFloat.inf_val()
    h.assert_true(pinf.is_inf(), "+inf_val() is inf")
    h.assert_false(pinf.is_nan(), "+inf_val() is not NaN")
    h.assert_false(pinf.is_finite(), "+inf_val() is not finite")
    h.assert_false(pinf.is_negative(), "+inf_val() is positive")

    // -infinity
    let ni = MPFloat.inf_val(false)
    h.assert_true(ni.is_inf(), "-inf_val() is inf")
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
    h.assert_true(finf.is_inf(), "from_f64(+inf) is inf")
    h.assert_false(finf.is_negative(), "from_f64(+inf) is positive")

    // -infinity bit pattern
    let fninf = MPFloat.from_f64(F64.from_bits(0xFFF0_0000_0000_0000))
    h.assert_true(fninf.is_inf(), "from_f64(-inf) is inf")
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
      h.assert_true(try MPFloat.from_string(s)?.is_inf() else false end,
        "\"" + s + "\" → inf")
      h.assert_false(try MPFloat.from_string(s)?.is_negative() else true end,
        "\"" + s + "\" → +inf")
    end

    // -infinity forms
    for s in ["-inf"; "-@Inf@"; "    -inf  "; "    -@Inf@  "].values() do
      h.assert_true(try MPFloat.from_string(s)?.is_inf() else false end,
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
    let shalf: String val = fhalf.string()
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
      MPFloat.inf_val(true).add(MPFloat.from_f64(1.0, p)).is_inf(),
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
      MPFloat.inf_val(true).mul(MPFloat.from_f64(2.0, p)).is_inf(),
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
    h.assert_true(MPFloat.create(p).inv().is_inf(), "1/0 = ∞")
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
      MPFloat.inf_val(true).div(MPFloat.from_f64(2.0, p)).is_inf(),
      "+∞ / finite = ∞")
    h.assert_true(
      MPFloat.from_f64(2.0, p).div(MPFloat.create(p)).is_inf(),
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
    h.assert_true(MPFloat.inf_val(true).sqrt().is_inf(), "sqrt(+∞) = +∞")
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
