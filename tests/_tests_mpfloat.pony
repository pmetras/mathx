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

use "random"
use "collections"
use "format"
use "debug"


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
    let z8 = MPFloat(64)
    h.assert_true(z8.string().at("0."), "MPFloat(64) starts with '0.'")


// ── Pi constant ────────────────────────────────────────────────────────────

class iso _DisabledTestMPFloatPi is UnitTest
  """
  MPFloat.pi(n) should converge to π = 3.14159…
  More bits of precision → more decimal digits in the output.
  """
  fun name(): String => "MPFloat/pi"

  fun apply(h: TestHelper) =>
    let pi4 = MPFloat.pi(32)
    h.assert_true(pi4.string().at("3."), "pi(32) starts with '3.'")

    let pi8 = MPFloat.pi(64)
    h.assert_true(pi8.string().at("3."), "pi(64) starts with '3.'")

    // Higher precision yields more output characters
    h.assert_true(
      pi8.string().size() > pi4.string().size(),
      "pi(64) has more characters than pi(32)")

    // First few digits of π (after stripping '_' separators)
    // string() outputs full decimal digits; 64 bits → ~19 digits
    // Expected prefix: "3.1415926"
    let s = _MPFStr(pi8.string())
    h.assert_true(s.at("3.1415926"), "pi(64) decimal prefix matches π")


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
    let zero = MPFloat(64)
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
    let zero = MPFloat(64)
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

    // All-zero digits → "0.0_000_000_0" for size 8 (64 bits)
    // string() returns String iso^; annotate as String val so it can be passed
    // as String box to _MPFStr.
    let z: String = MPFloat(64).string()
    h.assert_true(z.at("0."), "zero(64) starts with '0.'")
    // All decimal digits should be '0'
    let zs = _MPFStr(z)   // remove '_'  (String box accepted)
    h.assert_true(zs == "0.00000000", "zero(64) all decimal digits are 0")

    // Underscore placement: for size=4 (32 bits) the format is "<d>.<f0>_<f1><f2><f3>_"
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

    let z64 = MPFloat(64)
    h.assert_true(z64.is_zero(), "MPFloat(64) is zero")
    h.assert_true(z64.is_finite(), "MPFloat(64) is finite")
    h.assert_false(z64.is_negative(), "MPFloat(64) is not negative")


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


// ── from_f32 ─────────────────────────────────────────────────────────────────

class iso _TestMPFloatFromF32 is UnitTest
  """
  `from_f32` preserves sign, special values, and approximate magnitude.
  Digit values are verified indirectly through `string()`.
  """
  fun name(): String => "MPFloat/from_f32"

  fun apply(h: TestHelper) =>
    let p: ULong = 32
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // NaN propagation: quiet NaN bit pattern
    let fnan = MPFloat.from_f32(F32.from_bits(0x7FF8_0000))
    h.assert_true(fnan.is_nan(), "from_f32(NaN) is NaN")

    // See https://en.wikipedia.org/wiki/Single-precision_floating-point_format#Notable_single-precision_cases
    // +infinity bit pattern
    let finf = MPFloat.from_f32(F32.from_bits(0x7F80_0000))
    h.assert_true(finf.is_infinite(), "from_f32(+inf) is inf")
    h.assert_false(finf.is_negative(), "from_f32(+inf) is positive")

    // -infinity bit pattern
    let fninf = MPFloat.from_f32(F32.from_bits(0xFF80_0000))
    h.assert_true(fninf.is_infinite(), "from_f32(-inf) is inf")
    h.assert_true(fninf.is_negative(), "from_f32(-inf) is negative")

    // +0
    let fz = MPFloat.from_f32(0.0)
    h.assert_true(fz.is_zero(), "from_f32(0.0) is zero")
    h.assert_false(fz.is_negative(), "from_f32(0.0) is +0")

    // -0 detected via sign bit
    let fnz = MPFloat.from_f32(-0.0)
    h.assert_true(fnz.is_zero(), "from_f32(-0.0) is zero")
    h.assert_true(fnz.is_negative(), "from_f32(-0.0) is -0")

    // Positive value: string() must start with "3."
    let fpi = MPFloat.from_f32(3.14159, 64)
    h.assert_false(fpi.is_zero(), "from_f32(3.14) is not zero")
    h.assert_false(fpi.is_negative(), "from_f32(3.14) is positive")
    h.assert_true(fpi.is_finite(), "from_f32(3.14) is finite")
    h.assert_true(fpi.string().at("3."), "from_f32(3.14) string starts with \"3.\"")

    // Negative value
    let fneg = MPFloat.from_f32(-2.71828, 64)
    h.assert_true(fneg.is_negative(), "from_f32(-2.71) is negative")
    h.assert_true(fneg.string().at("-2."), "from_f32(-2.71) string starts with \"-2.\"")

    // Random values
    let rand = Rand
    for i in Range(0, 100) do
      let f = F32.from_bits(rand.u32())
      let mpf = MPFloat.from_f32(f, p)
      
      // 1. Round-trip check
      let rt_f = mpf.f32()
      if f.nan() then
        h.assert_true(rt_f.nan(), "Random floats [" + i.string() + "] NaN round-trip failed")
      elseif f.infinite() then
        h.assert_true(rt_f.infinite(), "Random floats [" + i.string() + "] Inf round-trip failed")
        h.assert_eq[Bool](f < 0, rt_f < 0, "Random floats [" + i.string() + "] Inf sign mismatch")
      else
        // Bit-accurate round-trip check
        h.assert_eq[F32](f, rt_f, "Random floats [" + i.string() + "] round-trip failed for " + f.string() + " got " + rt_f.string())
      end
      
      // 2. Semantic string representation check (only for finite parseable numbers)
      if f.finite() then
        try
          let s: String val = mpf.string()
          let f_parsed = s.f32()?
          let diff = (f - f_parsed).abs()
          let max_abs = f.abs().max(f_parsed.abs()).max(1.0)
          // Tighter tolerance: 2 ULPs for string parsing jitter
          h.assert_true(diff <= (max_abs * f.epsilon() * 2), 
            "Random floats [" + i.string() + "] string representation semantic mismatch: " + f.string() + " became " + s)
        else
          h.fail("Random floats [" + i.string() + "] MPFloat produced unparseable string")
        end
      end
    end


// ── from_f64 ─────────────────────────────────────────────────────────────────

class iso _TestMPFloatFromF64 is UnitTest
  """
  `from_f64` preserves sign, special values, and approximate magnitude.
  Digit values are verified indirectly through `string()`.
  """
  fun name(): String => "MPFloat/from_f64"

  fun apply(h: TestHelper) =>
    let p: ULong = 64
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // See https://en.wikipedia.org/wiki/Double-precision_floating-point_format#Double-precision_examples
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
    let fpi = MPFloat.from_f64(3.14159, 64)
    h.assert_false(fpi.is_zero(), "from_f64(3.14) is not zero")
    h.assert_false(fpi.is_negative(), "from_f64(3.14) is positive")
    h.assert_true(fpi.is_finite(), "from_f64(3.14) is finite")
    h.assert_true(fpi.string().at("3."), "from_f64(3.14) string starts with \"3.\"")

    // Negative value
    let fneg = MPFloat.from_f64(-2.71828, 64)
    h.assert_true(fneg.is_negative(), "from_f64(-2.71) is negative")
    h.assert_true(fneg.string().at("-2."), "from_f64(-2.71) string starts with \"-2.\"")

    // Random values
    let rand = Rand
    for i in Range(0, 100) do
      let f = F64.from_bits(rand.u64())
      let mpf = MPFloat.from_f64(f, p)
      
      // 1. Round-trip check
      let rt_f = mpf.f64()
      if f.nan() then
        h.assert_true(rt_f.nan(), "Random floats [" + i.string() + "] NaN round-trip failed")
      elseif f.infinite() then
        h.assert_true(rt_f.infinite(), "Random floats [" + i.string() + "] Inf round-trip failed")
        h.assert_eq[Bool](f < 0, rt_f < 0, "Random floats [" + i.string() + "] Inf sign mismatch")
      else
        // Bit-accurate round-trip check
        h.assert_eq[F64](f, rt_f, "Random floats [" + i.string() + "] round-trip failed for " + f.string() + " got " + rt_f.string())
      end
      
      // 2. Semantic string representation check (only for finite parseable numbers)
      if f.finite() then
        try
          let s: String val = mpf.string()
          let f_parsed = s.f64()?
          let diff = (f - f_parsed).abs()
          let max_abs = f.abs().max(f_parsed.abs()).max(1.0)
          // Tighter tolerance: 2 ULPs for string parsing jitter
          h.assert_true(diff <= (max_abs * f.epsilon() * 2), 
            "Random floats [" + i.string() + "] string representation semantic mismatch: " + f.string() + " became " + s)
        else
          h.fail("Random floats [" + i.string() + "] MPFloat produced unparseable string")
        end
      end
    end


// ── from_ulong ─────────────────────────────────────────────────────────────────

class iso _TestMPFloatFromULong is UnitTest
  """
  `from_ulong` preserves sign.
  Digit values are verified indirectly through `string()`.
  """
  fun name(): String => "MPFloat/from_ulong"

  fun apply(h: TestHelper) =>
    // 0
    let fz = MPFloat.from_ulong(0)
    h.assert_true(fz.is_zero(), "from_ulong(0) is zero")
    h.assert_true(fz == MPFloat.from_f32(0.0), "MPFloat.from_ulong(0) == MPFloat.from_f32(0.0)")
    h.assert_eq[MPFloat](fz, MPFloat.from_f64(0.0), "MPFloat.from_ulong(0) == MPFloat.from_f64(0.0)")
    h.assert_true(fz.is_integer(), "from_ulong(0) is integer")

    // Random numbers
    let rand = Rand
    for i in Range(1, 100) do
      let r: ULong = rand.ulong()
      let fr = MPFloat.from_ulong(r)
      h.assert_eq[String](fr.string().trim(0, try fr.string().find(".0")?.usize() else 0 end), r.string(), "Random from_ulong(" + r.string() + ")")
      h.assert_true(fr.is_integer(), "Random from_ulong(" + fr.string() + ") is integer")
    end

    // Maximum value
    let fpm = MPFloat.from_ulong(ULong.max_value())
    h.assert_eq[String](fpm.string().trim(0, try fpm.string().find(".0")?.usize() else 0 end), ULong.max_value().string(), "Maximum ULong value")
    h.assert_true(fpm.is_integer(), "from_ulong(ULong.max_value()) is integer")


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
  (not via F64), so the entire `prec` bits of precision are used regardless
  of how many significant digits the string contains.
  """
  fun name(): String => "MPFloat/from_string/decimal"

  fun apply(h: TestHelper) =>
    let p: ULong = 112
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // Integer: "3"
    let s3 = try MPFloat.from_string("3", p)? else MPFloat.nan_val() end
    h.assert_true(s3.is_finite(), "\"3\" is finite")
    h.assert_false(s3.is_negative(), "\"3\" is positive")
    ae(s3, MPFloat.from_f64(3.0, p), "\"3\" ≈ 3")

    // Negative: "-2.5" (exact in base 256)
    let sm = try MPFloat.from_string("-2.5", p)? else MPFloat.nan_val() end
    h.assert_true(sm.is_negative(), "\"-2.5\" is negative")
    h.assert_true(sm.is_finite(), "\"-2.5\" is finite")
    ae(sm, MPFloat.from_f64(-2.5, p), "\"-2.5\" ≈ -2.5")

    // Exponent notation: "314e-2" = 3.14
    let se = try MPFloat.from_string("314e-2", p)? else MPFloat.nan_val() end
    h.assert_true(se.is_finite(), "\"314e-2\" is finite")
    h.assert_false(se.is_negative(), "\"314e-2\" is positive")
    ae(se, MPFloat.from_f64(3.14, p), "\"314e-2\" ≈ 3.14")

    // '@' as exponent separator (GMP/MPFR style for base ≤ 10).
    let sat = try MPFloat.from_string("314@-2", p)? else MPFloat.nan_val() end
    h.assert_true(sat.is_finite(), "\"314@-2\" is finite")
    ae(sat, MPFloat.from_f64(3.14, p), "\"314@-2\" ≈ 3.14")

    // Large exponent: "1E3" = 1000 → finite
    let s1k = try MPFloat.from_string("1E3", 64)? else MPFloat.nan_val() end
    h.assert_true(s1k.is_finite(), "\"1E3\" is finite")
    h.assert_false(s1k.is_negative(), "\"1E3\" is positive")

    // Large exponent: "-5.E3859045" → finite (large _exponent, prec digits).
    let sle = try MPFloat.from_string("-5.E3859045", 64)? else MPFloat.nan_val() end
    h.assert_true(sle.is_finite(), "\"-5.E3859045\" is finite")
    h.assert_true(sle.is_negative(), "\"-5.E3859045\" is negative")

    // High-precision parsing: 17 significant digits exceed F64 precision (≈15).
    // "10000000000000002" → the final '2' must survive into the MPFloat.
    // With F64 this would round to 10000000000000000.
    // Round-trip: from_string × from_f64 comparison via string prefix check.
    let shp = try MPFloat.from_string("10000000000000002", 64)? else MPFloat.nan_val() end
    h.assert_true(shp.is_finite(), "17-digit string is finite")
    h.assert_false(shp.is_negative(), "17-digit string is positive")
    h.assert_true(
      shp.string().at("1000000000000000"),
      "17-digit string starts with \"1000000000000000\"")

    // Unsupported base errors.
    h.assert_false(
      try MPFloat.from_string("ff", 64, 16)?; true else false end,
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
    h.assert_true(MPFloat(64).string() == "0.0", "MPFloat(64) → \"0.0\"")

    // Positive value ≈ 3: string starts with "3."
    let f3 = MPFloat.from_f64(3.0, 32)
    h.assert_true(f3.string().at("3."), "from_f64(3.0, 32) starts with \"3.\"")

    // Negative value ≈ -2: string starts with "-2."
    let fneg = MPFloat.from_f64(-2.0, 32)
    h.assert_true(fneg.string().at("-2."), "from_f64(-2.0, 32) starts with \"-2.\"")

    // Value in (0, 1): string starts with "0."
    let fhalf = MPFloat.from_f64(0.5, 32)
    h.assert_true(fhalf.string().at("0."), "from_f64(0.5, 32) starts with \"0.\"")

    // Verify decimal accuracy: 0.5 → first fractional digit is 5
    let shalf: String = fhalf.string()
    h.assert_true(
      try shalf(2)? == '5' else false end,
      "from_f64(0.5, 32): first fractional digit is 5")


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
    let fp = MPFloat.from_f64(3.14, 32)
    h.assert_false(fp.is_negative(), "3.14 is positive")
    h.assert_true(fp.neg().is_negative(), "neg(3.14) is negative")
    h.assert_false(fp.neg().neg().is_negative(), "double neg restores positive")

    let fn2 = MPFloat.from_f64(-2.0, 32)
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
    let p: ULong = 160
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

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
    ae(MPFloat.from_f64(3.0, p).add(MPFloat.create(p)), MPFloat.from_f64(3.0, p), "x + 0 = x")
    ae(MPFloat.create(p).add(MPFloat.from_f64(3.0, p)), MPFloat.from_f64(3.0, p), "0 + x = x")

    // Same-sign addition.
    ae(MPFloat.from_f64(3.0, p).add(MPFloat.from_f64(2.0, p)), MPFloat.from_f64(5.0, p), "3 + 2 = 5")
    ae(MPFloat.from_f64(-3.0, p).add(MPFloat.from_f64(-2.0, p)), MPFloat.from_f64(-5.0, p), "-3 + -2 = -5")

    // Opposite-sign addition (subtraction of magnitudes).
    ae(MPFloat.from_f64(3.0, p).add(MPFloat.from_f64(-2.0, p)), MPFloat.from_f64(1.0, p), "3 + (-2) = 1")
    ae(MPFloat.from_f64(2.0, p).add(MPFloat.from_f64(-3.0, p)), MPFloat.from_f64(-1.0, p), "2 + (-3) = -1")

    // Exact cancellation → zero.
    h.assert_true(
      MPFloat.from_f64(3.0, p).add(MPFloat.from_f64(-3.0, p)).is_zero(),
      "3 + (-3) = 0")

    // Carry across digit boundary: 255/256 + 1/256 = 1.
    ae(MPFloat.from_f64(255.0 / 256.0, p).add(MPFloat.from_f64(1.0 / 256.0, p)),
      MPFloat.from_f64(1.0, p), "255/256 + 1/256 = 1")

    // Exponent alignment: 256 + 1 = 257.
    ae(MPFloat.from_f64(256.0, p).add(MPFloat.from_f64(1.0, p)),
      MPFloat.from_f64(257.0, p), "256 + 1 = 257")


class iso _TestMPFloatSub is UnitTest
  """
  `sub()` delegates to `add(that.neg())`, so it is a thin wrapper. The tests
  verify correct sign handling for both orderings and exact cancellation.
  """
  fun name(): String => "MPFloat/sub"

  fun apply(h: TestHelper) =>
    let p: ULong = 136
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    ae(MPFloat.from_f64(5.0, p).sub(MPFloat.from_f64(2.0, p)), MPFloat.from_f64(3.0, p), "5 - 2 = 3")
    ae(MPFloat.from_f64(2.0, p).sub(MPFloat.from_f64(5.0, p)), MPFloat.from_f64(-3.0, p), "2 - 5 = -3")
    h.assert_true(
      MPFloat.from_f64(5.0, p).sub(MPFloat.from_f64(5.0, p)).is_zero(),
      "5 - 5 = 0")
    ae(MPFloat.from_f64(-5.0, p).sub(MPFloat.from_f64(2.0, p)), MPFloat.from_f64(-7.0, p), "-5 - 2 = -7")
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
    let p: ULong = 144
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

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

    // Exact values.
    ae(MPFloat.from_f64(3.0, p).mul(MPFloat.from_f64(2.0, p)), MPFloat.from_f64(6.0, p), "3 × 2 = 6")
    ae(MPFloat.from_f64(10.0, p).mul(MPFloat.from_f64(10.0, p)), MPFloat.from_f64(100.0, p), "10 × 10 = 100")
    // Exponent propagation: (1/256) × 256 = 1.
    ae(MPFloat.from_f64(1.0 / 256.0, p).mul(MPFloat.from_f64(256.0, p)),
      MPFloat.from_f64(1.0, p), "(1/256) × 256 = 1")

    // Random numbers
    let rand = Rand
    for i in Range(0, 100) do
      let f1 = rand.real()
      let f2 = rand.real()
      ae(MPFloat.from_f64(f1, p) * MPFloat.from_f64(f2, p), MPFloat.from_f64(f1 * f2),
        "Random mul from F64 " + f1.string() + " × " + f2.string())

      let l1 = rand.ulong()
      let l2 = rand.ulong()
      ae(MPFloat.from_ulong(l1, p) * MPFloat.from_ulong(l2, p), MPFloat.from_mpint(MPInt.from[ULong](l1) * MPInt.from[ULong](l2)),
        "Random mul from ULong " + l1.string() + " × " + l2.string())
    end

    // Large random integers multiplications
    var mp1 = MPInt.from[ULong](rand.ulong())
    var mp2 = MPInt.from[ULong](rand.ulong())
    for i in Range(0, 30) do
      mp1 = mp1 * MPInt.from[ULong](rand.ulong())
      mp2 = mp2 * MPInt.from[ULong](rand.ulong())
      ae(MPFloat.from_mpint(mp1, p) * MPFloat.from_mpint(mp2, p), MPFloat.from_mpint(mp1 * mp2), "Random large integers mul " + mp1.string() + " × " + mp2.string())

      let f1 = rand.real()
      let f2 = rand.real()
      let mf1 = MPFloat.from_mpint(mp1, p) * MPFloat.from_f64(f1, p)
      let mf2 = MPFloat.from_mpint(mp2, p) * MPFloat.from_f64(f2, p)
      ae(mf1 * mf2, MPFloat.from_mpint(mp1, p) * MPFloat.from_mpint(mp2, p) * MPFloat.from_f64(f1 * f2),
        "Random large floats mul " + mf1.string() + " × " + mf2.string())
    end


// ── Inversion ──────────────────────────────────────────────────────────────

class iso _TestMPFloatInv is UnitTest
  """
  `inv` computes 1/x via Newton's method with correct sign and exponent.
  """
  fun name(): String => "MPFloat/inv"

  fun apply(h: TestHelper) =>
    let p: ULong = 80
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

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

    // Numeric values.
    ae(MPFloat.from_f64(1.0, p).inv(), MPFloat.from_f64(1.0, p), "1/1 = 1")
    ae(MPFloat.from_f64(2.0, p).inv(), MPFloat.from_f64(0.5, p), "1/2 = 0.5")
    // Round-trip: (1/2) × 2 = 1.
    let x = MPFloat.from_f64(2.0, p)
    ae(x.inv().mul(x), MPFloat.from_f64(1.0, p), "(1/2) × 2 = 1")
    // 1/256 and 1/(1/256) = 256.
    ae(MPFloat.from_f64(256.0, p).inv(), MPFloat.from_f64(1.0 / 256.0, p), "1/256")
    ae(MPFloat.from_f64(1.0 / 256.0, p).inv(), MPFloat.from_f64(256.0, p), "1/(1/256) = 256")


// ── Division ───────────────────────────────────────────────────────────────

class iso _TestMPFloatDiv is UnitTest
  """
  `div` computes this/that = this × (1/that) with correct sign and exponent.
  """
  fun name(): String => "MPFloat/div"

  fun apply(h: TestHelper) =>
    let p: ULong = 104
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

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
    ae(MPFloat.from_f64(6.0, p).div(MPFloat.from_f64(2.0, p)), MPFloat.from_f64(3.0, p), "6 / 2 = 3")
    ae(MPFloat.from_f64(1.0, p).div(MPFloat.from_f64(4.0, p)), MPFloat.from_f64(0.25, p), "1 / 4 = 0.25")
    // Round-trip: (6 / 2) × 2 = 6.
    let xrt = MPFloat.from_f64(6.0, p)
    let yrt = MPFloat.from_f64(2.0, p)
    ae(xrt.div(yrt).mul(yrt), xrt, "(6/2) × 2 = 6")


class iso _TestMPFloatSqrt is UnitTest
  """
  `sqrt` computes √this with correct sign, exponent, and IEEE 754 special
  cases. Uses the parity-split Newton reciprocal-square-root algorithm.
  """
  fun name(): String => "MPFloat/sqrt"

  fun apply(h: TestHelper) =>
    let p: ULong = 96
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

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

    // Exact values.
    ae(MPFloat.from_f64(4.0, p).sqrt(), MPFloat.from_f64(2.0, p), "sqrt(4) = 2")
    ae(MPFloat.from_f64(16.0, p).sqrt(), MPFloat.from_f64(4.0, p), "sqrt(16) = 4")
    ae(MPFloat.from_f64(64.0, p).sqrt(), MPFloat.from_f64(8.0, p), "sqrt(64) = 8")
    // √0.25 = 0.5 (even exponent, exact in base 256).
    ae(MPFloat.from_f64(0.25, p).sqrt(), MPFloat.from_f64(0.5, p), "sqrt(0.25) = 0.5")
    // √2 ≈ 1.41421… (irrational).
    ae(MPFloat.from_f64(2.0, p).sqrt(), MPFloat.from_f64(F64(2.0).sqrt(), p), "sqrt(2) ≈ 1.41421…")
    // Round-trip: sqrt(4)² = 4.
    let s4 = MPFloat.from_f64(4.0, p).sqrt()
    ae(s4.mul(s4), MPFloat.from_f64(4.0, p), "sqrt(4)² = 4")


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
    let p: ULong = 120
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
    let p: ULong = 192
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // Zero: always +0.
    let z = MPFloat.from_mpint(MPInt.from[ILong](0), p)
    h.assert_true(z.is_zero(),     "from_mpint(0) is zero")
    h.assert_false(z.is_negative(), "from_mpint(0) is +0")
    h.assert_false(z.is_nan(),     "from_mpint(0) is not NaN")
    h.assert_true(z.is_finite(),   "from_mpint(0) is finite")

    // Positive small integer: 1.
    let one = MPFloat.from_mpint(MPInt.from[ILong](1), p)
    h.assert_false(one.is_negative(), "from_mpint(1) is positive")
    h.assert_true(one.is_finite(),    "from_mpint(1) is finite")
    ae(one, MPFloat.from_f64(1.0, p), "from_mpint(1) ≈ 1")

    // Negative integer: -3.
    let neg3 = MPFloat.from_mpint(MPInt.from[ILong](-3), p)
    h.assert_true(neg3.is_negative(), "from_mpint(-3) is negative")
    ae(neg3, MPFloat.from_f64(-3.0, p), "from_mpint(-3) ≈ -3")

    // Larger integer: 1000.
    let thou = MPFloat.from_mpint(MPInt.from[ILong](1000), p)
    h.assert_false(thou.is_negative(), "from_mpint(1000) is positive")
    ae(thou, MPFloat.from_f64(1000.0, p), "from_mpint(1000) ≈ 1000")

    // Power-of-2: 256 and 65536.
    ae(MPFloat.from_mpint(MPInt.from[ILong](256), p), MPFloat.from_f64(256.0, p), "from_mpint(256) = 256")
    ae(MPFloat.from_mpint(MPInt.from[ILong](65536), p), MPFloat.from_f64(65536.0, p), "from_mpint(65536) = 65536")

    // Negative large: -65536.
    let ni65536 = MPFloat.from_mpint(MPInt.from[ILong](-65536), p)
    h.assert_true(ni65536.is_negative(), "from_mpint(-65536) is negative")
    ae(ni65536, MPFloat.from_f64(-65536.0, p), "from_mpint(-65536) ≈ -65536")

    // Round-trip sign: from_mpint(n).is_negative() == n.is_negative()
    let npos = MPInt.from[ILong](42)
    let nneg = MPInt.from[ILong](-42)
    h.assert_false(MPFloat.from_mpint(npos, p).is_negative(), "positive MPInt → positive MPFloat")
    h.assert_true( MPFloat.from_mpint(nneg, p).is_negative(), "negative MPInt → negative MPFloat")

    // Precision parameter is honoured: value must be finite and positive.
    let big = MPFloat.from_mpint(MPInt.from[ILong](1000000), 32)
    h.assert_true(big.is_finite(),  "from_mpint(1e6, prec=32) is finite")
    h.assert_false(big.is_negative(), "from_mpint(1e6, prec=32) is positive")
    ae(big, MPFloat.from_f64(1000000.0, 32), "from_mpint(1e6, prec=32) ≈ 1e6")

    // ── Large MPInt (beyond I64 / F64 range) ─────────────────────────────────

    // 10^20 > U64.max (≈1.8×10^19): requires MPInt arithmetic to construct.
    //
    // With p=64 bits (≈19.3 significant decimal digits), 10^20 (21 digits) cannot
    // be held exactly; the best 8-byte approximation is just below 10^20 and
    // its string representation starts with "9".  We only verify finiteness,
    // positiveness, and that the leading digit is "9" or "1" (order-of-
    // magnitude correct).
    let e20: MPInt = try MPInt.from_string("100000000000000000000")? else MPInt.from[ILong](0) end
    let fe20 = MPFloat.from_mpint(e20, p)
    h.assert_true(fe20.is_finite(),    "from_mpint(10^20, p=64 bits) is finite")
    h.assert_false(fe20.is_negative(), "from_mpint(10^20, p=64 bits) is positive")
    h.assert_true(
      fe20.string().at("9") or fe20.string().at("1"),
      "from_mpint(10^20, p=64 bits) leading digit is 9 or 1")

    // With p=96 bits (≈29 significant decimal digits), 10^20 (21 digits) fits
    // comfortably and the string should represent "1e+20".
    let fe20p12 = MPFloat.from_mpint(e20, 96)
    h.assert_true(fe20p12.is_finite(),    "from_mpint(10^20, p=96) is finite")
    h.assert_false(fe20p12.is_negative(), "from_mpint(10^20, p=96) is positive")
    h.assert_true(
      fe20p12.string().at("1e+20") or fe20p12.string().at("100000000000000000000"),
      "from_mpint(10^20, p=96) string represents 10^20")

    // 29-digit positive: precision must be ≥ ceil(29 / log10(2)) ≈ 96 bits
    // to represent all digits. With p=128 bits all 29 digits survive.
    let d29: MPInt =
      try MPInt.from_string("12345678901234567890123456789")?
      else MPInt.from[ILong](0) end
    let fd29 = MPFloat.from_mpint(d29, 128)
    h.assert_true(fd29.is_finite(),    "from_mpint(29-digit) is finite")
    h.assert_false(fd29.is_negative(), "from_mpint(29-digit) is positive")
    h.assert_true(
      fd29.string().at("12345678901234567890"),
      "from_mpint(29-digit) first 20 digits preserved with prec=128")

    // Negative large: -10^20 with p=96 for sufficient precision.
    let ne20: MPInt =
      try MPInt.from_string("-100000000000000000000")? else MPInt.from[ILong](0) end
    let fne20 = MPFloat.from_mpint(ne20, 96)
    h.assert_true(fne20.is_finite(),    "from_mpint(-10^20) is finite")
    h.assert_true(fne20.is_negative(),  "from_mpint(-10^20) is negative")
    h.assert_true(
      fne20.string().at("-1e+20") or fne20.string().at("-100000000000000000000"),
      "from_mpint(-10^20, p=96) string represents -10^20")

    // 39-digit positive with precision of 400 bits
    let d39: MPInt =
      try MPInt.from_string("123456789012345678901234567890123456789")?
      else MPInt.from[ILong](0) end
    let fd39 = MPFloat.from_mpint(d39, 400)
    h.assert_true(fd39.is_finite(),    "from_mpint(39-digit) is finite")
    h.assert_false(fd39.is_negative(), "from_mpint(39-digit) is positive")
    h.assert_true(
      fd39.string().at("12345678901234567890"),
      "from_mpint(39-digit) first 20 digits preserved with prec=400")
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
    let p: ULong = 240
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // Helpers: build MPFloat from a decimal string.
    let f = {(s: String): MPFloat =>
      try MPFloat.from_string(s, p)? else MPFloat.create(p) end
    }

    // ── Exact integer operands ────────────────────────────────────────────────

    // 7 / 2 = 3 rem 1
    (let q7_2, let r7_2) = f("7").divrem(f("2"))
    ae(q7_2, f("3"), "divrem(7,2) q = 3")
    ae(r7_2, f("1"), "divrem(7,2) r = 1")
    h.assert_false(r7_2.is_negative(), "divrem(7,2) r ≥ 0")

    // -7 / 2 = -3 rem -1  (truncation toward zero)
    (let qn7_2, let rn7_2) = f("-7").divrem(f("2"))
    ae(qn7_2, f("-3"), "divrem(-7,2) q = -3")
    ae(rn7_2, f("-1"), "divrem(-7,2) r = -1")
    h.assert_true(rn7_2.is_negative(), "divrem(-7,2) r ≤ 0")

    // 7 / -2 = -3 rem 1  (truncation toward zero)
    (let q7_n2, let r7_n2) = f("7").divrem(f("-2"))
    ae(q7_n2, f("-3"), "divrem(7,-2) q = -3")
    ae(r7_n2, f("1"),  "divrem(7,-2) r = 1")
    h.assert_false(r7_n2.is_negative(), "divrem(7,-2) r ≥ 0")

    // -7 / -2 = 3 rem -1
    (let qn7_n2, let rn7_n2) = f("-7").divrem(f("-2"))
    ae(qn7_n2, f("3"),  "divrem(-7,-2) q = 3")
    ae(rn7_n2, f("-1"), "divrem(-7,-2) r = -1")

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
    ae(q10_5, f("2"), "divrem(10,5) q = 2")
    h.assert_true(r10_5.is_zero(), "divrem(10,5) r = 0")

    // ── rem ───────────────────────────────────────────────────────────────────

    // rem sign follows dividend
    h.assert_false(f("7").rem(f("2")).is_negative(),   "rem(7,2) ≥ 0")
    h.assert_true( f("-7").rem(f("2")).is_negative(),  "rem(-7,2) ≤ 0")
    h.assert_false(f("7").rem(f("-2")).is_negative(),  "rem(7,-2) ≥ 0")
    h.assert_true( f("-7").rem(f("-2")).is_negative(), "rem(-7,-2) ≤ 0")

    // ── fld ──────────────────────────────────────────────────────────────────

    // Same-sign operands: fld = trunc
    ae(f("7").fld(f("2")),   f("3"),  "fld(7,2) = 3")
    ae(f("-7").fld(f("-2")), f("3"),  "fld(-7,-2) = 3")

    // Opposite-sign: fld = trunc - 1
    ae(f("-7").fld(f("2")),  f("-4"), "fld(-7,2) = -4")
    ae(f("7").fld(f("-2")),  f("-4"), "fld(7,-2) = -4")

    // Exact division: no adjustment needed.
    ae(f("6").fld(f("2")),   f("3"),  "fld(6,2) = 3")
    ae(f("-6").fld(f("2")),  f("-3"), "fld(-6,2) = -3")
    // Non-dyadic exact division (10/5 = 2; exercises post-correction via divrem).
    ae(f("10").fld(f("5")),  f("2"),  "fld(10,5) = 2")
    h.assert_true(f("10").mod(f("5")).is_zero(), "mod(10,5) = 0")

    // ── mod ──────────────────────────────────────────────────────────────────

    // mod sign follows divisor
    h.assert_false(f("7").mod(f("2")).is_negative(),   "mod(7,2) ≥ 0")
    h.assert_false(f("-7").mod(f("2")).is_negative(),  "mod(-7,2) ≥ 0  (sign = divisor)")
    h.assert_true( f("7").mod(f("-2")).is_negative(),  "mod(7,-2) ≤ 0  (sign = divisor)")
    h.assert_true( f("-7").mod(f("-2")).is_negative(), "mod(-7,-2) ≤ 0")

    // mod values: mod(-7, 2) = 1,  mod(7, -2) = -1
    ae(f("-7").mod(f("2")),  f("1"),  "mod(-7,2) = 1")
    ae(f("7").mod(f("-2")),  f("-1"), "mod(7,-2) = -1")

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
    ae(qu, f("3"), "divrem_unsafe(7,2) q = 3")
    ae(ru, f("1"), "divrem_unsafe(7,2) r = 1")

    ae(f("7").rem_unsafe(f("2")),   f("1"),  "rem_unsafe(7,2) = 1")
    ae(f("-7").rem_unsafe(f("2")),  f("-1"), "rem_unsafe(-7,2) = -1")
    ae(f("-7").fld_unsafe(f("2")),  f("-4"), "fld_unsafe(-7,2) = -4")
    ae(f("-7").mod_unsafe(f("2")),  f("1"),  "mod_unsafe(-7,2) = 1")


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
    let p: ULong = 72
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // from_string uses exact MPInt arithmetic for exact binary fractions (e.g.
    // "2.0", "0.5", "0.25") and best-rounded approximation for others ("0.1").
    let f = {(s: String): MPFloat =>
      try MPFloat.from_string(s, p)? else MPFloat.create(p) end
    }
    let fi = {(n: F64): MPFloat => MPFloat.from_f64(n, p) }

    // ── trunc ────────────────────────────────────────────────────────────────

    ae(f("2.7").trunc(),   fi(2.0),   "trunc(2.7) = 2")
    ae(f("0.9").trunc(),   fi(0.0),   "trunc(0.9) = 0")
    ae(f("2.0").trunc(),   fi(2.0),   "trunc(2.0) = 2")
    ae(f("256.0").trunc(), fi(256.0), "trunc(256.0) = 256")
    ae(f("-2.7").trunc(),  fi(-2.0),  "trunc(-2.7) = -2")
    h.assert_true(f("-0.9").trunc().is_zero(), "trunc(-0.9) = 0")
    ae(f("-2.0").trunc(),  fi(-2.0),  "trunc(-2.0) = -2")

    // Special values pass through.
    h.assert_true(MPFloat.nan_val().trunc().is_nan(),        "trunc(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().trunc().is_infinite(),   "trunc(+inf) = +inf")
    h.assert_false(MPFloat.inf_val().trunc().is_negative(),  "trunc(+inf) is positive")

    // ── floor ────────────────────────────────────────────────────────────────

    ae(f("2.7").floor(),   fi(2.0),  "floor(2.7) = 2")
    ae(f("2.0").floor(),   fi(2.0),  "floor(2.0) = 2")
    ae(f("0.9").floor(),   fi(0.0),  "floor(0.9) = 0")
    ae(f("-2.7").floor(),  fi(-3.0), "floor(-2.7) = -3")
    ae(f("-2.0").floor(),  fi(-2.0), "floor(-2.0) = -2")
    ae(f("-0.1").floor(),  fi(-1.0), "floor(-0.1) = -1")

    h.assert_true(MPFloat.nan_val().floor().is_nan(),            "floor(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().floor().is_infinite(),       "floor(+inf) = +inf")
    h.assert_true(MPFloat.inf_val().neg().floor().is_infinite(), "floor(-inf) = -inf")

    // ── ceil ─────────────────────────────────────────────────────────────────

    ae(f("2.1").ceil(),   fi(3.0),  "ceil(2.1) = 3")
    ae(f("2.0").ceil(),   fi(2.0),  "ceil(2.0) = 2")
    ae(f("0.1").ceil(),   fi(1.0),  "ceil(0.1) = 1")
    ae(f("-2.7").ceil(),  fi(-2.0), "ceil(-2.7) = -2")
    ae(f("-2.0").ceil(),  fi(-2.0), "ceil(-2.0) = -2")
    h.assert_true(f("-0.9").ceil().is_zero(), "ceil(-0.9) = 0")

    h.assert_true(MPFloat.nan_val().ceil().is_nan(),       "ceil(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().ceil().is_infinite(),  "ceil(+inf) = +inf")

    // ── round (half away from zero) ───────────────────────────────────────────

    ae(f("2.4").round(),  fi(2.0),  "round(2.4) = 2")
    ae(f("2.6").round(),  fi(3.0),  "round(2.6) = 3")
    ae(f("0.4").round(),  fi(0.0),  "round(0.4) = 0")
    ae(f("0.5").round(),  fi(1.0),  "round(0.5) = 1  (half away)")
    ae(f("2.5").round(),  fi(3.0),  "round(2.5) = 3  (half away)")
    ae(f("-2.4").round(), fi(-2.0), "round(-2.4) = -2")
    ae(f("-2.6").round(), fi(-3.0), "round(-2.6) = -3")
    ae(f("-0.4").round(), fi(0.0),  "round(-0.4) = 0")
    ae(f("-0.5").round(), fi(-1.0), "round(-0.5) = -1  (half away)")
    ae(f("-2.5").round(), fi(-3.0), "round(-2.5) = -3  (half away)")
    ae(f("3").round(),    fi(3.0),  "round(3) = 3")
    ae(f("-3").round(),   fi(-3.0), "round(-3) = -3")

    h.assert_true(MPFloat.nan_val().round().is_nan(),       "round(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().round().is_infinite(),  "round(+inf) = +inf")


// ── is_integer ────────────────────────────────────────────────────────────────

class iso _TestMPFloatIsInteger is UnitTest
  """
  `is_integer()` returns true iff the value is finite and has no non-zero
  fractional base-256 bytes.
  """

  fun name(): String => "MPFloat/is_integer"

  fun apply(h: TestHelper) =>
    let p: ULong = 168

    // Finite exact integers.
    h.assert_true(MPFloat.from_f64(0.0,   p).is_integer(), "0 is integer")
    h.assert_true(MPFloat.from_f64(1.0,   p).is_integer(), "1 is integer")
    h.assert_true(MPFloat.from_f64(-3.0,  p).is_integer(), "-3 is integer")
    h.assert_true(MPFloat.from_f64(256.0, p).is_integer(), "256 is integer")
    h.assert_true(MPFloat.from_f64(1024.0, p).is_integer(), "1024 is integer")

    // Finite values with a fractional part.
    h.assert_false(MPFloat.from_f64(0.5,  p).is_integer(), "0.5 not integer")
    h.assert_false(MPFloat.from_f64(1.5,  p).is_integer(), "1.5 not integer")
    h.assert_false(MPFloat.from_f64(-2.25, p).is_integer(), "-2.25 not integer")

    // Special values.
    h.assert_false(MPFloat.nan_val().is_integer(),  "NaN not integer")
    h.assert_false(MPFloat.inf_val().is_integer(),  "+inf not integer")
    h.assert_false(MPFloat.inf_val().neg().is_integer(), "-inf not integer")


// ── ln ────────────────────────────────────────────────────────────────────────

class iso _TestMPFloatLn is UnitTest
  """
  Natural logarithm `ln(x)`:
  - Special-case propagation (NaN, ±∞, non-positive).
  - Known values compared with `almost_eq` (rel_tol = 1e-14).
  - Round-trip: exp(ln(x)) ≈ x within tolerance.
  """

  fun name(): String => "MPFloat/ln"

  fun apply(h: TestHelper) =>
    let p: ULong = 160
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()

    // Helper: almost-equal assertion.
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // Special cases.
    h.assert_true(MPFloat.nan_val().ln().is_nan(),          "ln(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().ln().is_infinite(),     "ln(+inf) = +inf")
    h.assert_false(MPFloat.inf_val().ln().is_negative(),    "ln(+inf) > 0")
    h.assert_true(MPFloat.from_f64(-1.0, p).ln().is_nan(), "ln(-1) = NaN")
    h.assert_true(MPFloat.create(p).ln().is_infinite(),     "ln(0) = -inf")
    h.assert_true(MPFloat.create(p).ln().is_negative(),     "ln(0) is -inf (negative)")

    // ln(1) = 0 (exact).
    h.assert_true(MPFloat.from_f64(1.0, p).ln().is_zero(), "ln(1) = 0")

    // ln(2) ≈ 0.6931471805599453.
    ae(MPFloat.from_f64(2.0, p).ln(),
       MPFloat.from_f64(0.6931471805599453, p),
       "ln(2) ≈ 0.6931471805599453")

    // ln(8) = 3 × ln(2) ≈ 2.0794415416798357.
    ae(MPFloat.from_f64(8.0, p).ln(),
       MPFloat.from_f64(2.0794415416798357, p),
       "ln(8) ≈ 2.079441541679836")

    // Round-trip: exp(ln(x)) ≈ x for x = 2, 7.
    let x2: MPFloat = MPFloat.from_f64(2.0, p)
    ae(x2.ln().exp(), x2, "exp(ln(2)) ≈ 2")

    let x7: MPFloat = MPFloat.from_f64(7.0, p)
    ae(x7.ln().exp(), x7, "exp(ln(7)) ≈ 7")


// ── exp ───────────────────────────────────────────────────────────────────────

class iso _TestMPFloatExp is UnitTest
  """
  Natural exponential `e^x`:
  - Special-case propagation.
  - Known values compared with `almost_eq` (rel_tol = 1e-14).
  - Round-trip: ln(exp(x)) ≈ x within tolerance.
  """

  fun name(): String => "MPFloat/exp"

  fun apply(h: TestHelper) =>
    let p: ULong = 112
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()

    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // Special cases.
    h.assert_true(MPFloat.nan_val().exp().is_nan(),         "exp(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().exp().is_infinite(),    "exp(+inf) = +inf")
    h.assert_false(MPFloat.inf_val().exp().is_negative(),   "exp(+inf) > 0")
    h.assert_true(MPFloat.inf_val().neg().exp().is_zero(),  "exp(-inf) = 0")

    // exp(0) = 1 (exact).
    h.assert_true(MPFloat.create(p).exp().almost_eq(
      MPFloat.from_f64(1.0, p), tol, tol), "exp(0) = 1")

    // exp(1) ≈ 2.718281828459045235.
    ae(MPFloat.from_f64(1.0, p).exp(),
       MPFloat.from_f64(2.718281828459045, p),
       "exp(1) ≈ 2.718281828459045")

    // exp(-1) ≈ 0.36787944117144233.
    ae(MPFloat.from_f64(-1.0, p).exp(),
       MPFloat.from_f64(0.36787944117144233, p),
       "exp(-1) ≈ 0.36787944117144233")

    // exp(5) ≈ 148.4131591025766.
    ae(MPFloat.from_f64(5.0, p).exp(),
       MPFloat.from_f64(148.4131591025766, p),
       "exp(5) ≈ 148.4131591025766")

    // Round-trip: ln(exp(x)) ≈ x for x = 3.
    let x3: MPFloat = MPFloat.from_f64(3.0, p)
    ae(x3.exp().ln(), x3, "ln(exp(3)) ≈ 3")


// ── log2 ──────────────────────────────────────────────────────────────────────

class iso _TestMPFloatLog2 is UnitTest
  """
  Base-2 logarithm `log₂(x)`:
  - Special cases from ln.
  - Values compared with `almost_eq` (rel_tol = 1e-14).
  """

  fun name(): String => "MPFloat/log2"

  fun apply(h: TestHelper) =>
    let p: ULong = 144
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()

    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // Special cases.
    h.assert_true(MPFloat.nan_val().log2().is_nan(),          "log2(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().log2().is_infinite(),     "log2(+inf) = +inf")
    h.assert_true(MPFloat.from_f64(-1.0, p).log2().is_nan(), "log2(-1) = NaN")
    h.assert_true(MPFloat.create(p).log2().is_infinite(),     "log2(0) = -inf")
    h.assert_true(MPFloat.create(p).log2().is_negative(),     "log2(0) negative")

    // log2(1) = 0 (exact).
    h.assert_true(MPFloat.from_f64(1.0, p).log2().is_zero(), "log2(1) = 0")

    // log2(2) ≈ 1.
    ae(MPFloat.from_f64(2.0, p).log2(),
       MPFloat.from_f64(1.0, p),
       "log2(2) ≈ 1")

    // log2(8) ≈ 3.
    ae(MPFloat.from_f64(8.0, p).log2(),
       MPFloat.from_f64(3.0, p),
       "log2(8) ≈ 3")

    // log2(0.5) ≈ -1.
    ae(MPFloat.from_f64(0.5, p).log2(),
       MPFloat.from_f64(-1.0, p),
       "log2(0.5) ≈ -1")


// ── log10 ─────────────────────────────────────────────────────────────────────

class iso _TestMPFloatLog10 is UnitTest
  """
  Base-10 logarithm `log₁₀(x)`:
  - Special cases.
  - Values compared with `almost_eq` (rel_tol = 1e-14).
  """

  fun name(): String => "MPFloat/log10"

  fun apply(h: TestHelper) =>
    let p: ULong = 200
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()

    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // Special cases.
    h.assert_true(MPFloat.nan_val().log10().is_nan(),          "log10(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().log10().is_infinite(),     "log10(+inf) = +inf")
    h.assert_true(MPFloat.from_f64(-1.0, p).log10().is_nan(), "log10(-1) = NaN")
    h.assert_true(MPFloat.create(p).log10().is_infinite(),     "log10(0) = -inf")
    h.assert_true(MPFloat.create(p).log10().is_negative(),     "log10(0) negative")

    // log10(1) = 0 (exact).
    h.assert_true(MPFloat.from_f64(1.0, p).log10().is_zero(), "log10(1) = 0")

    // log10(10) ≈ 1.
    ae(MPFloat.from_f64(10.0,   p).log10(), MPFloat.from_f64(1.0, p), "log10(10) ≈ 1")

    // log10(100) ≈ 2.
    ae(MPFloat.from_f64(100.0,  p).log10(), MPFloat.from_f64(2.0, p), "log10(100) ≈ 2")

    // log10(1000) ≈ 3.
    ae(MPFloat.from_f64(1000.0, p).log10(), MPFloat.from_f64(3.0, p), "log10(1000) ≈ 3")

    // log10(0.1) ≈ -1.  Note: from_f64(0.1) is the nearest F64 to 1/10.
    ae(MPFloat.from_f64(0.1, p).log10(), MPFloat.from_f64(-1.0, p), "log10(0.1) ≈ -1")


// ── exp2 ──────────────────────────────────────────────────────────────────────

class iso _TestMPFloatExp2 is UnitTest
  """
  Base-2 exponential `2^x`:
  - Special cases.
  - Values compared with `almost_eq` (rel_tol = 1e-14).
  """

  fun name(): String => "MPFloat/exp2"

  fun apply(h: TestHelper) =>
    let p: ULong = 160
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()

    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    // Special cases.
    h.assert_true(MPFloat.nan_val().exp2().is_nan(),          "exp2(NaN) = NaN")
    h.assert_true(MPFloat.inf_val().exp2().is_infinite(),     "exp2(+inf) = +inf")
    h.assert_true(MPFloat.inf_val().neg().exp2().is_zero(),   "exp2(-inf) = 0")

    // exp2(0) = 1 (exact).
    ae(MPFloat.create(p).exp2(),            MPFloat.from_f64(1.0,    p), "exp2(0) = 1")

    // exp2(1) = 2.
    ae(MPFloat.from_f64(1.0,  p).exp2(),   MPFloat.from_f64(2.0,    p), "exp2(1) = 2")

    // exp2(3) = 8.
    ae(MPFloat.from_f64(3.0,  p).exp2(),   MPFloat.from_f64(8.0,    p), "exp2(3) = 8")

    // exp2(-1) = 0.5 (exact in base 256).
    ae(MPFloat.from_f64(-1.0, p).exp2(),   MPFloat.from_f64(0.5,    p), "exp2(-1) = 0.5")

    // exp2(10) = 1024.
    ae(MPFloat.from_f64(10.0, p).exp2(),   MPFloat.from_f64(1024.0, p), "exp2(10) = 1024")


// ── powi ──────────────────────────────────────────────────────────────────────

class iso _TestMPFloatPowi is UnitTest
  """
  Integer power `x^n` via binary exponentiation:
  - x^0 = 1 for any finite x.
  - Results compared with `almost_eq` (rel_tol = 1e-14).
  - Negative exponents via inv().
  - NaN and ±∞ propagation.
  """

  fun name(): String => "MPFloat/powi"

  fun apply(h: TestHelper) =>
    let p: ULong = 152
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()

    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    let two   = MPFloat.from_f64(2.0,  p)
    let three = MPFloat.from_f64(3.0,  p)

    // x^0 = 1.
    ae(two.powi(0),   MPFloat.from_f64(1.0,    p), "2^0 = 1")
    ae(three.powi(0), MPFloat.from_f64(1.0,    p), "3^0 = 1")

    // Positive powers.
    ae(two.powi(1),   MPFloat.from_f64(2.0,    p), "2^1 = 2")
    ae(two.powi(3),   MPFloat.from_f64(8.0,    p), "2^3 = 8")
    ae(two.powi(10),  MPFloat.from_f64(1024.0, p), "2^10 = 1024")
    ae(three.powi(3), MPFloat.from_f64(27.0,   p), "3^3 = 27")

    // Negative powers (computed via inv() which uses Newton iteration).
    ae(two.powi(-1),  MPFloat.from_f64(0.5,    p), "2^(-1) = 0.5")
    ae(two.powi(-2),  MPFloat.from_f64(0.25,   p), "2^(-2) = 0.25")

    // NaN and ±∞.
    h.assert_true(MPFloat.nan_val().powi(2).is_nan(),       "NaN^2 = NaN")
    h.assert_true(MPFloat.inf_val().powi(2).is_infinite(),  "+inf^2 = +inf")


// ── pow ───────────────────────────────────────────────────────────────────────

class iso _TestMPFloatPow is UnitTest
  """
  General real power `x^y = exp(y × ln(x))`:
  - x > 0: computed via exp/ln; results compared with `almost_eq`.
  - x = 0: 0^pos = 0, 0^neg = NaN.
  - x < 0 with integer y: delegates to powi (exact).
  - x < 0 with non-integer y: NaN.
  - x^0 = 1 for any finite x.
  """

  fun name(): String => "MPFloat/pow"

  fun apply(h: TestHelper) =>
    let p: ULong = 136
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()

    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    let two     = MPFloat.from_f64(2.0,  p)
    let three   = MPFloat.from_f64(3.0,  p)
    let four    = MPFloat.from_f64(4.0,  p)
    let zero    = MPFloat.create(p)
    let half    = MPFloat.from_f64(0.5,  p)
    let neg_one = MPFloat.from_f64(-1.0, p)
    let neg_two = MPFloat.from_f64(-2.0, p)

    // x^0 = 1.
    ae(two.pow(zero),  MPFloat.from_f64(1.0, p), "2^0 = 1")
    ae(zero.pow(zero), MPFloat.from_f64(1.0, p), "0^0 = 1 (convention)")

    // Positive base, integer exponents: exp(n×ln(x)).
    ae(two.pow(three),   MPFloat.from_f64(8.0, p), "2^3 = 8")
    ae(three.pow(two),   MPFloat.from_f64(9.0, p), "3^2 = 9")

    // Positive base, fractional exponent: 4^0.5 = sqrt(4) = 2.
    ae(four.pow(half),   MPFloat.from_f64(2.0, p), "4^0.5 = 2")

    // 0^positive = 0.
    h.assert_true(zero.pow(two).is_zero(), "0^2 = 0")

    // 0^negative = NaN.
    h.assert_true(zero.pow(neg_one).is_nan(), "0^(-1) = NaN")

    // Negative base, integer exponent (delegates to powi — exact).
    ae(neg_two.pow(three), MPFloat.from_f64(-8.0, p), "(-2)^3 = -8")
    ae(neg_two.pow(two),   MPFloat.from_f64(4.0,  p), "(-2)^2 = 4")

    // Negative base, non-integer exponent = NaN.
    h.assert_true(neg_two.pow(half).is_nan(), "(-2)^0.5 = NaN")

    // NaN propagation.
    h.assert_true(MPFloat.nan_val().pow(two).is_nan(), "NaN^2 = NaN")
    h.assert_true(two.pow(MPFloat.nan_val()).is_nan(),  "2^NaN = NaN")


// ── High-precision ────────────────────────────────────────────────────────────

class iso _TestMPFloatHighPrec is UnitTest
  """
  Arithmetic and transcendental functions at 100-byte (~240-digit) precision.

  Tests both |x| > 1 (large values) and |x| < 1 (small / fractional values).
  The tolerance 1e-100 checks at least 100 decimal digits of accuracy,
  which is far less than the theoretical maximum (~240 digits) but exercises
  convergence of all Newton / Taylor series at high precision.
  """

  fun name(): String => "MPFloat/high_precision"

  fun apply(h: TestHelper) =>
    let p: ULong = 800
    // 100-byte precision → ~240 decimal digits.  We demand 100-digit accuracy.
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    let one  = MPFloat.from_f64(1.0,  p)
    let two  = MPFloat.from_f64(2.0,  p)
    let four = MPFloat.from_f64(4.0,  p)
    let half = MPFloat.from_f64(0.5,  p)
    let qtr  = MPFloat.from_f64(0.25, p)

    // ── Basic arithmetic — large values (|x| > 1) ─────────────────────────

    // 2 × 0.5 = 1  (exact)
    ae(two.mul(half), one, "2 × 0.5 = 1 at p=100")
    // 1 / 2 = 0.5  (exact)
    ae(two.inv(), half, "1/2 = 0.5 at p=100")
    // 4 / 2 = 2    (exact)
    ae(four.div(two), two, "4/2 = 2 at p=100")
    // 2^50 (exact power of 2 — binary representation is exact in base 256)
    ae(two.powi(50), MPFloat.from_f64(F64.from[U64](1125899906842624), p), "2^50 at p=100")

    // ── Basic arithmetic — small values (|x| < 1) ─────────────────────────

    // 0.5 × 0.5 = 0.25  (exact)
    ae(half.mul(half), qtr, "0.5 × 0.5 = 0.25 at p=100")
    // 1 / 0.5 = 2       (exact)
    ae(half.inv(), two, "1/0.5 = 2 at p=100")
    // 0.25 / 0.5 = 0.5  (exact)
    ae(qtr.div(half), half, "0.25/0.5 = 0.5 at p=100")
    // (0.5)^10 = 2^{-10} = 1/1024  (exact power of 2)
    ae(half.powi(10), MPFloat.from_f64(1.0 / 1024.0, p), "0.5^10 = 1/1024 at p=100")

    // ── sqrt — large and small ─────────────────────────────────────────────

    ae(four.sqrt(), two, "sqrt(4) = 2 at p=100")
    ae(qtr.sqrt(),  half, "sqrt(0.25) = 0.5 at p=100")
    // sqrt(2) round-trip: sqrt(2)² ≈ 2
    let sqrt2 = two.sqrt()
    ae(sqrt2.mul(sqrt2), two, "sqrt(2)² ≈ 2 at p=100")
    // sqrt(0.5) round-trip: sqrt(0.5)² ≈ 0.5
    let sqrthalf = half.sqrt()
    ae(sqrthalf.mul(sqrthalf), half, "sqrt(0.5)² ≈ 0.5 at p=100")

    // ── ln / exp — large and small ─────────────────────────────────────────

    // exp(0) = 1
    ae(MPFloat.from_f64(0.0, p).exp(), one, "exp(0) = 1 at p=100")
    // exp(ln(x)) ≈ x for x > 1
    ae(two.ln().exp(), two, "exp(ln(2)) ≈ 2 at p=100")
    ae(four.ln().exp(), four, "exp(ln(4)) ≈ 4 at p=100")
    // exp(ln(x)) ≈ x for x < 1
    ae(half.ln().exp(), half, "exp(ln(0.5)) ≈ 0.5 at p=100")
    ae(qtr.ln().exp(),  qtr,  "exp(ln(0.25)) ≈ 0.25 at p=100")
    // exp(-1) × exp(1) = 1
    let e_val  = one.exp()
    let em_val = MPFloat.from_f64(-1.0, p).exp()
    ae(e_val.mul(em_val), one, "exp(1) × exp(-1) = 1 at p=100")
    // ln(a × b) = ln(a) + ln(b)
    ae(two.mul(four).ln(), two.ln().add(four.ln()), "ln(2×4) = ln(2)+ln(4) at p=100")

    // ── log2 / exp2 ────────────────────────────────────────────────────────

    // log2(2) = 1
    ae(two.log2(), one, "log2(2) = 1 at p=100")
    // exp2(log2(x)) ≈ x
    ae(two.log2().exp2(),  two,  "exp2(log2(2)) ≈ 2 at p=100")
    ae(half.log2().exp2(), half, "exp2(log2(0.5)) ≈ 0.5 at p=100")

    // ── log10 ──────────────────────────────────────────────────────────────

    // log10(10) = 1
    let ten = MPFloat.from_f64(10.0, p)
    ae(ten.log10(), one, "log10(10) = 1 at p=100")

    // ── pow ────────────────────────────────────────────────────────────────

    // 4^0.5 = 2  (via exp/ln)
    ae(four.pow(half), two, "4^0.5 = 2 at p=100")
    // 0.25^0.5 = 0.5
    ae(qtr.pow(half), half, "0.25^0.5 = 0.5 at p=100")


// ── Trigonometric functions ───────────────────────────────────────────────────

class iso _TestMPFloatTrig is UnitTest
  """
  Tests for `sin`, `cos`, `tan`, `csc`, `sec`, `cot`.
  All comparisons use `almost_eq` with matching precision.
  """
  fun name(): String => "MPFloat/trig"

  fun apply(h: TestHelper) =>
    let p: ULong = 224
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    let zero    = MPFloat.create(p)
    let one     = MPFloat.from_f64(1.0, p)
    let neg_one = MPFloat.from_f64(-1.0, p)
    let half    = MPFloat.from_f64(0.5, p)
    let two     = MPFloat.from_f64(2.0, p)

    // ── sin(0) = 0, cos(0) = 1, tan(0) = 0 ──────────────────────────────────
    ae(zero.sin(), zero, "sin(0) = 0")
    ae(zero.cos(), one,  "cos(0) = 1")
    ae(zero.tan(), zero, "tan(0) = 0")

    // ── Values at standard angles (F64 arguments keep size = p throughout) ────
    // Use F64.pi() so the argument MPFloat has exactly p bytes, ensuring the
    // result of sin/cos also has p bytes, matching the expected values.
    let pf: F64 = F64.pi()
    let pi6 = MPFloat.from_f64(pf / 6.0, p)
    let pi3 = MPFloat.from_f64(pf / 3.0, p)
    let pi4 = MPFloat.from_f64(pf / 4.0, p)
    let pi2 = MPFloat.from_f64(pf / 2.0, p)
    let pi  = MPFloat.from_f64(pf, p)

    ae(pi6.sin(), half,    "sin(π/6) = 0.5")
    ae(pi3.cos(), half,    "cos(π/3) = 0.5")
    ae(pi2.sin(), one,     "sin(π/2) = 1")
    ae(pi2.cos(), zero,    "cos(π/2) ≈ 0")
    ae(pi.cos(),  neg_one, "cos(π) = -1")
    let sin_pi = pi.sin()
    h.assert_true(sin_pi.almost_eq(zero, tol, tol),
      "sin(π) ≈ 0: got=" + sin_pi.string())
    ae(pi4.tan(), one,     "tan(π/4) = 1")

    // ── Compare against F64 reference for generic values ─────────────────────
    let v1: F64  = 1.0
    let v25: F64 = 2.5
    let x1 = MPFloat.from_f64(v1,  p)
    let x2 = MPFloat.from_f64(v25, p)
    ae(x1.sin(), MPFloat.from_f64(v1.sin(),  p), "sin(1) matches F64")
    ae(x1.cos(), MPFloat.from_f64(v1.cos(),  p), "cos(1) matches F64")
    let sin25 = x2.sin()
    let exp25 = MPFloat.from_f64(v25.sin(), p)
    h.assert_true(sin25.almost_eq(exp25, tol, tol),
      "sin(2.5) matches F64: got=" + sin25.string() + " exp=" + exp25.string())
    ae(x2.cos(), MPFloat.from_f64(v25.cos(), p), "cos(2.5) matches F64")

    // ── Pythagorean identity: sin²(x) + cos²(x) = 1 ─────────────────────────
    // Both sin and cos return p-byte results, so mul gives (2p-1)-byte
    // results; compare their sum to one (p bytes) via almost_eq.
    let s1 = x1.sin() ; let c1 = x1.cos()
    ae(s1.mul(s1).add(c1.mul(c1)), one, "sin²(1) + cos²(1) = 1")

    let s2 = x2.sin() ; let c2 = x2.cos()
    let pyth25 = s2.mul(s2).add(c2.mul(c2))
    h.assert_true(pyth25.almost_eq(one, tol, tol),
      "sin²(2.5) + cos²(2.5) = 1: got=" + pyth25.string())

    // ── Symmetry: sin(-x) = -sin(x), cos(-x) = cos(x) ───────────────────────
    ae(x1.neg().sin(), s1.neg(), "sin(-1) = -sin(1)")
    ae(x1.neg().cos(), c1,       "cos(-1) = cos(1)")

    // ── Reciprocals: csc, sec, cot ────────────────────────────────────────────
    ae(pi6.csc(), two, "csc(π/6) = 2")
    ae(pi3.sec(), two, "sec(π/3) = 2")
    ae(pi4.cot(), one, "cot(π/4) = 1")

    // ── Special values: NaN and ±∞ propagate to NaN ──────────────────────────
    let nan_val = MPFloat.nan_val()
    let pinf    = MPFloat.inf_val()
    let ninf    = MPFloat.inf_val(false)
    h.assert_true(nan_val.sin().is_nan(),  "sin(NaN) is NaN")
    h.assert_true(pinf.sin().is_nan(),     "sin(+∞) is NaN")
    h.assert_true(ninf.sin().is_nan(),     "sin(-∞) is NaN")
    h.assert_true(nan_val.cos().is_nan(),  "cos(NaN) is NaN")
    h.assert_true(pinf.cos().is_nan(),     "cos(+∞) is NaN")
    h.assert_true(nan_val.tan().is_nan(),  "tan(NaN) is NaN")
    h.assert_true(pinf.tan().is_nan(),     "tan(+∞) is NaN")
    h.assert_true(nan_val.csc().is_nan(),  "csc(NaN) is NaN")
    h.assert_true(nan_val.sec().is_nan(),  "sec(NaN) is NaN")
    h.assert_true(nan_val.cot().is_nan(),  "cot(NaN) is NaN")


// ── Hyperbolic functions ──────────────────────────────────────────────────────

class iso _TestMPFloatHyp is UnitTest
  """
  Tests for `sinh`, `cosh`, `tanh`, `csch`, `sech`, `coth`.
  All comparisons use `almost_eq` with matching precision.
  """
  fun name(): String => "MPFloat/hyp"

  fun apply(h: TestHelper) =>
    let p: ULong = 192
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    let zero = MPFloat.create(p)
    let one  = MPFloat.from_f64(1.0, p)

    // ── sinh(0) = 0, cosh(0) = 1, tanh(0) = 0 ───────────────────────────────
    ae(zero.sinh(), zero, "sinh(0) = 0")
    ae(zero.cosh(), one,  "cosh(0) = 1")
    ae(zero.tanh(), zero, "tanh(0) = 0")

    // ── Hyperbolic identity: cosh²(x) − sinh²(x) = 1 ────────────────────────
    let x1 = MPFloat.from_f64(1.0, p)
    let sh1 = x1.sinh() ; let ch1 = x1.cosh()
    ae(ch1.mul(ch1).sub(sh1.mul(sh1)), one,
      "cosh²(1) − sinh²(1) = 1")

    let x2 = MPFloat.from_f64(2.5, p)
    let sh2 = x2.sinh() ; let ch2 = x2.cosh()
    ae(ch2.mul(ch2).sub(sh2.mul(sh2)), one,
      "cosh²(2.5) − sinh²(2.5) = 1")

    // ── tanh(x) = sinh(x)/cosh(x) ────────────────────────────────────────────
    ae(x1.tanh(), sh1.div(ch1), "tanh(1) = sinh(1)/cosh(1)")

    // ── Known values from F64 ─────────────────────────────────────────────────
    // sinh(1) ≈ 1.1752011936438014, cosh(1) ≈ 1.5430806348152437
    ae(x1.sinh(), MPFloat.from_f64(F64(1.0).sinh(), p), "sinh(1) matches F64")
    ae(x1.cosh(), MPFloat.from_f64(F64(1.0).cosh(), p), "cosh(1) matches F64")
    ae(x1.tanh(), MPFloat.from_f64(F64(1.0).tanh(), p), "tanh(1) matches F64")

    // ── Symmetry: sinh(-x) = -sinh(x), cosh(-x) = cosh(x) ───────────────────
    ae(x1.neg().sinh(), sh1.neg(), "sinh(-1) = -sinh(1)")
    ae(x1.neg().cosh(), ch1,       "cosh(-1) = cosh(1)")

    // ── Reciprocals ───────────────────────────────────────────────────────────
    ae(x1.csch(), sh1.inv(), "csch(1) = 1/sinh(1)")
    ae(x1.sech(), ch1.inv(), "sech(1) = 1/cosh(1)")
    ae(x1.coth(), ch1.div(sh1), "coth(1) = cosh(1)/sinh(1)")

    // ── Special values ────────────────────────────────────────────────────────
    let nan_val = MPFloat.nan_val()
    let pinf    = MPFloat.inf_val()
    let ninf    = MPFloat.inf_val(false)

    h.assert_true(nan_val.sinh().is_nan(),     "sinh(NaN) is NaN")
    h.assert_true(pinf.sinh().is_infinite(),   "sinh(+∞) is +∞")
    h.assert_false(pinf.sinh().is_negative(),  "sinh(+∞) is positive")
    h.assert_true(ninf.sinh().is_infinite(),   "sinh(-∞) is -∞")
    h.assert_true(ninf.sinh().is_negative(),   "sinh(-∞) is negative")

    h.assert_true(nan_val.cosh().is_nan(),     "cosh(NaN) is NaN")
    h.assert_true(pinf.cosh().is_infinite(),   "cosh(+∞) is +∞")
    h.assert_false(pinf.cosh().is_negative(),  "cosh(+∞) is positive")
    h.assert_true(ninf.cosh().is_infinite(),   "cosh(-∞) is +∞")
    h.assert_false(ninf.cosh().is_negative(),  "cosh(-∞) is positive")

    h.assert_true(nan_val.tanh().is_nan(),     "tanh(NaN) is NaN")
    // tanh(±∞) = ±1
    ae(pinf.tanh(), one,      "tanh(+∞) = 1")
    ae(ninf.tanh(), one.neg(), "tanh(-∞) = -1")

    // sech(±∞) = 0
    ae(pinf.sech(), zero, "sech(+∞) = 0")
    ae(ninf.sech(), zero, "sech(-∞) = 0")

    // csch(0) = NaN (1/0)
    h.assert_true(zero.csch().is_nan(), "csch(0) is NaN")
    // coth(0) = NaN (1/0)
    h.assert_true(zero.coth().is_nan(), "coth(0) is NaN")


class iso _TestMPFloatUnsafeIntConversions is UnitTest
  """
  Verify unsafe integer conversion methods.
  """
  fun name(): String => "MPFloat/unsafe_int_conversions"

  fun apply(h: TestHelper) =>
    // I128
    let big: I128 = 0x1234567890ABCDEF1234567890ABCDEF
    try
      let mp_big = MPFloat.from_mpint(MPInt.from_string("1234567890ABCDEF1234567890ABCDEF", 16)?, 128)
      h.assert_eq[I128](big, mp_big.i128_unsafe(), "i128_unsafe(big)")
    end

    // I64
    h.assert_eq[I64](123, MPFloat.from_f64(123.456).i64_unsafe(), "i64_unsafe(123.456)")
    h.assert_eq[I64](-123, MPFloat.from_f64(-123.456).i64_unsafe(), "i64_unsafe(-123.456)")

    // I32
    h.assert_eq[I32](123, MPFloat.from_f64(123.456).i32_unsafe(), "i32_unsafe(123.456)")
    h.assert_eq[I32](-123, MPFloat.from_f64(-123.456).i32_unsafe(), "i32_unsafe(-123.456)")

    // I16
    h.assert_eq[I16](123, MPFloat.from_f64(123.456).i16_unsafe(), "i16_unsafe(123.456)")
    h.assert_eq[I16](-123, MPFloat.from_f64(-123.456).i16_unsafe(), "i16_unsafe(-123.456)")

    // I8
    h.assert_eq[I8](123, MPFloat.from_f64(123.456).i8_unsafe(), "i8_unsafe(123.456)")
    h.assert_eq[I8](-123, MPFloat.from_f64(-123.456).i8_unsafe(), "i8_unsafe(-123.456)")


class iso _TestMPFloatI8 is UnitTest
  """
  Verify i8 conversion methods.
  """
  fun name(): String => "MPFloat/i8"

  fun apply(h: TestHelper) =>
    // I8
    h.assert_eq[I8](123, MPFloat.from_f64(123.456).i8(), "i8(123.456)")
    h.assert_eq[I8](I8.max_value(), MPFloat.from_f64(200).i8(), "i8(200) saturates")
    h.assert_eq[I8](I8.min_value(), MPFloat.from_f64(-200).i8(), "i8(-200) saturates")
    h.assert_eq[I8](0, MPFloat.from_f64(0.123).i8(), "i8(0.123)")
    h.assert_eq[I8](0, MPFloat.from_f64(-0.123).i8(), "i8(-0.123)")

    // Edge cases around 2^(N-1)
    // I8: [-128, 127]
    h.assert_eq[I8](127, MPFloat.from_f64(127).i8(), "i8(127)")
    h.assert_eq[I8](I8.max_value(), MPFloat.from_f64(128).i8(), "i8(128) saturates")
    h.assert_eq[I8](-128, MPFloat.from_f64(-128).i8(), "i8(-128)")
    h.assert_eq[I8](I8.min_value(), MPFloat.from_f64(-129).i8(), "i8(-129) saturates")

    // Special values
    h.assert_eq[I8](0, MPFloat.nan_val().i8(), "i8(NaN)")
    h.assert_eq[I8](I8.max_value(), MPFloat.inf_val(true).i8(), "i8(+Inf)")
    h.assert_eq[I8](I8.min_value(), MPFloat.inf_val(false).i8(), "i8(-Inf)")


class iso _TestMPFloatI16 is UnitTest
  """
  Verify i16 conversion methods.
  """
  fun name(): String => "MPFloat/i16"

  fun apply(h: TestHelper) =>
    // I16
    h.assert_eq[I16](123, MPFloat.from_f64(123.456).i16(), "i16(123.456)")
    h.assert_eq[I16](I16.max_value(), MPFloat.from_f64(40000).i16(), "i16(40000) saturates")
    h.assert_eq[I16](I16.min_value(), MPFloat.from_f64(-40000).i16(), "i16(-40000) saturates")

    // Special values
    h.assert_eq[I16](0, MPFloat.nan_val().i16(), "i16(NaN)")
    h.assert_eq[I16](I16.max_value(), MPFloat.inf_val(true).i16(), "i16(+Inf)")
    h.assert_eq[I16](I16.min_value(), MPFloat.inf_val(false).i16(), "i16(-Inf)")


class iso _TestMPFloatI32 is UnitTest
  """
  Verify i32 conversion methods.
  """
  fun name(): String => "MPFloat/i32"

  fun apply(h: TestHelper) =>
    // I32
    h.assert_eq[I32](123, MPFloat.from_f64(123.456).i32(), "i32(123.456)")
    h.assert_eq[I32](I32.max_value(), MPFloat.from_f64(3e9).i32(), "i32(3e9) saturates")
    h.assert_eq[I32](I32.min_value(), MPFloat.from_f64(-3e9).i32(), "i32(-3e9) saturates")

    // Special values
    h.assert_eq[I32](0, MPFloat.nan_val().i32(), "i32(NaN)")
    h.assert_eq[I32](I32.max_value(), MPFloat.inf_val(true).i32(), "i32(+Inf)")
    h.assert_eq[I32](I32.min_value(), MPFloat.inf_val(false).i32(), "i32(-Inf)")


class iso _TestMPFloatI64 is UnitTest
  """
  Verify i64 conversion methods.
  """
  fun name(): String => "MPFloat/i64"

  fun apply(h: TestHelper) =>
    // I64
    h.assert_eq[I64](123, MPFloat.from_f64(123.456).i64(), "i64(123.456)")
    h.assert_eq[I64](I64.max_value(), MPFloat.from_f64(2e19).i64(), "i64(2e19) saturates")
    h.assert_eq[I64](I64.min_value(), MPFloat.from_f64(-2e19).i64(), "i64(-2e19) saturates")

    // Truncation of decimals
    h.assert_eq[I64](1, MPFloat.from_f64(1.9).i64(), "i64(1.9)")
    h.assert_eq[I64](-1, MPFloat.from_f64(-1.9).i64(), "i64(-1.9)")
    h.assert_eq[I64](0, MPFloat.from_f64(0.5).i64(), "i64(0.5)")
    h.assert_eq[I64](0, MPFloat.from_f64(-0.5).i64(), "i64(-0.5)")

    // Large values (within I64 range)
    let big: I64 = 0x1234567890ABCDEF
    try
      let mp_big = MPFloat.from_mpint(MPInt.from_string("1234567890ABCDEF", 16)?, 64)
      h.assert_eq[I64](big, mp_big.i64(), "i64(big)")

      let neg_big: I64 = -big
      let mp_neg_big = MPFloat.from_mpint(MPInt.from_string("-1234567890ABCDEF", 16)?, 64)
      h.assert_eq[I64](neg_big, mp_neg_big.i64(), "i64(neg_big)")

      // Overflow behavior (saturation like F64)
      // 2^64 + 1 saturates to I64.max_value()
      let overflow = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](64))) + MPInt.from[ILong](1)
      let mp_overflow = MPFloat.from_mpint(overflow, 68)
      h.assert_eq[I64](I64.max_value(), mp_overflow.i64(), "i64(2^64 + 1) saturates")

      // Boundary checks around 2^63
      let p2_63_minus_1 = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](63))) - MPInt.from[ILong](1)
      h.assert_eq[I64](I64.max_value(), MPFloat.from_mpint(p2_63_minus_1, 64).i64(), "i64(2^63 - 1)")

      let p2_63 = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](63)))
      h.assert_eq[I64](I64.max_value(), MPFloat.from_mpint(p2_63, 64).i64(), "i64(2^63) saturates")

      let n2_63 = -p2_63
      h.assert_eq[I64](I64.min_value(), MPFloat.from_mpint(n2_63, 64).i64(), "i64(-2^63)")

      let n2_63_minus_1 = n2_63 - MPInt.from[ILong](1)
      h.assert_eq[I64](I64.min_value(), MPFloat.from_mpint(n2_63_minus_1, 64).i64(), "i64(-2^63 - 1) saturates")

      // Huge values
      let huge = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](1000)))
      h.assert_eq[I64](I64.max_value(), MPFloat.from_mpint(huge).i64(), "i64(2^1000) saturates")
      h.assert_eq[I64](I64.min_value(), MPFloat.from_mpint(-huge).i64(), "i64(-2^1000) saturates")
    else
      h.fail("MPInt arithmetic failed in tests")
    end

    // Special values
    h.assert_eq[I64](0, MPFloat.nan_val().i64(), "i64(NaN)")
    h.assert_eq[I64](I64.max_value(), MPFloat.inf_val(true).i64(), "i64(+Inf)")
    h.assert_eq[I64](I64.min_value(), MPFloat.inf_val(false).i64(), "i64(-Inf)")


class iso _TestMPFloatI128 is UnitTest
  """
  Verify i128() conversion method.
  """
  fun name(): String => "MPFloat/i128"

  fun apply(h: TestHelper) =>
    // Basic integers
    h.assert_eq[I128](0, MPFloat.from_f64(0.0).i128(), "i128(0.0)")
    h.assert_eq[I128](1, MPFloat.from_f64(1.0).i128(), "i128(1.0)")
    h.assert_eq[I128](-1, MPFloat.from_f64(-1.0).i128(), "i128(-1.0)")
    h.assert_eq[I128](123456789, MPFloat.from_f64(123456789.0).i128(), "i128(123456789.0)")

    // Truncation of decimals
    h.assert_eq[I128](1, MPFloat.from_f64(1.9).i128(), "i128(1.9)")
    h.assert_eq[I128](-1, MPFloat.from_f64(-1.9).i128(), "i128(-1.9)")
    h.assert_eq[I128](0, MPFloat.from_f64(0.5).i128(), "i128(0.5)")
    h.assert_eq[I128](0, MPFloat.from_f64(-0.5).i128(), "i128(-0.5)")

    // Large values (within I128 range)
    let big: I128 = 0x1234567890ABCDEF1234567890ABCDEF
    try
      let mp_big = MPFloat.from_mpint(MPInt.from_string("1234567890ABCDEF1234567890ABCDEF", 16)?, 128)
      h.assert_eq[I128](big, mp_big.i128(), "i128(big)")

      let neg_big: I128 = -big
      let mp_neg_big = MPFloat.from_mpint(MPInt.from_string("-1234567890ABCDEF1234567890ABCDEF", 16)?, 128)
      h.assert_eq[I128](neg_big, mp_neg_big.i128(), "i128(neg_big)")

      // Overflow behavior (saturation like F64)
      // 2^128 + 1 saturates to I128.max_value()
      let overflow = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](128))) + MPInt.from[ILong](1)
      let mp_overflow = MPFloat.from_mpint(overflow, 136)
      h.assert_eq[I128](I128.max_value(), mp_overflow.i128(), "i128(2^128 + 1) saturates")

      // Boundary checks around 2^127
      let p2_127_minus_1 = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](127))) - MPInt.from[ILong](1)
      h.assert_eq[I128](I128.max_value(), MPFloat.from_mpint(p2_127_minus_1, 128).i128(), "i128(2^127 - 1)")

      let p2_127 = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](127)))
      h.assert_eq[I128](I128.max_value(), MPFloat.from_mpint(p2_127, 128).i128(), "i128(2^127) saturates")

      let n2_127 = -p2_127
      h.assert_eq[I128](I128.min_value(), MPFloat.from_mpint(n2_127, 128).i128(), "i128(-2^127)")

      let n2_127_minus_1 = n2_127 - MPInt.from[ILong](1)
      h.assert_eq[I128](I128.min_value(), MPFloat.from_mpint(n2_127_minus_1, 128).i128(), "i128(-2^127 - 1) saturates")

      // Huge values
      let huge = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](1000)))
      h.assert_eq[I128](I128.max_value(), MPFloat.from_mpint(huge).i128(), "i128(2^1000) saturates")
      h.assert_eq[I128](I128.min_value(), MPFloat.from_mpint(-huge).i128(), "i128(-2^1000) saturates")
    else
      h.fail("MPInt arithmetic failed in tests")
    end

    // Special values
    h.assert_eq[I128](0, MPFloat.nan_val().i128(), "i128(NaN)")
    h.assert_eq[I128](I128.max_value(), MPFloat.inf_val(true).i128(), "i128(+Inf)")
    h.assert_eq[I128](I128.min_value(), MPFloat.inf_val(false).i128(), "i128(-Inf)")

/* TODO RESTORE
    // Random values
    let rand = Rand
    for i in Range(0, 100) do
      let t = rand.i128()
      h.assert_eq[I128](t, MPFloat.from_ilong(t.ilong()).i128(), "Random [" + i.string() + "] i128(" + t.string() + ")")
    end
*/

class iso _TestMPFloatPi is UnitTest
  """
  Tests for π calculations.
  Machin's formula agains Chudnovsky's.
  All comparisons use `almost_eq` with matching precision.
  """
  fun name(): String => "MPFloat/pi"

  fun apply(h: TestHelper) =>
    var p: ULong = 800
    // 100-byte precision → ~240 decimal digits.  We demand 100-digit accuracy.
    let tol: MPFloat = MPFloat.epsilon(p).sqrt()
    let ae = {(got: MPFloat, expected: MPFloat, msg: String) =>
      h.assert_true(got.almost_eq(expected, tol, tol), msg)
    }

    let pi_machin = MPFloat.pi(p)
    (let machin_m, let machin_e, let machin_i) = pi_machin.exact_string()
    let mm: String val = consume machin_m
    h.log("Pi Machin = " + pi_machin.string())

    let pi_chudnovsky = MPFloat.pi_chudnovsky(p)
    (let chudnovsky_m, let chudnovsky_e, let chudnovsky_i) = pi_chudnovsky.exact_string()
    let cm: String val = consume chudnovsky_m
    h.log("Pi Chudnovsky = " + pi_chudnovsky.string())

    let diff = pi_chudnovsky - pi_machin
    h.log("Difference = " + diff.string())
    h.log("")

    let pi_bbp = MPFloat.pi_bbp(p)
    (let bbp_m, let bbp_e, let bbp_i) = pi_bbp.exact_string()
    let bm: String val = consume bbp_m
    h.log("Pi BBP = " + pi_bbp.string())

    let diff2 = pi_bbp - pi_machin
    h.log("Difference = " + diff2.string())
    h.log("")

    h.assert_eq[String](mm, cm, "Pi mantissas are different:\nMachin     = " +
      mm + "\nChudnovsky = " + cm)
    h.assert_eq[I64](machin_e, chudnovsky_e, "Pi exponents are different:\nMachin     = " +
      machin_e.string() + "\nChudnovsky = " + chudnovsky_e.string())


class iso _TestMPFloatConversions is UnitTest
  """
  Verify f64() and f32() conversion methods.
  """
  fun name(): String => "MPFloat/conversions"

  fun apply(h: TestHelper) =>
    // Specific values round-trip
    let vals64: Array[F64] = [0.0; 1.0; -1.0; 0.5; -0.5; 0.25; 123.456; 1e10; 1e-10; F64.min_normalised(); F64.epsilon(); F64.min_value(); F64.max_value(); F64.pi(); F64.e()]
    for v in vals64.values() do
      let mpf = MPFloat.from_f64(v)
      h.assert_eq[F64](v, mpf.f64(), "F64 round-trip failed for " + v.string())
    end

    let vals32: Array[F32] = [0.0; 1.0; -1.0; 0.5; -0.5; 0.25; 123.456; 1e5; 1e-5; F32.min_normalised(); F32.epsilon(); F32.min_value(); F32.max_value(); F32.pi(); F32.e()]
    for v in vals32.values() do
      let mpf = MPFloat.from_f32(v)
      h.assert_eq[F32](v, mpf.f32(), "F32 round-trip failed for " + v.string())
    end

    // Special values
    h.assert_true(MPFloat.nan_val().f64().nan(), "NaN round-trip F64")
    h.assert_true(MPFloat.inf_val(true).f64().infinite(), "Inf round-trip F64")
    h.assert_true(MPFloat.inf_val(true).f64() > 0, "Inf sign F64")
    h.assert_true(MPFloat.inf_val(false).f64().infinite(), "-Inf round-trip F64")
    h.assert_true(MPFloat.inf_val(false).f64() < 0, "-Inf sign F64")

    h.assert_true(MPFloat.nan_val().f32().nan(), "NaN round-trip F32")
    h.assert_true(MPFloat.inf_val(true).f32().infinite(), "Inf round-trip F32")
    h.assert_true(MPFloat.inf_val(true).f32() > 0, "Inf sign F32")
    h.assert_true(MPFloat.inf_val(false).f32().infinite(), "-Inf round-trip F32")
    h.assert_true(MPFloat.inf_val(false).f32() < 0, "-Inf sign F32")

    // Negative zero
    let neg_zero_f64 = MPFloat.from_f64(-0.0)
    h.assert_eq[F64](-0.0, neg_zero_f64.f64(), "-0.0 round-trip F64")
    h.assert_true(neg_zero_f64.f64().bits() == 0x8000000000000000, "-0.0 sign bit F64")

    let neg_zero_f32 = MPFloat.from_f32(-0.0)
    h.assert_eq[F32](-0.0, neg_zero_f32.f32(), "-0.0 round-trip F32")
    h.assert_true(neg_zero_f32.f32().bits() == 0x80000000, "-0.0 sign bit F32")

    // Random values
    let rand = Rand
    for i in Range(0, 100) do
      let f = F64.from_bits(rand.u64())
      let rt_f = MPFloat.from_f64(f).f64()
      if f.nan() then
        h.assert_true(rt_f.nan(), "Random round-trip [" + i.string() + "] failed for F64=NaN")
      elseif f.infinite() then
        h.assert_true(rt_f.infinite(), "Random round-trip [" + i.string() + "] failed for F64=" + f.string())
        h.assert_eq[Bool](f < 0, rt_f < 0, "Random round-trip [" + i.string() + "] sign mismatch for F64=" + f.string())
      elseif f == 0.0 then
        h.assert_eq[U64](f.bits(), rt_f.bits(), "Random round-trip [" + i.string() + "] sign mismatch for F64=0")
      else
        // Test on 16 ulp
        h.assert_true((f - rt_f).abs() <= (f.abs().max(1.0) * f.epsilon() * 16), "Random round-trip [" + i.string() + "] failed for F64=" + f.string() + " got " + rt_f.string())
      end
    end
    for i in Range(0, 100) do
      let f = F32.from_bits(rand.u32())
      let rt_f = MPFloat.from_f32(f).f32()
      if f.nan() then
        h.assert_true(rt_f.nan(), "Random round-trip [" + i.string() + "] failed for F32=NaN")
      elseif f.infinite() then
        h.assert_true(rt_f.infinite(), "Random round-trip [" + i.string() + "] failed for F32=" + f.string())
        h.assert_eq[Bool](f < 0, rt_f < 0, "Random round-trip [" + i.string() + "] sign mismatch for F32=" + f.string())
      elseif f == 0 then
        h.assert_eq[U32](f.bits(), rt_f.bits(), "Random round-trip [" + i.string() + "] sign mismatch for F32=0")
      else
        // Test on 16 ulp
        h.assert_true((f - rt_f).abs() <= (f.abs().max(1.0) * f.epsilon() * 16), "Random round-trip [" + i.string() + "] failed for F32=" + f.string() + " got " + rt_f.string())
      end
    end
    for i in Range(0, 100) do
      let mpf = MPFloat.from_f64(F64.from_bits(rand.u64()))
      let rt_mpf = MPFloat.from_f64(mpf.f64())
      if mpf.is_nan() then
        h.assert_true(rt_mpf.is_nan(), "Random round-trip [" + i.string() + "] failed for MPFloat.from_f64=NaN")
      elseif mpf.is_infinite() then
        h.assert_true(rt_mpf.is_infinite(), "Random round-trip [" + i.string() + "] failed for MPFloat.from_f64=" + mpf.string())
        h.assert_eq[Bool](mpf.is_negative(), rt_mpf.is_negative(), "Random round-trip [" + i.string() + "] sign mismatch for MPFloat.from_f64=" + mpf.string())
      else
        // Use almost_eq for finite numbers
        h.assert_true(mpf.almost_eq(rt_mpf), "Random round-trip [" + i.string() + "] failed for MPFloat.from_f64=" + mpf.string() + " got " + rt_mpf.string())
      end
    end
    for i in Range(0, 100) do
      let mpf = MPFloat.from_f32(F32.from_bits(rand.u32()))
      let rt_mpf = MPFloat.from_f32(mpf.f32())
      if mpf.is_nan() then
        h.assert_true(rt_mpf.is_nan(), "Random round-trip [" + i.string() + "] failed for MPFloat.from_f32=NaN")
      elseif mpf.is_infinite() then
        h.assert_true(rt_mpf.is_infinite(), "Random round-trip [" + i.string() + "] failed for MPFloat.from_f32=" + mpf.string())
        h.assert_eq[Bool](mpf.is_negative(), rt_mpf.is_negative(), "Random round-trip [" + i.string() + "] sign mismatch for MPFloat.from_f32=" + mpf.string())
      else
        h.assert_true(mpf.almost_eq(rt_mpf), "Random round-trip [" + i.string() + "] failed for MPFloat.from_f32=" + mpf.string() + " got " + rt_mpf.string())
      end
    end

/* TODO RESTORE
class iso _TestMPFloatU8 is UnitTest
  """
  Verify u8 conversion methods.
  """
  fun name(): String => "MPFloat/u8"

  fun apply(h: TestHelper) =>
    // U8
    h.assert_eq[U8](123, MPFloat.from_f64(123.456).u8(), "u8(123.456)")
    h.assert_eq[U8](U8.max_value(), MPFloat.from_f64(400).u8(), "u8(400) saturates")
    h.assert_eq[U8](U8.min_value(), MPFloat.from_f64(-200).u8(), "u8(-200) saturates")
    h.assert_eq[U8](0, MPFloat.from_f64(0.123).u8(), "u8(0.123)")
    h.assert_eq[U8](0, MPFloat.from_f64(-0.123).u8(), "u8(-0.123)")

    // Edge cases around 2^N
    // U8: [0, 255]
    h.assert_eq[U8](255, MPFloat.from_f64(255).u8(), "u8(255)")
    h.assert_eq[U8](U8.max_value(), MPFloat.from_f64(256).u8(), "u8(256) saturates")
    h.assert_eq[U8](0, MPFloat.from_f64(0).u8(), "u8(0)")
    h.assert_eq[U8](U8.min_value(), MPFloat.from_f64(-1).u8(), "u8(-1) saturates")

    // Special values
    h.assert_eq[U8](0, MPFloat.nan_val().u8(), "u8(NaN)")
    h.assert_eq[U8](U8.max_value(), MPFloat.inf_val(true).u8(), "u8(+Inf)")
    h.assert_eq[U8](U8.min_value(), MPFloat.inf_val(false).u8(), "u8(-Inf)")


class iso _TestMPFloatU16 is UnitTest
  """
  Verify u16 conversion methods.
  """
  fun name(): String => "MPFloat/u16"

  fun apply(h: TestHelper) =>
    // U16
    h.assert_eq[U16](123, MPFloat.from_f64(123.456).u16(), "u16(123.456)")
    h.assert_eq[U16](U16.max_value(), MPFloat.from_f64(400000).u16(), "u16(400000) saturates")
    h.assert_eq[U16](U16.min_value(), MPFloat.from_f64(-40000).u16(), "u16(-40000) saturates")

    // Special values
    h.assert_eq[U16](0, MPFloat.nan_val().u16(), "u16(NaN)")
    h.assert_eq[U16](U16.max_value(), MPFloat.inf_val(true).u16(), "u16(+Inf)")
    h.assert_eq[U16](U16.min_value(), MPFloat.inf_val(false).u16(), "u16(-Inf)")


class iso _TestMPFloatU32 is UnitTest
  """
  Verify u32 conversion methods.
  """
  fun name(): String => "MPFloat/u32"

  fun apply(h: TestHelper) =>
    // U32
    h.assert_eq[U32](123, MPFloat.from_f64(123.456).u32(), "u32(123.456)")
    h.assert_eq[U32](U32.max_value(), MPFloat.from_f64(3e12).u32(), "u32(3e12) saturates")
    h.assert_eq[U32](U32.min_value(), MPFloat.from_f64(-3e9).u32(), "u32(-3e9) saturates")

    // Special values
    h.assert_eq[U32](0, MPFloat.nan_val().u32(), "u32(NaN)")
    h.assert_eq[U32](U32.max_value(), MPFloat.inf_val(true).u32(), "u32(+Inf)")
    h.assert_eq[U32](U32.min_value(), MPFloat.inf_val(false).u32(), "u32(-Inf)")


class iso _TestMPFloatU64 is UnitTest
  """
  Verify u64 conversion methods.
  """
  fun name(): String => "MPFloat/u64"

  fun apply(h: TestHelper) =>
    // U64
    h.assert_eq[U64](123, MPFloat.from_f64(123.456).u64(), "u64(123.456)")
    h.assert_eq[U64](U64.max_value(), MPFloat.from_f64(2e24).u64(), "u64(2e19) saturates")
    h.assert_eq[U64](U64.min_value(), MPFloat.from_f64(-2e19).u64(), "u64(-2e19) saturates")

    // Truncation of decimals
    h.assert_eq[U64](1, MPFloat.from_f64(1.9).u64(), "u64(1.9)")
    h.assert_eq[U64](0, MPFloat.from_f64(-1.9).u64(), "u64(-1.9)")
    h.assert_eq[U64](0, MPFloat.from_f64(0.5).u64(), "u64(0.5)")
    h.assert_eq[U64](0, MPFloat.from_f64(-0.5).u64(), "u64(-0.5)")

    // Large values (within U64 range)
    let big: U64 = 0x1234567890ABCDEF
    try
      let mp_big = MPFloat.from_mpint(MPInt.from_string("1234567890ABCDEF", 16)?, 64)
      h.assert_eq[U64](big, mp_big.u64(), "u64(big)")

      // Overflow behavior (saturation like F64)
      // 2^64 + 1 saturates to U64.max_value()
      let overflow = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](64))) + MPInt.from[ILong](1)
      let mp_overflow = MPFloat.from_mpint(overflow, 68)
      h.assert_eq[U64](U64.max_value(), mp_overflow.u64(), "u64(2^64 + 1) saturates")

      // Boundary checks around 2^64
      let p2_64_minus_1 = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](64))) - MPInt.from[ILong](1)
      h.assert_eq[U64](U64.max_value(), MPFloat.from_mpint(p2_64_minus_1, 68).u64(), "u64(2^64 - 1)")

      let p2_64 = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](64)))
      h.assert_eq[U64](U64.max_value(), MPFloat.from_mpint(p2_64, 64).u64(), "u64(2^64) max value")

      let n2_64_plus_1 = p2_64 + MPInt.from[ILong](1)
      h.assert_eq[U64](U64.max_value(), MPFloat.from_mpint(n2_64_plus_1, 68).u64(), "u64(2^64 + 1)")

      // Huge values
      let huge = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](1000)))
      h.assert_eq[U64](U64.max_value(), MPFloat.from_mpint(huge).u64(), "u64(2^1000) saturates")
    else
      h.fail("MPInt arithmetic failed in tests")
    end

    // Special values
    h.assert_eq[U64](0, MPFloat.nan_val().u64(), "u64(NaN)")
    h.assert_eq[U64](U64.max_value(), MPFloat.inf_val(true).u64(), "u64(+Inf)")
    h.assert_eq[U64](U64.min_value(), MPFloat.inf_val(false).u64(), "u64(-Inf)")


class iso _TestMPFloatU128 is UnitTest
  """
  Verify u128() conversion method.
  """
  fun name(): String => "MPFloat/u128"

  fun apply(h: TestHelper) =>
    // Basic integers
    h.assert_eq[U128](0, MPFloat.from_f64(0.0).u128(), "u128(0.0)")
    h.assert_eq[U128](1, MPFloat.from_f64(1.0).u128(), "u128(1.0)")
    h.assert_eq[U128](U128.min_value(), MPFloat.from_f64(-1.0).u128(), "u128(-1.0)")
    h.assert_eq[U128](123456789, MPFloat.from_f64(123456789.0).u128(), "u128(123456789.0)")

    // Truncation of decimals
    h.assert_eq[U128](1, MPFloat.from_f64(1.9).u128(), "u128(1.9)")
    h.assert_eq[U128](U128.min_value(), MPFloat.from_f64(-1.9).u128(), "u128(-1.9)")
    h.assert_eq[U128](0, MPFloat.from_f64(0.5).u128(), "u128(0.5)")
    h.assert_eq[U128](U128.min_value(), MPFloat.from_f64(-0.5).u128(), "u128(-0.5)")

    // Large values (within U128 range)
    let big: U128 = 0x1234567890ABCDEF1234567890ABCDEF
    try
      let mp_big = MPFloat.from_mpint(MPInt.from_string("1234567890ABCDEF1234567890ABCDEF", 16)?, 128)
      h.assert_eq[U128](big, mp_big.u128(), "u128(big)")

      let mp_neg_big = MPFloat.from_mpint(MPInt.from_string("-1234567890ABCDEF1234567890ABCDEF", 16)?, 128)
      h.assert_eq[U128](U128.min_value(), mp_neg_big.u128(), "u128(neg_big)")

      // Overflow behavior (saturation like F64)
      // 2^128 + 1 saturates to U128.max_value()
      let overflow = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](128))) + MPInt.from[ILong](1)
      let mp_overflow = MPFloat.from_mpint(overflow, 136)
      h.assert_eq[U128](U128.max_value(), mp_overflow.u128(), "u128(2^128 + 1) saturates")

      // Boundary checks around 2^128
      let p2_128 = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](128)))
      h.assert_eq[U128](U128.max_value(), MPFloat.from_mpint(p2_128, 128).u128(), "u128(2^128 saturates)")

      let p2_128_minus_1 = p2_128 - MPInt.from[ILong](1)
      h.assert_eq[U128](U128.max_value(), MPFloat.from_mpint(p2_128_minus_1, 128).u128(), "u128(2^128 - 1)")

      let p2_128_minus_2 = p2_128_minus_1 - MPInt.from[ILong](1)
      h.assert_eq[U128](U128.max_value() - 1, MPFloat.from_mpint(p2_128_minus_2, 128).u128(), "u128(-2^128 - 2)")

      // Huge values
      let huge = (MPInt.from[ULong](1).bit_shl(MPInt.from[ILong](1000)))
      h.assert_eq[U128](U128.max_value(), MPFloat.from_mpint(huge).u128(), "u128(2^1000) saturates")
      h.assert_eq[U128](U128.min_value(), MPFloat.from_mpint(-huge).u128(), "u128(-2^1000) saturates")
    else
      h.fail("MPInt arithmetic failed in tests")
    end

    // Special values
    h.assert_eq[U128](0, MPFloat.nan_val().u128(), "u128(NaN)")
    h.assert_eq[U128](U128.max_value(), MPFloat.inf_val(true).u128(), "u128(+Inf)")
    h.assert_eq[U128](U128.min_value(), MPFloat.inf_val(false).u128(), "u128(-Inf)")

    // Random values
    let rand = Rand
    for i in Range(0, 100) do
      let t = rand.u128()
      h.assert_eq[U128](t, MPFloat.from[U128](t).u128(), "Random [" + i.string() + "] u128(" + t.string() + ")")
    end
*/
