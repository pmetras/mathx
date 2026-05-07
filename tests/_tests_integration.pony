use "time"

use "../mathx"
use "../pony_testx"


// First 1000 significant decimal digits of π (no decimal point).
primitive _PiKnown1000
  fun apply(): String =>
    "3141592653589793238462643383279502884197169399375105820974944592307816406286" +
    "2089986280348253421170679821480865132823066470938446095505822317253594081284" +
    "8111745028410270193852110555964462294895493038196442881097566593344612847564" +
    "8233786783165271201909145648566923460348610454326648213393607260249141273724" +
    "5870066063155881748815209209628292540917153643678925903600113305305488204665" +
    "2138414695194151160943305727036575959195309218611738193261179310511854807446" +
    "2379962749567351885752724891227938183011949129833673362440656643086021394946" +
    "3952247371907021798609437027705392171762931767523846748184676694051320005681" +
    "2714526356082778577134275778960917363717872146844090122495343014654958537105" +
    "0792279689258923542019956112129021960864034418159813629774771309960518707211" +
    "3499999983729780499510597317328160963185950244594553469083026425223082533446" +
    "8503526193118817101000313783875288658753320838142061717766914730359825349042" +
    "8755468731159562863882353787593751957781857780532171226806613001927876611195" +
    "909216420198"


// Precision required for 1000 decimal digits:
// ceil(1000 / log10(256)) ≈ 416 bytes = 3328 bits; use 3500 for guard digits.
// With the proportional exact_string gate (gate = _size() + 150 = 587) and
// k_bytes = _exponent − _size() = 1 − 437 = −436, the exact MPInt path is used.
primitive _Pi1000Prec
  fun apply(): USize => 3500


class iso _TestMPFloatPi1000Machin is UnitTest
  """
  High-precision π via Machin's formula to 1000 decimal digits.
  """
  fun name(): String => "integration/pi1000/machin"
  fun label(): String => "integration"

  fun apply(h: TestHelper) =>
    h.long_test(600_000_000_000) // 10-minute timeout

    let t0 = Time.nanos()
    let result = MPFMathLib.pi_machin(_Pi1000Prec())
    let elapsed_ms = (Time.nanos() - t0).f64() / 1_000_000.0
    h.log("Pi Machin elapsed: " + elapsed_ms.string() + " ms")

    (let m_iso, let exp, let inexact) = result.exact_string()
    let m: String val = consume m_iso
    h.log("Pi Machin (first 50): " + m.substring(0, 50))
    h.assert_false(inexact, "Pi Machin: exact_string must use exact path (inexact = false)")
    h.assert_true(
      m.at(_PiKnown1000(), 0),
      "Pi Machin 1000-digit mismatch: got \"" + m.substring(0, 50) +
      "\", dec_exp=" + exp.string())
    h.complete(true)


class iso _TestMPFloatPi1000BBP is UnitTest
  """
  High-precision π via Bailey–Borwein–Plouffe (BBP) formula to 1000 decimal digits.
  """
  fun name(): String => "integration/pi1000/bbp"
  fun label(): String => "integration"

  fun apply(h: TestHelper) =>
    h.long_test(600_000_000_000)

    let t0 = Time.nanos()
    let result = MPFMathLib.pi_bbp(_Pi1000Prec())
    let elapsed_ms = (Time.nanos() - t0).f64() / 1_000_000.0
    h.log("Pi BBP elapsed: " + elapsed_ms.string() + " ms")

    (let m_iso, let exp, let inexact) = result.exact_string()
    let m: String val = consume m_iso
    h.log("Pi BBP (first 50): " + m.substring(0, 50))
    h.assert_false(inexact, "Pi BBP: exact_string must use exact path (inexact = false)")
    h.assert_true(
      m.at(_PiKnown1000(), 0),
      "Pi BBP 1000-digit mismatch: got \"" + m.substring(0, 50) +
      "\", dec_exp=" + exp.string())
    h.complete(true)


class iso _TestMPFloatPi1000Chudnovsky is UnitTest
  """
  High-precision π via Chudnovsky series to 1000 decimal digits.
  """
  fun name(): String => "integration/pi1000/chudnovsky"
  fun label(): String => "integration"

  fun apply(h: TestHelper) =>
    h.long_test(600_000_000_000)

    let t0 = Time.nanos()
    let result = MPFloat.pi(_Pi1000Prec())
    let elapsed_ms = (Time.nanos() - t0).f64() / 1_000_000.0
    h.log("Pi Chudnovsky elapsed: " + elapsed_ms.string() + " ms")

    (let m_iso, let exp, let inexact) = result.exact_string()
    let m: String val = consume m_iso
    h.log("Pi Chudnovsky (first 50): " + m.substring(0, 50))
    h.assert_false(inexact, "Pi Chudnovsky: exact_string must use exact path (inexact = false)")
    h.assert_true(
      m.at(_PiKnown1000(), 0),
      "Pi Chudnovsky 1000-digit mismatch: got \"" + m.substring(0, 50) +
      "\", dec_exp=" + exp.string())
    h.complete(true)
