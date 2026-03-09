
use "../mathx"
use "../pony_testx"
use "collections"
use "random"

class iso _TestMPIntFldMod is UnitTest
  """
  Tests for floored division (fld) and modulo (mod).
  Property: a = fld(a, b) * b + mod(a, b)
  where mod(a, b) has the same sign as b.
  """

  fun name(): String =>
    "MPInt/fld_mod"

  fun apply(h: TestHelper) =>
    let rand = Rand()
    let zero = MPInt.from_ilong(0)

    for _ in Range(0, 100) do
      let a = _random_mpint(rand, 10)
      let b = _random_mpint(rand, 5)
      
      if b.is_zero() then continue end
      
      let q = a.fld(b)
      let r = a.mod(b)
      
      // Property check: a = q*b + r
      let check = (q * b) + r
      h.assert_true(check == a, "fld/mod property: a = q*b + r")
      
      // Sign check: r is 0 or sign(r) == sign(b)
      if not r.is_zero() then
        h.assert_true((r < zero) == (b < zero), "mod sign matches divisor sign")
      end
      
      // Magnitude check: |r| < |b|
      h.assert_true(r.abs_lt(b.abs()), "mod magnitude: |r| < |b|")
    end

    // Specific cases
    let a1 = MPInt.from_ilong(10)
    let b1 = MPInt.from_ilong(3)
    h.assert_true(a1.fld(b1) == MPInt.from_ilong(3), "10 fld 3 = 3")
    h.assert_true(a1.mod(b1) == MPInt.from_ilong(1), "10 mod 3 = 1")

    let a2 = MPInt.from_ilong(-10)
    let b2 = MPInt.from_ilong(3)
    h.assert_true(a2.fld(b2) == MPInt.from_ilong(-4), "-10 fld 3 = -4")
    h.assert_true(a2.mod(b2) == MPInt.from_ilong(2), "-10 mod 3 = 2")

    let a3 = MPInt.from_ilong(10)
    let b3 = MPInt.from_ilong(-3)
    h.assert_true(a3.fld(b3) == MPInt.from_ilong(-4), "10 fld -3 = -4")
    h.assert_true(a3.mod(b3) == MPInt.from_ilong(-2), "10 mod -3 = -2")

    let a4 = MPInt.from_ilong(-10)
    let b4 = MPInt.from_ilong(-3)
    h.assert_true(a4.fld(b4) == MPInt.from_ilong(3), "-10 fld -3 = 3")
    h.assert_true(a4.mod(b4) == MPInt.from_ilong(-1), "-10 mod -3 = -1")

  fun _random_mpint(rand: Rand, size: USize): MPInt =>
    var res = MPInt.from_ilong(rand.ilong())
    for i in Range(1, size) do
      res = res + (MPInt.from_ilong(rand.ilong()).abs().digit_shl(i))
    end
    if (rand.next() % 2) == 0 then
      res.neg()
    else
      res
    end

class iso _TestMPIntPopcount is UnitTest
  """
  Tests for popcount method.
  """

  fun name(): String =>
    "MPInt/popcount"

  fun apply(h: TestHelper) =>
    h.assert_true(MPInt.from_ilong(0).popcount() == MPInt.from_ilong(0), "popcount(0)")
    h.assert_true(MPInt.from_ilong(1).popcount() == MPInt.from_ilong(1), "popcount(1)")
    h.assert_true(MPInt.from_ilong(7).popcount() == MPInt.from_ilong(3), "popcount(7)")
    h.assert_true(MPInt.from_ilong(255).popcount() == MPInt.from_ilong(8), "popcount(255)")
    h.assert_true(MPInt.from_ilong(-7).popcount() == MPInt.from_ilong(3), "popcount(-7) uses absolute value")
    // 0b1010_1100 = 172: 4 bits set
    h.assert_true(MPInt.from_ilong(172).popcount() == MPInt.from_ilong(4), "popcount(172) = 4")

    // Large value: 0xFFFF_FFFF_FFFF_FFFF
    try
      let large = MPInt.from_string("18446744073709551615")?
      h.assert_true(large.popcount() == MPInt.from_ilong(64), "popcount(2^64 - 1)")
    end

    // ctz: trailing-zero count (operates on absolute value)
    let zero = MPInt.from_ilong(0)
    h.assert_true(zero.ctz().is_zero(),                                   "ctz(0) = 0")
    h.assert_true(MPInt.from_ilong(1).ctz()      == zero,                 "ctz(1) = 0")
    h.assert_true(MPInt.from_ilong(2).ctz()      == MPInt.from_ilong(1),  "ctz(2) = 1")
    h.assert_true(MPInt.from_ilong(4).ctz()      == MPInt.from_ilong(2),  "ctz(4) = 2")
    h.assert_true(MPInt.from_ilong(8).ctz()      == MPInt.from_ilong(3),  "ctz(8) = 3")
    h.assert_true(MPInt.from_ilong(12).ctz()     == MPInt.from_ilong(2),  "ctz(12) = 2")
    // 65536 = 2^16: crosses digit boundary
    h.assert_true(MPInt.from_ilong(65536).ctz()  == MPInt.from_ilong(16), "ctz(65536) = 16")
    h.assert_true((MPInt.from_ilong(65536) * MPInt.from_ilong(4)).ctz() == MPInt.from_ilong(18), "ctz(4*65536) = 18")
    // Absolute value
    h.assert_true(MPInt.from_ilong(-4).ctz() == MPInt.from_ilong(4).ctz(), "ctz(-x) = ctz(x)")
    // ctz_unsafe == ctz
    h.assert_true(MPInt.from_ilong(12).ctz_unsafe() == MPInt.from_ilong(12).ctz(), "ctz_unsafe == ctz")

class iso _TestMPIntBytewidth is UnitTest
  """
  Tests for bytewidth method.
  """

  fun name(): String =>
    "MPInt/bytewidth"

  fun apply(h: TestHelper) =>
    // 0 -> 0 bits -> 0 bytes
    h.assert_true(MPInt.from_ilong(0).bytewidth() == 0, "bytewidth(0)")
    
    // 1 -> 1 bit + 1 sign bit = 2 bits -> 1 byte
    h.assert_true(MPInt.from_ilong(1).bytewidth() == 1, "bytewidth(1)")
    
    // 127 (0x7F) -> 7 bits + 1 sign bit = 8 bits -> 1 byte
    h.assert_true(MPInt.from_ilong(127).bytewidth() == 1, "bytewidth(127)")
    
    // 128 (0x80) -> 8 bits + 1 sign bit = 9 bits -> 2 bytes
    h.assert_true(MPInt.from_ilong(128).bytewidth() == 2, "bytewidth(128)")
    
    // 65535 (0xFFFF) -> 16 bits + 1 sign bit = 17 bits -> 3 bytes
    h.assert_true(MPInt.from_ilong(65535).bytewidth() == 3, "bytewidth(65535)")
    // 255 (0xFF) -> 8 bits + 1 sign bit = 9 bits -> 2 bytes
    h.assert_true(MPInt.from_ilong(255).bytewidth() == 2, "bytewidth(255)")
    h.assert_true(MPInt.from_ilong(-1).bytewidth() == 1, "bytewidth(-1)")

    // bitwidth edge cases
    h.assert_true(MPInt.from_ilong(0).bitwidth().is_zero(),                            "bitwidth(0) = 0")
    h.assert_true(MPInt.from_ilong(1).bitwidth()     == MPInt.from_ilong(2),           "bitwidth(1) = 2")
    h.assert_true(MPInt.from_ilong(-1).bitwidth()    == MPInt.from_ilong(2),           "bitwidth(-1) = 2")
    h.assert_true(MPInt.from_ilong(127).bitwidth()   == MPInt.from_ilong(8),           "bitwidth(127) = 8")
    h.assert_true(MPInt.from_ilong(128).bitwidth()   == MPInt.from_ilong(9),           "bitwidth(128) = 9")
    h.assert_true(MPInt.from_ilong(255).bitwidth()   == MPInt.from_ilong(9),           "bitwidth(255) = 9")
    h.assert_true(MPInt.from_ilong(65535).bitwidth() == MPInt.from_ilong(17),          "bitwidth(65535) = 17")
    h.assert_true(MPInt.from_ilong(65536).bitwidth() == MPInt.from_ilong(18),          "bitwidth(65536) = 18")
    h.assert_true(MPInt.from_ilong(-127).bitwidth()  == MPInt.from_ilong(127).bitwidth(), "bitwidth(-x) = bitwidth(x)")
