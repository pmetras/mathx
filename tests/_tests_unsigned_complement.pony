// Tests for complimentary fonctions for unsigned integers


use "collections"
use "random"

use "../mathx"
use "../pony_testx"


class iso _TestUnsignedComplementPow[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on UnsignedComplement pow.
  """

  fun name(): String =>
    "UnsignedComplement[U]/pow"

  fun apply(h: TestHelper) =>
    let zero = U.from[USize](0)
    let one = U.from[USize](1)

    // 0^0 = 1 by convention (loop never runs, returns initial res = 1)
    h.assert_true(UnsignedComp[U].pow(zero, zero) == one)
    // 0^n = 0 for n > 0
    h.assert_true(UnsignedComp[U].pow(zero, one) == zero)
    h.assert_true(UnsignedComp[U].pow(zero, U.from[USize](5)) == zero)
    // 1^n = 1 for any n
    h.assert_true(UnsignedComp[U].pow(one, U.from[USize](12)) == one)
    // a^0 = 1
    h.assert_true(UnsignedComp[U].pow(U.from[USize](12), zero) == one)
    // a^1 = a
    h.assert_true(UnsignedComp[U].pow(U.from[USize](12), one) == U.from[USize](12))
    h.assert_true(UnsignedComp[U].pow(U.from[USize](12), U.from[USize](2)) == U.from[USize](144))

    let eight = U.from[USize](8)
    if eight.bitwidth() > eight then
      h.assert_true(UnsignedComp[U].pow(U.from[USize](12), U.from[USize](3)) == U.from[USize](1728))
    end


class iso _TestUnsignedComplementPowPartial[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on UnsignedComplement pow_partial (raises on overflow).
  """

  fun name(): String =>
    "UnsignedComplement[U]/pow_partial"

  fun apply(h: TestHelper) =>
    let zero = U.from[USize](0)
    let one = U.from[USize](1)

    // Basic non-overflowing cases
    h.assert_true(
      try UnsignedComp[U].pow_partial(zero, zero)? == one else false end)
    h.assert_true(
      try UnsignedComp[U].pow_partial(U.from[USize](12), one)? == U.from[USize](12) else false end)
    h.assert_true(
      try UnsignedComp[U].pow_partial(U.from[USize](12), U.from[USize](2))? == U.from[USize](144) else false end)

    // 2 ^ bitwidth overflows every unsigned type
    h.assert_true(
      try
        UnsignedComp[U].pow_partial(U.from[USize](2), zero.bitwidth())?
        false  // should not reach here
      else
        true   // error raised as expected
      end)


class iso _TestUnsignedComplementLog2[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on pow2 and log2
  """

  fun name(): String =>
    "UnsignedComplement[U]/pow2 & log2"

  fun apply(h: TestHelper) =>
    h.assert_true(UnsignedComp[U].pow2(U.from[USize](0)) == U.from[USize](1))
    h.assert_true(UnsignedComp[U].pow2(U.from[USize](1)) == U.from[USize](2))
    h.assert_true(UnsignedComp[U].pow2(U.from[USize](7)) == U.from[USize](128))
    h.assert_true(UnsignedComp[U].pow2(U.from[USize](6)) == U.from[USize](64))

    let zero = U.from[USize](0)

    // log2 is defined as popcount(a - 1); equals floor(log2) only for exact powers of 2
    h.assert_true(UnsignedComp[U].log2(U.from[USize](1)) == zero)
    h.assert_true(UnsignedComp[U].log2(U.from[USize](2)) == U.from[USize](1))
    h.assert_true(UnsignedComp[U].log2(U.from[USize](4)) == U.from[USize](2))
    h.assert_true(UnsignedComp[U].log2(U.from[USize](128)) == U.from[USize](7))
    // log2(0) = bitwidth (0 - 1 wraps to U.max_value() which has bitwidth set bits)
    h.assert_true(UnsignedComp[U].log2(zero) == zero.bitwidth())
    // non-powers of 2: popcount(a - 1)
    h.assert_true(UnsignedComp[U].log2(U.from[USize](123)) == U.from[USize](5))

    let eight = U.from[USize](8)
    let bitwidth = eight.bitwidth()
    if bitwidth > eight then
      // popcount(12344) = popcount(0b11000000111000) = 5
      h.assert_true(UnsignedComp[U].log2(U.from[USize](12345)) == U.from[USize](5))
    end
    if bitwidth > (eight + eight) then
      // popcount(12345677) = popcount(0xBC614D) = 12
      h.assert_true(UnsignedComp[U].log2(U.from[USize](12345678)) == U.from[USize](12))
    end

    // round-trip: log2(pow2(i)) == i for i in [0, 7]
    for i in Range(0, 8) do
      h.assert_true(UnsignedComp[U].log2(UnsignedComp[U].pow2(U.from[USize](i))) == U.from[USize](i))
    end
    // pow2 consistency with pow
    for i in Range(0, 8) do
      h.assert_true(UnsignedComp[U].pow2(U.from[USize](i)) == UnsignedComp[U].pow(U.from[USize](2), U.from[USize](i)))
    end


class iso _TestUnsignedComplementHamming[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on Hamming distance
  """

  fun name(): String =>
    "UnsignedComplement[U]/hamming_distance"

  fun apply(h: TestHelper) =>
    let zero = U.from[USize](0)

    // specific values (all fit in U8, so valid for any unsigned type)
    h.assert_true(UnsignedComp[U].hamming_distance(zero, U.from[USize](123)) == U.from[USize](6))
    h.assert_true(UnsignedComp[U].hamming_distance(U.from[USize](1), U.from[USize](123)) == U.from[USize](5))
    h.assert_true(UnsignedComp[U].hamming_distance(U.from[USize](1), U.from[USize](255)) == U.from[USize](7))
    // 170 = 0b10101010, 85 = 0b01010101: all 8 low bits differ
    h.assert_true(UnsignedComp[U].hamming_distance(U.from[USize](170), U.from[USize](85)) == U.from[USize](8))
    // reflexive: distance from a to itself is 0
    h.assert_true(UnsignedComp[U].hamming_distance(zero, zero) == zero)
    h.assert_true(UnsignedComp[U].hamming_distance(U.from[USize](42), U.from[USize](42)) == zero)
    // all bits differ: distance between 0 and U.max_value() equals the bitwidth
    h.assert_true(UnsignedComp[U].hamming_distance(zero, U.max_value()) == zero.bitwidth())
    // symmetry
    h.assert_true(
      UnsignedComp[U].hamming_distance(U.from[USize](170), U.from[USize](85)) ==
      UnsignedComp[U].hamming_distance(U.from[USize](85), U.from[USize](170)))


class iso _TestUnsignedComplementFloorLog2[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on floor_log2: the true floor(log2(a)), position of the highest set bit.
  """

  fun name(): String =>
    "UnsignedComplement[U]/floor_log2"

  fun apply(h: TestHelper) =>
    let zero = U.from[USize](0)
    let one = U.from[USize](1)

    // Exact powers of 2
    h.assert_true(UnsignedComp[U].floor_log2(one) == zero)
    h.assert_true(UnsignedComp[U].floor_log2(U.from[USize](2)) == one)
    h.assert_true(UnsignedComp[U].floor_log2(U.from[USize](4)) == U.from[USize](2))
    h.assert_true(UnsignedComp[U].floor_log2(U.from[USize](128)) == U.from[USize](7))
    // Non-powers of 2: floor(log2) equals position of highest bit
    h.assert_true(UnsignedComp[U].floor_log2(U.from[USize](3)) == one)    // 2^1=2 <= 3 < 4=2^2
    h.assert_true(UnsignedComp[U].floor_log2(U.from[USize](5)) == U.from[USize](2))  // 4 <= 5 < 8
    h.assert_true(UnsignedComp[U].floor_log2(U.from[USize](123)) == U.from[USize](6)) // 64 <= 123 < 128
    h.assert_true(UnsignedComp[U].floor_log2(U.from[USize](255)) == U.from[USize](7)) // 128 <= 255 < 256
    // Sentinel for 0
    h.assert_true(UnsignedComp[U].floor_log2(zero) == zero.bitwidth())
    // round-trip: floor_log2(2^k) == k for k in [0, 7]
    for k in Range(0, 8) do
      h.assert_true(UnsignedComp[U].floor_log2(UnsignedComp[U].pow2(U.from[USize](k))) == U.from[USize](k))
    end
    // pow2(floor_log2(a)) <= a for a > 0
    for a in Range(1, 256) do
      let fl = UnsignedComp[U].floor_log2(U.from[USize](a))
      h.assert_true(UnsignedComp[U].pow2(fl) <= U.from[USize](a))
    end


class iso _TestUnsignedComplementIsPowerOfTwo[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on is_power_of_two.
  """

  fun name(): String =>
    "UnsignedComplement[U]/is_power_of_two"

  fun apply(h: TestHelper) =>
    // 0 is not a power of 2
    h.assert_false(UnsignedComp[U].is_power_of_two(U.from[USize](0)))
    // exact powers of 2
    h.assert_true(UnsignedComp[U].is_power_of_two(U.from[USize](1)))
    h.assert_true(UnsignedComp[U].is_power_of_two(U.from[USize](2)))
    h.assert_true(UnsignedComp[U].is_power_of_two(U.from[USize](4)))
    h.assert_true(UnsignedComp[U].is_power_of_two(U.from[USize](64)))
    h.assert_true(UnsignedComp[U].is_power_of_two(U.from[USize](128)))
    // non-powers
    h.assert_false(UnsignedComp[U].is_power_of_two(U.from[USize](3)))
    h.assert_false(UnsignedComp[U].is_power_of_two(U.from[USize](5)))
    h.assert_false(UnsignedComp[U].is_power_of_two(U.from[USize](6)))
    h.assert_false(UnsignedComp[U].is_power_of_two(U.from[USize](123)))
    h.assert_false(UnsignedComp[U].is_power_of_two(U.from[USize](255)))
    // consistency: pow2(k) is always a power of 2
    for k in Range(0, 8) do
      h.assert_true(UnsignedComp[U].is_power_of_two(UnsignedComp[U].pow2(U.from[USize](k))))
    end


class iso _TestUnsignedComplementNextPow2[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on next_power_of_two.
  """

  fun name(): String =>
    "UnsignedComplement[U]/next_power_of_two"

  fun apply(h: TestHelper) =>
    let zero = U.from[USize](0)
    let one = U.from[USize](1)

    // Boundary: a <= 1 returns 1
    h.assert_true(UnsignedComp[U].next_power_of_two(zero) == one)
    h.assert_true(UnsignedComp[U].next_power_of_two(one) == one)
    // Exact powers of 2 map to themselves
    h.assert_true(UnsignedComp[U].next_power_of_two(U.from[USize](2)) == U.from[USize](2))
    h.assert_true(UnsignedComp[U].next_power_of_two(U.from[USize](4)) == U.from[USize](4))
    h.assert_true(UnsignedComp[U].next_power_of_two(U.from[USize](128)) == U.from[USize](128))
    // Non-powers round up
    h.assert_true(UnsignedComp[U].next_power_of_two(U.from[USize](3)) == U.from[USize](4))
    h.assert_true(UnsignedComp[U].next_power_of_two(U.from[USize](5)) == U.from[USize](8))
    h.assert_true(UnsignedComp[U].next_power_of_two(U.from[USize](100)) == U.from[USize](128))
    // result is always a power of 2 (for non-overflowing inputs)
    for a in Range(1, 129) do
      h.assert_true(UnsignedComp[U].is_power_of_two(UnsignedComp[U].next_power_of_two(U.from[USize](a))))
    end
    // result >= a
    for a in Range(1, 129) do
      h.assert_true(UnsignedComp[U].next_power_of_two(U.from[USize](a)) >= U.from[USize](a))
    end


class iso _TestUnsignedComplementIsqrt[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on isqrt: integer square root.
  """

  fun name(): String =>
    "UnsignedComplement[U]/isqrt"

  fun apply(h: TestHelper) =>
    let zero = U.from[USize](0)
    let one = U.from[USize](1)

    // Base cases
    h.assert_true(UnsignedComp[U].isqrt(zero) == zero)
    h.assert_true(UnsignedComp[U].isqrt(one) == one)
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](2)) == one)
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](3)) == one)
    // Perfect squares
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](4)) == U.from[USize](2))
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](9)) == U.from[USize](3))
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](16)) == U.from[USize](4))
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](25)) == U.from[USize](5))
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](144)) == U.from[USize](12))
    // Non-perfect squares: floor(sqrt)
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](2)) == one)   // sqrt(2)≈1.41
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](10)) == U.from[USize](3))  // sqrt(10)≈3.16
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](15)) == U.from[USize](3))  // sqrt(15)≈3.87
    h.assert_true(UnsignedComp[U].isqrt(U.from[USize](255)) == U.from[USize](15)) // sqrt(255)≈15.97
    // Defining property: k = isqrt(a)  =>  k^2 <= a < (k+1)^2
    // Upper bound up to 224: isqrt(224)=14 and 15^2=225 fits in U8 (no overflow)
    for a in Range(0, 225) do
      let k = UnsignedComp[U].isqrt(U.from[USize](a))
      h.assert_true((k * k) <= U.from[USize](a))
      h.assert_true(U.from[USize](a) < ((k + one) * (k + one)))
    end


class iso _TestUnsignedComplementParity[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on parity: odd/even number of set bits.
  """

  fun name(): String =>
    "UnsignedComplement[U]/parity"

  fun apply(h: TestHelper) =>
    // Even parity (even number of set bits)
    h.assert_false(UnsignedComp[U].parity(U.from[USize](0)))   // 0 set bits
    h.assert_false(UnsignedComp[U].parity(U.from[USize](3)))   // 0b11: 2 set bits
    h.assert_false(UnsignedComp[U].parity(U.from[USize](15)))  // 0b1111: 4 set bits
    h.assert_false(UnsignedComp[U].parity(U.from[USize](255))) // 8 set bits
    // Odd parity (odd number of set bits)
    h.assert_true(UnsignedComp[U].parity(U.from[USize](1)))   // 1 set bit
    h.assert_true(UnsignedComp[U].parity(U.from[USize](2)))   // 1 set bit
    h.assert_true(UnsignedComp[U].parity(U.from[USize](4)))   // 1 set bit
    h.assert_true(UnsignedComp[U].parity(U.from[USize](7)))   // 0b111: 3 set bits
    h.assert_true(UnsignedComp[U].parity(U.from[USize](127))) // 7 set bits
    // Consistency: parity(a xor b) == parity(a) xor parity(b)
    let a = U.from[USize](85)   // 0b01010101
    let b = U.from[USize](170)  // 0b10101010
    h.assert_true(UnsignedComp[U].parity(a xor b) == (UnsignedComp[U].parity(a) xor UnsignedComp[U].parity(b)))


class iso _TestUnsignedComplementIlog[U: UnsignedInteger[U] val] is UnitTest
  """
  Tests on ilog: floor of log in an arbitrary base.
  """

  fun name(): String =>
    "UnsignedComplement[U]/ilog"

  fun apply(h: TestHelper) =>
    let zero = U.from[USize](0)
    let one = U.from[USize](1)

    // Base 2: should match floor_log2
    h.assert_true(UnsignedComp[U].ilog(one, U.from[USize](2)) == zero)
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](2), U.from[USize](2)) == one)
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](7), U.from[USize](2)) == U.from[USize](2))
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](8), U.from[USize](2)) == U.from[USize](3))
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](255), U.from[USize](2)) == U.from[USize](7))
    // Base 10
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](9), U.from[USize](10)) == zero)
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](10), U.from[USize](10)) == one)
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](99), U.from[USize](10)) == one)
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](100), U.from[USize](10)) == U.from[USize](2))
    // Base 3
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](8), U.from[USize](3)) == one)   // 3^1=3<=8<9=3^2
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](9), U.from[USize](3)) == U.from[USize](2))
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](26), U.from[USize](3)) == U.from[USize](2)) // 9<=26<27
    h.assert_true(UnsignedComp[U].ilog(U.from[USize](27), U.from[USize](3)) == U.from[USize](3))
    // Undefined cases return 0
    h.assert_true(UnsignedComp[U].ilog(zero, U.from[USize](2)) == zero)
    h.assert_true(UnsignedComp[U].ilog(one, U.from[USize](0)) == zero)
    h.assert_true(UnsignedComp[U].ilog(one, one) == zero)
    // Consistency with floor_log2 for base 2
    for a in Range(1, 256) do
      h.assert_true(
        UnsignedComp[U].ilog(U.from[USize](a), U.from[USize](2)) ==
        UnsignedComp[U].floor_log2(U.from[USize](a)))
    end


