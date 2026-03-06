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


