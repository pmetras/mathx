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
    h.assert_true(UnsignedComp[U].pow(U.from[USize](1), U.from[USize](12)) == U.from[USize](1))
    h.assert_true(UnsignedComp[U].pow(U.from[USize](12), U.from[USize](0)) == U.from[USize](1))
    h.assert_true(UnsignedComp[U].pow(U.from[USize](12), U.from[USize](1)) == U.from[USize](12))
    h.assert_true(UnsignedComp[U].pow(U.from[USize](12), U.from[USize](2)) == U.from[USize](144))

    let eight = U.from[USize](8)
    if eight.bitwidth() > eight then
      h.assert_true(UnsignedComp[U].pow(U.from[USize](12), U.from[USize](3)) == U.from[USize](1860867))
    end


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

    h.assert_true(UnsignedComp[U].log2(U.from[USize](1)) == U.from[USize](0))
    h.assert_true(UnsignedComp[U].log2(U.from[USize](123)) == U.from[USize](5))

    let eight = U.from[USize](8)
    let bitwidth = eight.bitwidth()
    if bitwidth > eight then
      h.assert_true(UnsignedComp[U].log2(U.from[USize](12345)) == U.from[USize](12))
    end
    if bitwidth > (eight + eight) then
      h.assert_true(UnsignedComp[U].log2(U.from[USize](12345678)) == U.from[USize](22))
    end

    // TODO: Can't manage to write the test for in in [0, bitwidth)
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
    h.assert_true(UnsignedComp[U].hamming_distance(U.from[USize](0), U.from[USize](123)) == U.from[USize](6))
    h.assert_true(UnsignedComp[U].hamming_distance(U.from[USize](1), U.from[USize](123)) == U.from[USize](5))
    h.assert_true(UnsignedComp[U].hamming_distance(U.from[USize](1), U.from[USize](255)) == U.from[USize](7))
    h.assert_true(UnsignedComp[U].hamming_distance(U.from[USize](170), U.from[USize](85)) == U.from[USize](8))


