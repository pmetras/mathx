"""
Tests for bignum/gmp package (requires libgmp and libmpfr).

Compiled separately from tests_bignum/ to avoid name clash between
bignum.MPFloat and bignum/gmp.MPFloat.
"""

use "../bignum/gmp"
use "../pony_testx"

actor Main is TestList
  new create(env: Env) =>
    PonyTest(env, this)

  fun tag tests(test: PonyTest) =>
    // BigReal[T] generic interface tests against the GMP/MPFR backend
    test(_TestBigRealGMPFloatPredicates)
    test(_TestBigRealGMPFloatComparisons)
    test(_TestBigRealGMPFloatArithmetic)
    test(_TestBigRealGMPFloatRoots)
    test(_TestBigRealGMPFloatMinMax)
    test(_TestBigRealGMPFloatRounding)
    test(_TestBigRealGMPFloatLogExp)
    test(_TestBigRealGMPFloatTrig)
    test(_TestBigRealGMPFloatInverseTrig)
    test(_TestBigRealGMPFloatHyp)
    test(_TestBigRealGMPFloatConversions)
