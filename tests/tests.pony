"""
Tests for math package.
"""

use "collections"
use "random"

use "../mathx"
use "../pony_testx"

actor Main is TestList
  new create(env: Env) =>
    PonyTest(env, this)

  fun tag tests(test: PonyTest) =>
    """
    The list of tests on `maths`.
    """
    // No tests on platform limits!

    // Primes
    test(_TestPrimeTest)
    test(_TestPrimeTestExtended)
    test(_TestNextPrime)
    test(_TestPrimeFactorsStatic)
    test(_TestPrimeFactorsRandom)
    test(_TestPrimeIterator)
    test(_TestCoprime)
    test(_TestProbablyPrime)

    // Modular arithmetic
    test(_TestGCDLCM)
    test(_TestModularInverse)

    // Complex numbers
    test(_TestComplexAdd[F32])
    test(_TestComplexSub[F32])
    test(_TestComplexMul[F32])
    test(_TestComplexDiv[F32])
    test(_TestComplexAbs[F32])
    test(_TestComplexSqrt[F32])
    test(_TestComplexPow[F32])
    test(_TestComplexMisc[F32])
    test(_TestComplexTrigo[F32])
    test(_TestComplexLog[F32])
    test(_TestComplexIdentities[F32])
    test(_TestComplexRandom[F32, U32])
    test(_TestComplexAdd[F64])
    test(_TestComplexSub[F64])
    test(_TestComplexMul[F64])
    test(_TestComplexDiv[F64])
    test(_TestComplexAbs[F64])
    test(_TestComplexSqrt[F64])
    test(_TestComplexPow[F64])
    test(_TestComplexMisc[F64])
    test(_TestComplexTrigo[F64])
    test(_TestComplexLog[F64])
    test(_TestComplexIdentities[F64])
    test(_TestComplexRandom[F64, U64])

    // Miscellaneous fonctions for unsigned
    test(_TestUnsignedComplementPow[U8])
    test(_TestUnsignedComplementLog2[U8])
    test(_TestUnsignedComplementHamming[U8])

    // FFT: Fast Fourier Transform
    test(_TestFFT[F64])
    test(_TestFFT[F32])
    test(_TestFFT2[F64])
    test(_TestFFT2[F32])
    test(_TestFFTComplex[F64])
    test(_TestFFTComplex[F32])
    test(_TestFFTReal[F64])
    test(_TestFFTReal[F32])
    test(_TestFFTConvolution[F64])
    test(_TestFFTConvolution[F32])
    test(_TestFFTBluestein[F64])
    test(_TestFFTBluestein[F32])

    // MPFloat numbers
    test(_TestMPFloatCreate)
    test(_TestMPFloatFromF64)
    test(_TestMPFloatFrom)
    test(_TestMPFloatFromMPFloat)
    test(_TestMPFloatConversionF64)
    test(_TestMPFloatConversionF32)
    test(_TestMPFloatConversionISize)
    test(_TestMPFloatConversionUSize)
    test(_TestMPFloatConversionILong)
    test(_TestMPFloatConversionULong)
    test(_TestMPFloatConversionString)
    test(_TestMPFloatConversionExactString)
    test(_TestMPFloatConstants)
    test(_TestMPFloatMisc)
    test(_TestMPFloatComparisons)
    test(_TestMPFloatArithmetic)
    test(_TestMPFloatRoots)
    test(_TestMPFloatAbs)
    test(_TestMPFloatExponential)
    test(_TestMPFloatTrigonometric)
    test(_TestMPFloatHyperbolic)
    test(_TestMPFloatInteger)
    test(_TestMPFloatPrecision)

    // MPInt numbers
    test(_TestMPIntCreate)
    test(_TestMPIntFromILong)
    test(_TestMPIntComparisons)
    test(_TestMPIntAbsComparisons)
    test(_TestMPIntConversionILong)
    test(_TestMPIntMiscellaneous)
    test(_TestMPIntArithmetic)
    test(_TestMPIntKaratsuba)
    test(_TestMPIntFastMultiplication)

