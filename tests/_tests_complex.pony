// Tests for complex functions

// Equality with floats is tricky. These tests were written to allow `F32` pass.
// When strict equality `==` was not possible, it was replaced with a near
// equality with `almost_eq()` and a tolerance factor. As the platform use
// IEEE 754  representation for float, it allowed to detect where numeric errors
// occur. This was done to spot where some operations create rounding errors.
// Then they were tested with `F64`, and new `==` tests were found failing...
//
// Conclusion, be very cautious when using `==` with floats.
// Also some tests could fail on other hardware or platforms not using IEEE 754
// floats.

use "collections"
use "random"

use "../mathx"
use "../pony_testx"


class iso _TestComplexAdd[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Additions on Complex[F].
  """
  fun name(): String =>
    "Complex[F]/Additions"

  fun apply(h: TestHelper) =>
    let a = Complex[F](F.from[F64](-2.0), F.from[F64](1.0))
    let b = Complex[F](F.from[F64](1.0), F.from[F64](-2.0))
    let c = Complex[F](F.from[F64](-1.0), F.from[F64](-1.0))

    h.assert_true((a + b) == c)
    h.assert_true((a +~ b) ==~ c)

    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    let four = F.from[ISize](4)
    let one_plus_i = Complex[F](one, one)
    let one_minus_i = Complex[F](one, -one)
    
    // Safe
    h.assert_true((one_plus_i * one_minus_i) == Complex[F](two, zero))
    h.assert_true((one_plus_i * one_plus_i) == (Complex[F](two, zero) * Complex[F].j()))
    h.assert_true((one_minus_i * one_minus_i) == (Complex[F](-two, zero) * Complex[F].j()))
    h.assert_true((one_plus_i * one_plus_i * one_plus_i) == Complex[F](-two, two))
    h.assert_true((one_minus_i * one_minus_i * one_minus_i) == (Complex[F](-two, -two)))
    h.assert_true((one_plus_i / one_minus_i) == Complex[F].j())
    h.assert_true((one_plus_i.powi(2) / one_minus_i) == Complex[F](-one, one))
    h.assert_true((one_plus_i.powi(2) / one_minus_i.powi(2)) == Complex[F](-one, zero))
    h.assert_true((one_plus_i.powi(2) + one_minus_i.powi(2)) == Complex[F](zero, zero))
    h.assert_true((one_plus_i.powi(3) + one_minus_i.powi(3)) == Complex[F](-four, zero))
    
    // Unsafe
    h.assert_true((one_plus_i *~ one_minus_i) ==~ Complex[F](two, zero))
    h.assert_true((one_plus_i *~ one_plus_i) ==~ (Complex[F](two, zero) *~ Complex[F].j()))
    h.assert_true((one_minus_i *~ one_minus_i) ==~ (Complex[F](-~two, zero) *~ Complex[F].j()))
    h.assert_true((one_plus_i *~ one_plus_i *~ one_plus_i) ==~ Complex[F](-~two, two))
    h.assert_true((one_minus_i *~ one_minus_i *~ one_minus_i) ==~ (Complex[F](-~two, -~two)))
    h.assert_true((one_plus_i /~ one_minus_i) ==~ Complex[F].j())
    h.assert_true((one_plus_i.powi_unsafe(2) /~ one_minus_i) ==~ Complex[F](-~one, one))
    h.assert_true((one_plus_i.powi_unsafe(2) /~ one_minus_i.powi_unsafe(2)) ==~ Complex[F](-~one, zero))
    h.assert_true((one_plus_i.powi_unsafe(2) +~ one_minus_i.powi_unsafe(2)) ==~ Complex[F](zero, zero))
    h.assert_true((one_plus_i.powi_unsafe(3) +~ one_minus_i.powi_unsafe(3)) ==~ Complex[F](-~four, zero))


class iso _TestComplexSub[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Substractions on Complex[F].
  """
  fun name(): String =>
    "Complex[F]/Substractions"

  fun apply(h: TestHelper) =>
    let a = Complex[F](F.from[F64](-2.0), F.from[F64](1.0))
    let b = Complex[F](F.from[F64](1.0), F.from[F64](-2.0))
    let c = Complex[F](F.from[F64](-3.0), F.from[F64](3.0))

    h.assert_true((a - b) == c)
    h.assert_true((a -~ b) ==~ c)


class iso _TestComplexMul[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Multiplications on Complex[F].
  """
  fun name(): String =>
    "Complex[F]/Multiplications"

  fun apply(h: TestHelper) =>
    let a = Complex[F](F.from[F64](-2.0), F.from[F64](1.0))
    let b = Complex[F](F.from[F64](1.0), F.from[F64](-2.0))
    let c = Complex[F](F.from[F64](0.0), F.from[F64](5.0))

    h.assert_true((a * b) == c)
    h.assert_true((a *~ b) ==~ c)


class iso _TestComplexDiv[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Divisions on Complex[F].
  """
  fun name(): String =>
    "Complex[F]/Divisions"

  fun apply(h: TestHelper) =>
    let a = Complex[F](F.from[F64](-2.0), F.from[F64](1.0))
    let b = Complex[F](F.from[F64](1.0), F.from[F64](-2.0))
    let c = Complex[F](F.from[F64](-4.0) / F.from[F64](5.0), F.from[F64](-3.0) / F.from[F64](5.0))
    let d = Complex[F](F.from[F64](-5.9), F.from[F64](0.0))
    let e = Complex[F](F.from[F64](0.0), F.from[F64](-8.3))

    let one = Complex[F](F.from[F64](1.0), F.from[F64](0.0))

    // Division
    h.assert_true((a / b) == c)
    h.assert_true((a /~ b) ==~ c)

    // Inverse
    // Safe
    h.log("a.invert() = " + a.invert().string() + ", 1 / a = " + (one / a).string())
    h.assert_true(a.invert() == (one / a))
    h.log("b.invert() = " + b.invert().string() + ", 1 / b = " + (one / b).string())
    h.assert_true(b.invert() == (one / b))
    h.log("c.invert() = " + c.invert().string() + ", 1 / c = " + (one / c).string())
    h.assert_true(c.invert() == (one / c))
    h.log("d.invert() = " + d.invert().string() + ", 1 / d = " + (one / d).string())
    h.assert_true(d.invert() == (one / d))
    h.log("e.invert() = " + e.invert().string() + ", 1 / e = " + (one / e).string())
    h.assert_true(e.invert() == (one / e))
    
    // Unsafe
    h.log("a.invert_unsafe() = " + a.invert_unsafe().string() + ", 1 /~ a = " + (one /~ a).string())
    h.assert_true(a.invert_unsafe() ==~ (one /~ a))
    h.log("b.invert() = " + b.invert_unsafe().string() + ", 1 /~ b = " + (one /~ b).string())
    h.assert_true(b.invert_unsafe() ==~ (one /~ b))
    h.log("c.invert() = " + c.invert_unsafe().string() + ", 1 /~ c = " + (one /~ c).string())
    h.assert_true(c.invert_unsafe() ==~ (one /~ c))
    h.log("d.invert() = " + d.invert_unsafe().string() + ", 1 /~ d = " + (one /~ d).string())
    h.assert_true(d.invert_unsafe() ==~ (one /~ d))
    h.log("e.invert() = " + e.invert_unsafe().string() + ", 1 /~ e = " + (one /~ e).string())
    h.assert_true(e.invert_unsafe() ==~ (one /~ e))


class iso _TestComplexAbs[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Modulo of Complex[F].
  """
  fun name(): String =>
    "Complex[F]/Modulus"

  fun apply(h: TestHelper) =>
    let a = Complex[F](F.from[F64](-2.0), F.from[F64](1.0))
    let b = Complex[F](F.from[F64](3.0), F.from[F64](-4.0))
    let c = Complex[F](F.from[F64](-1.0), F.from[F64](0.0))
    let d = Complex[F](F.from[F64](0.0), F.from[F64](-1.0))

    // Safe
    h.assert_true(a.abs() == F.from[F64](5.0).sqrt().abs())
    h.assert_true(b.abs() == F.from[F64](5.0))
    h.assert_true(c.abs() == F.from[F64](1.0))
    h.assert_true(d.abs() == F.from[F64](1.0))
    h.assert_true(d.abs() == c.abs())

    // Unsafe
    h.assert_true(a.abs_unsafe() ==~ F.from[F64](5.0).sqrt_unsafe().abs())
    h.assert_true(b.abs_unsafe() ==~ F.from[F64](5.0))
    h.assert_true(c.abs_unsafe() ==~ F.from[F64](1.0))
    h.assert_true(d.abs_unsafe() ==~ F.from[F64](1.0))
    h.assert_true(d.abs_unsafe() ==~ c.abs_unsafe())


class iso _TestComplexSqrt[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Square root on Complex[F].
  """
  fun name(): String =>
    "Complex[F]/Square root"

  fun apply(h: TestHelper) =>
    let a = Complex[F](F.from[F64](3.0), F.from[F64](4.0))
    let b = Complex[F](F.from[F64](-2.0), F.from[F64](-1.0))
    let c = Complex[F](F.from[F64](2.0), F.from[F64](1.0))
    let d = Complex[F](F.from[F64](8.0), F.from[F64](-6.0))
    let e = Complex[F](F.from[F64](3.0), F.from[F64](-1.0))
    let f = Complex[F](F.from[F64](-3.0), F.from[F64](1.0))

    // Safe
    h.log("a = " + a.string())
    h.log("a.sqrt() = " + a.sqrt().string())
    h.assert_true((a.sqrt() == b) or (a.sqrt() == c))
    h.log("b = " + b.string() + ", b * b = " + (b * b).string())
    h.assert_true((b * b) == a)
    h.log("c = " + c.string() + ", c * c = " + (c * c).string())
    h.assert_true((c * c) == a)

    h.log("d = " + d.string())
    h.log("d.sqrt() = " + d.sqrt().string())
    // Never trust equality with floats. This can be seen
    // when square root of d is printed "3-1i" instead of "3-i"
    h.assert_true(d.sqrt().almost_eq(e) or d.sqrt().almost_eq(f))
    h.log("e = " + e.string() + ", e * e = " + (e * e).string())
    h.assert_true((e * e) == d)
    h.log("f = " + f.string() + ", f * f = " + (f * f).string())
    h.assert_true((f * f) == d)

    // Unsafe
    h.log("a = " + a.string())
    h.log("a.sqrt_unsafe() = " + a.sqrt_unsafe().string())
    h.assert_true((a.sqrt_unsafe() ==~ b) or (a.sqrt_unsafe() ==~ c))
    h.log("b = " + b.string() + ", b *~ b = " + (b *~ b).string())
    h.assert_true((b *~ b) ==~ a)
    h.log("c = " + c.string() + ", c *~ c = " + (c *~ c).string())
    h.assert_true((c *~ c) ==~ a)

    h.log("d = " + d.string())
    h.log("d.sqrt_unsafe() = " + d.sqrt_unsafe().string())
    // Never trust equality with floats. This can be seen
    // when square root of d is printed "3-1i" instead of "3-i"
    h.assert_true(d.sqrt_unsafe().almost_eq(e) or
                  d.sqrt_unsafe().almost_eq(f))
    h.log("e = " + e.string() + ", e *~ e = " + (e *~ e).string())
    h.assert_true((e *~ e) ==~ d)
    h.log("f = " + f.string() + ", f *~ f = " + (f *~ f).string())
    h.assert_true((f *~ f) ==~ d)


class iso _TestComplexPow[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Power of complex
  """
  fun name(): String =>
    "Complex[F]/Power"

  fun apply(h: TestHelper) =>
    let a = Complex[F](F.from[F64](-2.8), F.from[F64](1.9))
    let b = Complex[F](F.from[F64](1.1), F.from[F64](-2.5))
    let c = Complex[F](F.from[F64](-6.1), F.from[F64](0.0))
    let d = Complex[F](F.from[F64](0.0), F.from[F64](5.3))

    let one = Complex[F](F.from[F64](1.0), F.from[F64](0.0))

    // Safe
    h.log("a = " + a.string() + ", b = " + b.string() + ", c = " + c.string() +
          ", d = " + d.string())
    h.log("a.powi(2) = " + a.powi(2).string() + ", a * a = " + (a * a).string() +
          ", delta = " + (a.powi(2) - (a * a)).string())
    h.assert_true(a.powi(2) == (a * a))
    h.log("a.powi(-5) = " + a.powi(-5).string() + ", 1 / a.powi(5) = " +
          (one / a.powi(5)).string() + ", delta = " +
          (a.powi(-5) - (one / a.powi(5))).string())
    h.assert_true(a.powi(-5) == (one / a.powi(5)))
    h.log("a.powi(7) = " + a.powi(7).string() + ", a *..7..* a = " +
          (a * a * a * a * a * a * a).string() + ", delta = " +
          (a.powi(7) - (a * a * a * a * a * a * a)).string())
    h.assert_true(a.powi(7).almost_eq(a * a * a * a * a * a * a))
    h.log("b.powi(2) = " + b.powi(2).string() + ", b * b = " + (b * b).string() +
          ", delta = " + (b.powi(2) - (b * b)).string())
    // This multiplication is not accurate with F64
    h.assert_true(b.powi(2).almost_eq(b * b))
    h.log("b.powi(-5) = " + b.powi(-5).string() + ", 1 / b.powi(5) = " +
          (one / b.powi(5)).string() + ", delta = " +
          (b.powi(-5) - (one / b.powi(5))).string())
    h.assert_true(b.powi(-5) == (one / b.powi(5)))
    h.log("b.powi(7) = " + b.powi(7).string() + ", b *..7..* b = " +
          (b * b * b * b * b * b * b).string() + ", delta = " +
          (b.powi(7) - (b * b * b * b * b * b * b)).string())
    h.assert_true(b.powi(7).almost_eq(b * b * b * b * b * b * b))
    h.log("c.powi(2) = " + c.powi(2).string() + ", c * c = " + (c * c).string() +
          ", delta = " + (c.powi(2) - (c * c)).string())
    h.assert_true(c.powi(2) == (c * c))
    h.log("c.powi(-5) = " + c.powi(-5).string() + ", 1 / c.powi(5) = " +
          (one / c.powi(5)).string() + ", delta = " +
          (c.powi(-5) - (one / c.powi(5))).string())
    h.assert_true(c.powi(-5) == (one / c.powi(5)))
    h.log("c.powi(7) = " + c.powi(7).string() + ", c *..7..* c = " +
          (c * c * c * c * c * c * c).string() + ", delta = " +
          (c.powi(7) - (c * c * c * c * c * c * c)).string())
    h.assert_true(c.powi(7).almost_eq(c * c * c * c * c * c * c))
    h.log("d.powi(2) = " + d.powi(2).string() + ", d * d = " + (d * d).string() +
          ", delta = " + (d.powi(2) - (d * d)).string())
    h.assert_true(d.powi(2) == (d * d))
    h.log("d.powi(-5) = " + d.powi(-5).string() + ", 1 / d.powi(5) = " +
          (one / d.powi(5)).string() + ", delta = " +
          (d.powi(-5) - (one / d.powi(5))).string())
    h.assert_true(d.powi(-5) == (one / d.powi(5)))
    h.log("d.powi(7) = " + d.powi(7).string() + ", d *..7..* d = " +
          (d * d * d * d * d * d * d).string() + ", delta = " +
          (d.powi(7) - (d * d * d * d * d * d * d)).string())
    h.assert_true(d.powi(7).almost_eq(d * d * d * d * d * d * d))

    // Unsafe
    h.log("a = " + a.string() + ", b = " + b.string() + ", c = " + c.string() +
          ", d = " + d.string())
    h.log("a.powi_unsafe(2) = " + a.powi_unsafe(2).string() + ", a *~ a = " +
          (a *~ a).string() + ", delta = " + (a.powi_unsafe(2) -~ (a *~ a)).string())
    h.assert_true(a.powi_unsafe(2) ==~ (a *~ a))
    h.log("a.powi_unsafe(-5) = " + a.powi_unsafe(-5).string() +
          ", 1 /~ a.powi_unsafe(5) = " + (one /~ a.powi_unsafe(5)).string() +
          ", delta = " + (a.powi_unsafe(-5) -~ (one /~ a.powi_unsafe(5))).string())
    h.assert_true(a.powi_unsafe(-5) ==~ (one /~ a.powi_unsafe(5)))
    h.log("a.powi_unsafe(7) = " + a.powi_unsafe(7).string() + ", a *~..7..*~ a = " +
          (a *~ a *~ a *~ a *~ a *~ a *~ a).string() + ", delta = " +
          (a.powi_unsafe(7) -~ (a *~ a *~ a *~ a *~ a *~ a *~ a)).string())
    h.assert_true(a.powi_unsafe(7).almost_eq(a *~ a *~ a *~ a *~ a *~ a *~ a))
    h.log("b.powi_unsafe(2) = " + b.powi_unsafe(2).string() + ", b *~ b = " +
          (b *~ b).string() + ", delta = " + (b.powi_unsafe(2) -~ (b *~ b)).string())
    // This multiplication is not accurate with F64
    h.assert_true(b.powi_unsafe(2).almost_eq(b *~ b))
    h.log("b.powi_unsafe(-5) = " + b.powi_unsafe(-5).string() +
          ", 1 /~ b.powi_unsafe(5) = " + (one /~ b.powi_unsafe(5)).string() +
          ", delta = " + (b.powi_unsafe(-5) -~ (one /~ b.powi_unsafe(5))).string())
    h.assert_true(b.powi_unsafe(-5) ==~ (one /~ b.powi_unsafe(5)))
    h.log("b.powi_unsafe(7) = " + b.powi_unsafe(7).string() + ", b *~..7..*~ b = " +
          (b *~ b *~ b *~ b *~ b *~ b *~ b).string() + ", delta = " +
          (b.powi_unsafe(7) -~ (b *~ b *~ b *~ b *~ b *~ b *~ b)).string())
    h.assert_true(b.powi_unsafe(7).almost_eq(b *~ b *~ b *~ b *~ b *~ b *~ b))
    h.log("c.powi_unsafe(2) = " + c.powi_unsafe(2).string() + ", c *~ c = " +
          (c *~ c).string() + ", delta = " + (c.powi_unsafe(2) -~ (c *~ c)).string())
    h.assert_true(c.powi_unsafe(2) ==~ (c *~ c))
    h.log("c.powi_unsafe(-5) = " + c.powi_unsafe(-5).string() +
          ", 1 /~ c.powi_unsafe(5) = " + (one /~ c.powi_unsafe(5)).string() +
          ", delta = " + (c.powi_unsafe(-5) -~ (one /~ c.powi_unsafe(5))).string())
    h.assert_true(c.powi_unsafe(-5) ==~ (one /~ c.powi_unsafe(5)))
    h.log("c.powi_unsafe(7) = " + c.powi_unsafe(7).string() + ", c *~..7..*~ c = " +
          (c *~ c *~ c *~ c *~ c *~ c *~ c).string() + ", delta = " +
          (c.powi_unsafe(7) -~ (c *~ c *~ c *~ c *~ c *~ c *~ c)).string())
    h.assert_true(c.powi_unsafe(7).almost_eq(c *~ c *~ c *~ c *~ c *~ c *~ c))
    h.log("d.powi_unsafe(2) = " + d.powi_unsafe(2).string() + ", d *~ d = " +
          (d *~ d).string() + ", delta = " + (d.powi_unsafe(2) -~ (d *~ d)).string())
    h.assert_true(d.powi_unsafe(2) ==~ (d *~ d))
    h.log("d.powi_unsafe(-5) = " + d.powi_unsafe(-5).string() +
          ", 1 /~ d.powi_unsafe(5) = " + (one /~ d.powi_unsafe(5)).string() +
          ", delta = " + (d.powi_unsafe(-5) -~ (one /~ d.powi_unsafe(5))).string())
    h.assert_true(d.powi_unsafe(-5) ==~ (one /~ d.powi_unsafe(5)))
    h.log("d.powi_unsafe(7) = " + d.powi_unsafe(7).string() + ", d *~..7..*~ d = " +
          (d *~ d *~ d *~ d *~ d *~ d *~ d).string() + ", delta = " +
          (d.powi_unsafe(7) -~ (d *~ d *~ d *~ d *~ d *~ d *~ d)).string())
    h.assert_true(d.powi_unsafe(7).almost_eq(d *~ d *~ d *~ d *~ d *~ d *~ d))


class iso _TestComplexMisc[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Miscellaneous operations on Complex[F].
  """
  fun name(): String =>
    "Complex[F]/Miscellaneous operation"

  fun apply(h: TestHelper) =>
    let a = Complex[F](F.from[F64](-2.0), F.from[F64](1.0))
    let b = Complex[F](F.from[F64](1.0), F.from[F64](-2.0))
    let c = Complex[F](F.from[F64](-4.0) / F.from[F64](5.0), F.from[F64](-3.0) / F.from[F64](5.0))
    let d = Complex[F](F.from[F64](163.0), F.from[F64](0.0))
    let e = Complex[F](F.from[F64](0.0), F.from[F64](-510.883))

    h.assert_true(a.real() == F.from[F64](-2.0))
    h.assert_true(a.imag() == F.from[F64](1.0))

    h.log("a = " + a.string())
    h.assert_true(a.string() == "-2+i")
    h.log("b = " + b.string())
    h.assert_true(b.string() == "1-2i")
    h.log("c = " + c.string())
    h.assert_true(c.string() == "-0.8-0.6i")
    h.log("d = " + d.string())
    h.assert_true(d.string() == "163")
    h.log("e = " + e.string())
    h.assert_true(e.string() == "-510.883i")
    h.log("i = " + Complex[F].j().string())
    h.assert_true(Complex[F].j().string() == "i")

    h.log("d.is_real() = " + d.is_real().string())
    h.assert_true(d.is_real())
    h.log("e.is_imag() = " + e.is_imag().string())
    h.assert_true(e.is_imag())
    h.log("((a / b) - c).is_real() = " + ((a / b) - c).is_real().string())
    h.assert_true(((a / b) - c).is_real())
    h.log("((a / b) - c).is_imag() = " + ((a / b) - c).is_imag().string())
    h.assert_true(((a / b) - c).is_imag())
    h.log("(a / b) - c = " + ((a / b) - c).string())
    h.assert_true(((a / b) - c) == Complex[F](F.from[F64](0.0), F.from[F64](0.0)))
    h.log("(a - a).is_null() = " + (a - a).is_null().string())
    h.assert_true((a - a).is_null())
    h.log("(a - b).is_null() = " + (a -b).is_null().string())
    h.assert_false((a - b).is_null())
    h.log("(a - c).finite() = " + (a - c).finite().string())
    h.assert_true((a - c).finite())

    let inf = F.from[ISize](1) / F.from[ISize](0)
    let nan = F.from[ISize](0) / F.from[ISize](0)

    h.log("(inf+i).infinite() = " + Complex[F](inf, F.from[ISize](1)).infinite().string())
    h.assert_true(Complex[F](inf, F.from[ISize](1)).infinite())
    h.log("(1+infi).finite() = " + Complex[F](F.from[ISize](1), inf).finite().string())
    h.assert_false(Complex[F](F.from[ISize](1), inf).finite())
    h.log("(nan+infi).finite() = " + Complex[F](nan, inf).finite().string())
    h.assert_false(Complex[F](nan, inf).finite())
    h.log("(nan+infi).infinite() = " + Complex[F](nan, inf).infinite().string())
    h.assert_true(Complex[F](nan, inf).infinite())
    h.log("(nan+infi).nan() = " + Complex[F](nan, inf).nan().string())
    h.assert_true(Complex[F](nan, inf).nan())


class iso _TestComplexTrigo[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Trigonometry of complex
  """
  fun name(): String =>
    "Complex[F]/Trigonometry"

  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let pi = F.from[F64](F64.pi())
    let one_r = Complex[F](one, zero)
    let zero_r = Complex[F](zero, zero)
    let i = Complex[F].j()

    // At the origin: cos(0)=1, sin(0)=0, cosh(0)=1, sinh(0)=0, tanh(0)=0
    h.assert_true(zero_r.cos().almost_eq(one_r))
    h.assert_true(zero_r.sin().almost_eq(zero_r))
    h.assert_true(zero_r.cosh().almost_eq(one_r))
    h.assert_true(zero_r.sinh().almost_eq(zero_r))
    h.assert_true(zero_r.tanh().almost_eq(zero_r))

    // cos(π/2) ≈ 0, sin(π/2) = 1
    let pi_2 = Complex[F](pi / F.from[ISize](2), zero)
    h.log("cos(π/2) = " + pi_2.cos().string())
    h.assert_true(pi_2.cos().almost_eq(zero_r))
    h.log("sin(π/2) = " + pi_2.sin().string())
    h.assert_true(pi_2.sin().almost_eq(one_r))

    // Cross-identities between trig and hyperbolic:
    // cos(iz) = cosh(z) → cos(i) = cosh(1)
    h.log("cos(i) = " + i.cos().string() + ", cosh(1) = " + one_r.cosh().string())
    h.assert_true(i.cos().almost_eq(one_r.cosh()))
    // sin(iz) = i*sinh(z) → sin(i) = i*sinh(1)
    let i_sinh1 = i * one_r.sinh()
    h.log("sin(i) = " + i.sin().string() + ", i * sinh(1) = " + i_sinh1.string())
    h.assert_true(i.sin().almost_eq(i_sinh1))
    // tan(iz) = i*tanh(z) → tan(i) = i*tanh(1)
    let i_tanh1 = i * one_r.tanh()
    h.log("tan(i) = " + i.tan().string() + ", i * tanh(1) = " + i_tanh1.string())
    h.assert_true(i.tan().almost_eq(i_tanh1))

    // Pythagorean identity: sin²(z) + cos²(z) = 1
    let z = Complex[F](F.from[F64](1.0), F.from[F64](0.5))
    let sz = z.sin()
    let cz = z.cos()
    let pyth = (sz * sz) + (cz * cz)
    h.log("sin²(z) + cos²(z) = " + pyth.string())
    h.assert_true(pyth.almost_eq(one_r))

    // Hyperbolic Pythagorean: cosh²(z) - sinh²(z) = 1
    let chz = z.cosh()
    let shz = z.sinh()
    let hyp_pyth = (chz * chz) - (shz * shz)
    h.log("cosh²(z) - sinh²(z) = " + hyp_pyth.string())
    h.assert_true(hyp_pyth.almost_eq(one_r))


class iso _TestComplexLog[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Logarithm of complex
  """
  fun name(): String =>
    "Complex[F]/Logarithm"

  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    let three = F.from[ISize](3)
    let four = F.from[ISize](4)
    let five = F.from[ISize](5)
    let pi = F.from[F64](F64.pi())

    h.log("log(1) = 0 = " + Complex[F](one, zero).log().string() + ", delta = " +
          (Complex[F](one, zero).log() - Complex[F](zero, zero)).string())
    h.assert_true(Complex[F](one, zero).log().is_null())
    h.log("log(1+i) = " + Complex[F](one, one).log().string() +
          ", log(2)/2 + i * pi/4 = " + Complex[F](two.log() / two, pi / four).string() +
          ", delta = " + (Complex[F](one, one).log() - Complex[F](two.log() / two, pi / four)).string())
    h.assert_true(Complex[F](one, one).log().almost_eq(Complex[F](two.log() / two, pi / four)))
    h.log("log(2+i) = " + Complex[F](two, one).log().string() +
          ", log(5)/2 + i * atan(1/2) = " + Complex[F](five.log() / two, (one / two).atan()).string() +
          ", delta = " + (Complex[F](two, one).log() - Complex[F](five.log() / two, (one / two).atan())).string())
    h.assert_true(Complex[F](two, one).log().almost_eq(Complex[F](five.log() / two, (one / two).atan())))
    h.log("log(3+i) = " + Complex[F](three, one).log().string() +
          ", log(10)/2 + i * atan(1/3) = " + Complex[F]((two * five).log() / two, (one / three).atan()).string() +
          ", delta = " + (Complex[F](three, one).log() - Complex[F]((two * five).log() / two, (one / three).atan())).string())
    h.assert_true(Complex[F](three, one).log() == Complex[F]((two * five).log() / two, (one / three).atan()))
    h.log("log(1+2i) = " + Complex[F](one, two).log().string() +
          ", log(5)/2 + i * atan(2) = " + Complex[F](five.log() / two, two.atan()).string() +
          ", delta = " + (Complex[F](one, two).log() - Complex[F](five.log() / two, two.atan())).string())
    h.assert_true(Complex[F](one, two).log().almost_eq(Complex[F](five.log() / two, two.atan())))
    h.log("log(i) = " + Complex[F].j().log().string() + ", i * pi/2 = " +
          Complex[F](zero, pi / two).string() + ", delta = " +
          (Complex[F].j().log() - Complex[F](zero, pi / two)).string())
    h.assert_true(Complex[F].j().log() == Complex[F](zero, pi / two))


class iso _TestComplexIdentities[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Identities on Complex[F].
  """
  fun name(): String =>
    "Complex[F]/Identities"

  fun apply(h: TestHelper) =>
    let a = Complex[F](F.from[F64](-2.0), F.from[F64](1.0))
    let b = Complex[F](F.from[F64](1.0), F.from[F64](-2.0))
    let c = Complex[F](F.from[F64](-4.0) / F.from[F64](5.0), F.from[F64](-3.0) / F.from[F64](5.0))
    let d = Complex[F](F.from[F64](5.0), F.from[F64](-3.0))
    
    let one = Complex[F](F.from[F64](1.0), F.from[F64](0.0))
    let two = Complex[F](F.from[F64](2.0), F.from[F64](0.0))

    h.log("a = " + a.string() + ", b = " + b.string() + ", c = " + c.string() +
          ", d = " + d.string())

    // Real part is linear
    h.assert_true((a + b).real() == (a.real() + b.real()))
    h.assert_true((a +~ b).real() ==~ (a.real() +~ b.real()))

    // Imaginary part is linear
    h.assert_true((a + b).imag() == (a.imag() + b.imag()))
    h.assert_true((a +~ b).imag() ==~ (a.imag() +~ b.imag()))

    // Modulus of difference of conjugate
    h.assert_true((a.conj() - b.conj()).abs() == (a - b).abs())
    h.assert_true((a.conj() -~ b.conj()).abs() ==~ (a -~ b).abs())

    // Real part of products
    h.assert_true((a * b).real() == ((a.real() * b.real()) - (a.imag() * b.imag())))
    h.assert_true((a *~ b).real() ==~ ((a.real() *~ b.real()) -~ (a.imag() *~ b.imag())))

    // Imaginary part of products
    h.assert_true((a * b).imag() == ((a.real() * b.imag()) + (a.imag() * b.real())))
    h.assert_true((a *~ b).imag() ==~ ((a.real() *~ b.imag()) +~ (a.imag() *~ b.real())))

    // Conjugate of product
    h.assert_true((a * b).conj() == (a.conj() * b.conj()))
    h.assert_true((a *~ b).conj_unsafe() ==~ (a.conj_unsafe() *~ b.conj_unsafe()))

    // Equality fails with F64
    h.assert_true((a * a.conj()).almost_eq(Complex[F](a.abs() * a.abs())))
    h.assert_true((a *~ a.conj()).almost_eq(Complex[F](a.abs() *~ a.abs())))

    // Conjugate of quotient
    h.assert_true(b.is_null() or ((a / b).conj() == (a.conj() / b.conj())))
    h.assert_true(b.is_null() or ((a /~ b).conj_unsafe() ==~ (a.conj_unsafe() /~ b.conj_unsafe())))

    // Conjugate properties
    h.assert_true((a + a.conj()).is_real())
    h.assert_true(a.real() == ((a + a.conj()) / two).real())
    h.assert_true((a - a.conj()).is_imag())
    h.log("a.imag() = " + a.imag().string() + ", ((a - a.conj()) / 2i).imag() = " +
          ((a - a.conj()) / two).imag().string())
    h.assert_true(a.imag() == ((a - a.conj()) / two).imag())
    
    h.assert_true((a +~ a.conj_unsafe()).is_real())
    h.assert_true(a.real() ==~ ((a +~ a.conj_unsafe()) / two).real())
    h.assert_true((a -~ a.conj_unsafe()).is_imag())
    h.assert_true(a.imag() ==~ ((a -~ a.conj_unsafe()) / two).imag())

    // Equality - Ptolemy's theorem
    let a' = a / Complex[F](a.abs()) // Use default argument for imaginary part
    let b' = b / Complex[F](b.abs())
    let c' = c / Complex[F](c.abs())
    let d' = d / Complex[F](d.abs())
    // a', b', c' and d' are on the same circle. Generalisation of Pythagorean theorem
    h.log("a' = " + a'.string() + ", b' = " + b'.string() + ", c' = " + c'.string() +
          ", d' = " + d'.string())
    h.log("(((a' - b') * (c' - d')) + ((a' - d') * (b' - c'))) = " +
          (((a' - b') * (c' - d')) + ((a' - d') * (b' - c'))).string() +
          ", ((a' - c') * (b' - d')) = " + ((a' - c') * (b' - d')).string() +
          ", delta = " + ((((a' - b') * (c' - d')) + ((a' - d') * (b' - c'))) - ((a' - c') * (b' - d'))).string())
    h.assert_true((((a' - b') * (c' - d')) + ((a' - d') * (b' - c'))).almost_eq(((a' - c') * (b' - d'))))
    h.assert_true((((a' -~ b') *~ (c' -~ d')) +~ ((a' -~ d') *~ (b' -~ c'))).almost_eq(((a' -~ c') *~ (b' -~ d'))))
    
    // Inequality on modulus
    h.assert_true((((a - b).abs() * (c - d).abs()) + ((a - d).abs() * (b - c).abs())) >= ((a - c).abs() * (b - d).abs()))
    h.assert_true((((a -~ b).abs_unsafe() *~ (c -~ d).abs_unsafe()) +~ ((a -~ d).abs_unsafe() *~ (b -~ c).abs_unsafe())) >=~ ((a -~ c).abs_unsafe() *~ (b -~ d).abs_unsafe()))

    // Product of complex
    h.log("(a * a.conj()).real() = " + (a * a.conj()).real().string() +
          ", a.abs() * a.abs() = " + a.abs().powi(2).string())
    // Equality fails with F64
    h.assert_true(Complex[F]((a * a.conj()).real()).almost_eq(Complex[F](a.abs().powi(2))))
    h.assert_true((a * a.conj()).almost_eq(Complex[F](a.abs().powi(2))))
    h.log("(a * a.conj()).imag() = " + (a + a.conj()).imag().string() + " = 0")
    h.assert_true((a * a.conj()).imag() == F.from[USize](0))
    h.log("(a * b) - (b * a) = " + ((a * b) - (b * a)).string())
    h.assert_true((a * b) == (b * a))
    h.log("a * (b + c) = " + (a * (b + c)).string() + ", ((a * b) + (a * c)) = " +
          ((a * b) + (a * c)).string() + ", delta = " +
          ((a * (b + c)) - ((a * b) + (a * c))).string())
    h.assert_true((a * (b + c)).almost_eq((a * b) + (a * c)))
    
    // Scalar product
    // Equality fails with F64
    h.assert_true(Complex[F](a.dotp(a)).almost_eq(Complex[F](a.abs().powi(2))))

    // Power
    h.assert_true((a.powi(4) + one) == ((a.powi(2) + Complex[F].j()) * (a.powi(2) - Complex[F].j())))
    h.assert_true((a.powi(6) + one) == ((a.powi(3) + Complex[F].j()) * (a.powi(3) - Complex[F].j())))
    h.assert_true((a.powi_unsafe(4) +~ one) ==~ ((a.powi_unsafe(2) +~ Complex[F].j()) *~ (a.powi(2) -~ Complex[F].j())))
    h.assert_true((a.powi_unsafe(6) +~ one) ==~ ((a.powi_unsafe(3) +~ Complex[F].j()) *~ (a.powi(3) -~ Complex[F].j())))

    // Square root
    h.log("(2i).sqrt() = " + (two * Complex[F].j()).sqrt().string() +
          " = 1+i, delta = " + ((two * Complex[F].j()).sqrt() - (one + Complex[F].j())).string())
    h.assert_true((two * Complex[F].j()).sqrt().almost_eq(one + Complex[F].j()))
    h.log("(-2i).sqrt() = " + (-two * Complex[F].j()).sqrt().string() +
          " = 1-i, delta = " + ((-two * Complex[F].j()).sqrt() - (one - Complex[F].j())).string())
    h.assert_true((-two * Complex[F].j()).sqrt().almost_eq(one - Complex[F].j()))


interface _FromBits[U: Any val]
  new val from_bits(b: U)
     "Create instance from a bits pattern"

  fun val bits(): U val
     "Give the bits pattern behind the instance"


class iso _TestComplexRandom[F: (_FromBits[U] val & Float & FloatingPoint[F]),
                             U: UnsignedInteger[U] val] is UnitTest
  """
  Random tests on Complex[F].
  """
  fun name(): String =>
    "Complex[F]/Random tests"

  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    //let two = F.from[ISize](2)
    //let three = F.from[ISize](3)
    //let four = F.from[ISize](4)
    //let five = F.from[ISize](5)
    //let pi = F.from[F64](F64.pi())
    
    let rand = Rand
    for i in Range(0, 10) do
      let x = F.from_bits(U.from[U64](rand.next()))
      let y = F.from_bits(U.from[U64](rand.next()))

      h.log("x = " + x.string() + ", y = " + y.string())

      let s1 = Complex[F](x, -y) * Complex[F](x, y)
      let s1' = Complex[F]((x * x) + (y * y), zero)
      // Don't do the test for NaN
      if (not s1.nan()) and (not s1'.nan()) then
        h.log("(x - iy) * (z + iy) = " + s1.string() + ", x * x + y * y = " +
              s1'.string() + ", delta = " + (s1 - s1').string())
        // Equality fails with F64
        h.assert_true(s1.almost_eq(s1'))
      end

      let s2 = Complex[F](x, -one) * Complex[F](x, one)
      let s2' = Complex[F]((x * x) + one, zero)
      if (not s2.nan()) and (not s2'.nan()) then
        h.log("(x - i) * ( x + i) = " + s2.string() + ", x * x + 1 = " +
              s2'.string() + ", delta = " + (s2 - s2').string())
        h.assert_true(s2.almost_eq(s2'))
      end

      let s3 = Complex[F](x, -one) * Complex[F](x, one) * Complex[F](x + one, zero)
      let s3' = Complex[F]((x * x * x) + (x * x) + x + one, zero)
      if (not s3.nan()) and (not s3'.nan()) then
        h.log("(x - i) * (x + i) (x + 1) = " + s3.string() +
              ", x^3 + x^2 + x + 1 = " + s3'.string() + ", delta = " +
              (s3 - s3').string())
        h.assert_true(s3.almost_eq(s3'))
      end

      let s4 = Complex[F](x, -one) * Complex[F](x, one) * Complex[F](x - one, zero)
      let s4' = Complex[F]((x * x * x) - ((x * x) - (x - one)), zero)
      if (not s4.nan()) and (not s4'.nan()) then
        h.log("(x - i) * (x + i) * (x - 1) = " + s4.string() +
              ", x^3 - x^2 + x - 1 = " + s4'.string() + ", delta = " +
              (s4 - s4').string())
        h.assert_true(s4.almost_eq(s4'))
      end

      let s5 = Complex[F](x, -one) * Complex[F](x, one) * Complex[F](x - one, zero) * Complex[F](x + one, zero)
      let s5' = Complex[F]((x * x * x * x) - one, zero)
      if (not s5.nan()) and (not s5'.nan()) then
        h.log("(x - i) * (x + i) * (x - 1) * (x + 1) = " + s5.string() +
            ", x^4 - 1 = " + s5'.string() + ", delta = " + (s5 - s5').string())
        h.assert_true(s5.almost_eq(s5'))
      end

      let a = Complex[F](F.from_bits(U.from[U64](rand.next())), F.from_bits(U.from[U64](rand.next())))
      let b = Complex[F](F.from_bits(U.from[U64](rand.next())), F.from_bits(U.from[U64](rand.next())))

      h.log("a = " + a.string() + ", b = " + b.string())

      // Triangle inequality
      let t1 = (a + b).abs()
      let t1' = a.abs() + b.abs()
      if (not t1.nan()) and (not t1'.nan()) then
        h.log("(a + b).abs() = " + t1.string() + ", (a.abs() + b.abs()) = " + t1'.string())
        h.assert_true(t1 <= t1')

        let t2 = (a +~ b).abs_unsafe()
        let t2' = a.abs_unsafe() +~ b.abs_unsafe()
        h.assert_true(t2 <=~ t2')
      end

      // Difference of modulus
      let t3 = a.abs() - b.abs()
      let t3' = (a - b).abs()
      let t3'' = (a.abs() - b.abs()).abs()
      if (not t3.nan()) and (not t3'.nan()) then
        h.log("a.abs() - b.abs() = " + t3.string() + ", (a - b).abs() = " +
              t3'.string() + ", (a.abs() - b.abs()).abs() = " + t3''.string())
        h.assert_true(t3 <= t3')
        h.assert_true(t3'' <= t3')

        let t4 = a.abs_unsafe() -~ b.abs_unsafe()
        let t4' = (a -~ b).abs_unsafe()
        let t4'' = (a.abs_unsafe() -~ b.abs_unsafe()).abs()
        h.assert_true(t4 <=~ t4')
        h.assert_true(t4'' <=~ t4')
      end
    end


class iso _TestComplexArg[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Test arg() (complex argument / angle) across all quadrants.
  """
  fun name(): String =>
    "Complex[F]/Arg"

  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let pi = F.from[F64](F64.pi())
    let pi_2 = pi / F.from[ISize](2)
    let pi_4 = pi / F.from[ISize](4)
    let three_pi_4 = (F.from[ISize](3) * pi) / F.from[ISize](4)

    // arg(1) = 0  (exact, positive real axis)
    h.assert_true(Complex[F](one, zero).arg() == zero)

    // arg(i) ≈ π/2  (positive imaginary axis)
    h.log("arg(i) = " + Complex[F](zero, one).arg().string() + ", π/2 = " + pi_2.string())
    h.assert_true(Complex[F](Complex[F](zero, one).arg()).almost_eq(Complex[F](pi_2)))

    // arg(-1) ≈ π  (negative real axis — this was the atan() bug)
    h.log("arg(-1) = " + Complex[F](-one, zero).arg().string() + ", π = " + pi.string())
    h.assert_true(Complex[F](Complex[F](-one, zero).arg()).almost_eq(Complex[F](pi)))

    // arg(-i) ≈ -π/2  (negative imaginary axis)
    h.log("arg(-i) = " + Complex[F](zero, -one).arg().string() + ", -π/2 = " + (-pi_2).string())
    h.assert_true(Complex[F](Complex[F](zero, -one).arg()).almost_eq(Complex[F](-pi_2)))

    // arg(1+i) ≈ π/4  (first quadrant)
    h.log("arg(1+i) = " + Complex[F](one, one).arg().string() + ", π/4 = " + pi_4.string())
    h.assert_true(Complex[F](Complex[F](one, one).arg()).almost_eq(Complex[F](pi_4)))

    // arg(-1+i) ≈ 3π/4  (second quadrant — also wrong with old atan())
    h.log("arg(-1+i) = " + Complex[F](-one, one).arg().string() + ", 3π/4 = " + three_pi_4.string())
    h.assert_true(Complex[F](Complex[F](-one, one).arg()).almost_eq(Complex[F](three_pi_4)))

    // arg(-1-i) ≈ -3π/4  (third quadrant)
    h.log("arg(-1-i) = " + Complex[F](-one, -one).arg().string() + ", -3π/4 = " + (-three_pi_4).string())
    h.assert_true(Complex[F](Complex[F](-one, -one).arg()).almost_eq(Complex[F](-three_pi_4)))

    // arg(1-i) ≈ -π/4  (fourth quadrant)
    h.log("arg(1-i) = " + Complex[F](one, -one).arg().string() + ", -π/4 = " + (-pi_4).string())
    h.assert_true(Complex[F](Complex[F](one, -one).arg()).almost_eq(Complex[F](-pi_4)))

    // arg and abs recover z: z = abs * exp(i * arg)
    let z = Complex[F](F.from[F64](-1.5), F.from[F64](2.0))
    let rho = z.abs()
    let theta = z.arg()
    let recovered = Complex[F](rho * theta.cos(), rho * theta.sin())
    h.log("z = " + z.string() + ", recovered = " + recovered.string())
    h.assert_true(z.almost_eq(recovered))


class iso _TestComplexExp[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Test exp() and exp_unsafe().
  """
  fun name(): String =>
    "Complex[F]/Exp"

  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let pi = F.from[F64](F64.pi())
    let pi_2 = pi / F.from[ISize](2)

    // exp(0) = 1
    h.log("exp(0) = " + Complex[F](zero, zero).exp().string())
    h.assert_true(Complex[F](zero, zero).exp().almost_eq(Complex[F](one, zero)))

    // exp(1+0i) = e  (matches real exponential)
    h.log("exp(1) = " + Complex[F](one, zero).exp().string() + ", e = " + one.exp().string())
    h.assert_true(Complex[F](one, zero).exp().almost_eq(Complex[F](one.exp(), zero)))

    // Euler's identity: exp(iπ) = -1
    h.log("exp(iπ) = " + Complex[F](zero, pi).exp().string())
    h.assert_true(Complex[F](zero, pi).exp().almost_eq(Complex[F](-one, zero)))

    // exp(iπ/2) = i
    h.log("exp(iπ/2) = " + Complex[F](zero, pi_2).exp().string())
    h.assert_true(Complex[F](zero, pi_2).exp().almost_eq(Complex[F](zero, one)))

    // exp(1 + iπ/2) = e*i
    let e = one.exp()
    h.log("exp(1 + iπ/2) = " + Complex[F](one, pi_2).exp().string() + ", e*i = " + Complex[F](zero, e).string())
    h.assert_true(Complex[F](one, pi_2).exp().almost_eq(Complex[F](zero, e)))

    // exp_unsafe gives the same results for normal inputs
    h.assert_true(Complex[F](zero, zero).exp_unsafe().almost_eq(Complex[F](one, zero)))
    h.assert_true(Complex[F](zero, pi).exp_unsafe().almost_eq(Complex[F](-one, zero)))
    h.assert_true(Complex[F](zero, pi_2).exp_unsafe().almost_eq(Complex[F](zero, one)))


class iso _TestComplexCPow[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Test pow() (complex exponentiation with complex exponent).
  """
  fun name(): String =>
    "Complex[F]/CPow"

  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    let three = F.from[ISize](3)
    let half = one / two

    // 2^3 = 8
    let two_r = Complex[F](two, zero)
    let three_r = Complex[F](three, zero)
    h.log("2.pow(3) = " + two_r.pow(three_r).string())
    h.assert_true(two_r.pow(three_r).almost_eq(Complex[F](F.from[ISize](8), zero)))

    // (-1)^(1/2) = i  (square root of -1 via complex pow)
    let minus_one = Complex[F](-one, zero)
    let half_r = Complex[F](half, zero)
    h.log("(-1).pow(0.5) = " + minus_one.pow(half_r).string())
    h.assert_true(minus_one.pow(half_r).almost_eq(Complex[F](zero, one)))

    // i^2 = -1
    let i = Complex[F](zero, one)
    let two_r2 = Complex[F](two, zero)
    h.log("i.pow(2) = " + i.pow(two_r2).string())
    h.assert_true(i.pow(two_r2).almost_eq(Complex[F](-one, zero)))

    // i^4 = 1
    let four_r = Complex[F](F.from[ISize](4), zero)
    h.log("i.pow(4) = " + i.pow(four_r).string())
    h.assert_true(i.pow(four_r).almost_eq(Complex[F](one, zero)))


class iso _TestComplexEdge[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Edge cases: NaN/Inf in string(), dotp_unsafe(), div by zero, log(0).
  """
  fun name(): String =>
    "Complex[F]/Edge cases"

  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let inf = one / zero
    let nan = zero / zero

    // string() with special float components — just verify no crash, non-empty output
    h.log("Complex(inf, 0).string() = " + Complex[F](inf, zero).string())
    h.log("Complex(nan, 0).string() = " + Complex[F](nan, zero).string())
    h.log("Complex(0, inf).string() = " + Complex[F](zero, inf).string())
    h.assert_true(Complex[F](inf, zero).string().size() > 0)
    h.assert_true(Complex[F](nan, zero).string().size() > 0)
    h.assert_true(Complex[F](zero, inf).string().size() > 0)

    // is_real with explicit tolerance: tiny imaginary part
    // epsilon()/2 is below the default tol (epsilon()), but above zero
    let tiny_im = Complex[F](one, F.epsilon() / F.from[ISize](2))
    h.assert_true(tiny_im.is_real())                           // |im|=eps/2 < default tol: passes
    h.assert_false(tiny_im.is_real(F.from[ISize](0)))          // strict zero tol: fails

    // div by zero: result should be infinite or NaN, not a valid finite number
    let one_r = Complex[F](one, zero)
    let zero_r = Complex[F](zero, zero)
    let d = one_r / zero_r
    h.log("1 / 0 = " + d.string())
    h.assert_false(d.finite())

    // log(0): abs(0).log() = -inf
    let log0 = zero_r.log()
    h.log("log(0) = " + log0.string())
    h.assert_true(log0.infinite())

    // dotp_unsafe: re*re + im*im (2D dot product)
    let a = Complex[F](F.from[ISize](3), F.from[ISize](4))
    let dot_self = a.dotp_unsafe(a)
    h.log("(3+4i).dotp_unsafe(3+4i) = " + dot_self.string() + " (expect 25)")
    h.assert_true(dot_self == F.from[ISize](25))

    let b = Complex[F](one, one)
    let c = Complex[F](one, -one)
    let dot_orth = b.dotp_unsafe(c)
    h.log("(1+i).dotp_unsafe(1-i) = " + dot_orth.string() + " (expect 0)")
    h.assert_true(dot_orth == zero)

