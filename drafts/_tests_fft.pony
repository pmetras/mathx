// Tests for Fast Fourier Transform


use "collections"
use "random"

use "../mathx"
use "../pony_testx"


primitive Misc
  """
  Miscelaneous functions
  """

  fun copy_carray[F: FloatingPoint[F] val](a: Array[Complex[F val] val]): Array[Complex[F val] val] =>
    """
    Create a copy of the array of complex numbers.
    """
    let res = Array[Complex[F val] val].create(a.size())
    try
      for i in Range(0, a.size()) do
        let c = Complex[F](a(i)?.real(), a(i)?.imag())
        res.push(c)
      end
    end
    res

  fun copy_darray[F: FloatingPoint[F] val](a: Array[F val] val): Array[F val] iso^ =>
    """
    Create a copy of the array of floating points.
    """
    let res = Array[F].create(a.size())
    try
      for i in Range(0, a.size()) do
        res.push(a(i)?)
      end
    end
    res


class iso _TestFFT[F: FloatingPoint[F] val] is UnitTest
  """
  FFT on Complex[F].
  """
  fun name(): String =>
    "FFT[Complex[F]]"

  fun apply(h: TestHelper) =>
    let a: Array[Complex[F] val] = Array[Complex[F]]
    for i in Range(1, 5) do
      a.push(Complex[F](F.from[USize](i), F.from[USize](0)))
    end
    h.assert_true(a.size() == 4)

    let b = Misc[Complex[F] val].copy_carray(a)

    // Calculate Fourier transform
    let a' = recover a end
    FFT[F].fourier(a')

    let zero = F.from[F64](0.0)
    let half = F.from[F64](0.5)
    let c = [Complex[F](F.from[F64](2.5), zero); Complex[F](-half, -half)
             Complex[F](-half, zero); Complex[F](-half, half)]
    h.assert_array_eq[Complex[F]](a', c)

    // Calculate inverse
    FFT[F].fournier(a', -1)
    // Back to original...
    h.assert_array_eq[Complex[F]](a', b)


