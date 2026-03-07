// Tests of FFT implementations
//
// Some of these tests compare the results of FFT class with those of a naive DFT
// implementations. The comparisons of both values, arrays of complex numbers,
// is done relative to the modulus of the complex numbers, using appromimated
// equality `almost_eq` with the default tolerance values (`F.epsilon().sqrt()`).

use "format"

use "collections"
use "random"
use "time"

use "../mathx"
use "../pony_testx"


primitive _FFT[F: (Float & FloatingPoint[F])]
  """
  A reference DFT that does not use FFT algorithm. Also includes a naive
  convolution.
  """
  
  fun naive_fft(a: Array[Complex[F]], inverse: Bool = false): Array[Complex[F]]^ =>
    """
    Naive DFT implementation. Don't use it for real DFT calculations. `O(a.size^2)` performance.
    """
    let size = a.size()
    let coef = if inverse then
        F.from[USize](2) * F.from[F64](F64.pi())
      else
        -F.from[USize](2) * F.from[F64](F64.pi())
      end
    let result = Array[Complex[F]](size)
    try
      for k in Range[USize](0, size) do
        var sum: Complex[F] = Complex[F]
        for t in Range[USize](0, size) do
          // More accurate to take modulo before division
          let theta = (coef * F.from[USize]((t * k) % size)) / F.from[USize](size)
          sum = sum + Complex[F]((a(t)?.real() * theta.cos()) - (a(t)?.imag() * theta.sin()),
                                 (a(t)?.real() * theta.sin()) + (a(t)?.imag() * theta.cos()))
        end
        result.push(sum)
      end
    end

    // If inverse, scale the result
    if inverse then
      let scale = Complex[F](F.from[USize](size))
      try
        for i in Range(0, size) do
          result.update(i, result(i)? / scale)?
        end
      end
    end
    result


  fun naive_convolve(a: Array[Complex[F]], b: Array[Complex[F]]): Array[Complex[F]]^ =>
    """
    Naive convolution implementation.
    
    Don't use it. `O(a.size * b.size)` performance.
    `a` and `b` must have the same size.
    """
    let size = a.size()
    let result = Array[Complex[F]](size)
    for i in Range[USize](0, size) do
      result.push(Complex[F])
    end
    try
      for i in Range[USize](0, size) do
        for j in Range[USize](0, size) do
          let k = (i + j) % size
          result.update(k, result(k)? + Complex[F]((a(i)?.real() * b(j)?.real()) - (a(i)?.imag() * b(j)?.imag()),
                                                   (a(i)?.real() * b(j)?.imag()) + (a(i)?.imag() * b(j)?.real())))?
        end
      end
    end
    result


class iso _TestFFT[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Tests on FFT with Complex[F]
  """

  fun name(): String =>
    "fft/fourier"

  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    // FFT of impulse is a constant
    let impulse: Array[Complex[F]] = []
    impulse.push(Complex[F](one))
    for i in Range(1, 256) do
      impulse.push(Complex[F](zero))
    end

    let impulse_fft: Array[Complex[F]] = []
    for i in Range(0, 256) do
      impulse_fft.push(Complex[F](one))
    end

    let imp = FFT[F].fourier(impulse)
    h.assert_array_eq[Complex[F]](impulse_fft, imp, "FFT of impulse is constant")

    // FFT of constant is impulse
    let constant: Array[Complex[F]] = []
    for i in Range(0, 512) do
      constant.push(Complex[F](one))
    end

    let constant_fft: Array[Complex[F]] = []
    constant_fft.push(Complex[F](F.from[ISize](512)))
    for i in Range(1, 512) do
      constant_fft.push(Complex[F](zero))
    end

    let const = FFT[F].fourier(constant)
    h.assert_array_eq[Complex[F]](constant_fft, const, "FFT of constant is impulse")

    // FFT is linear
    let a: Array[Complex[F]] = []
    let b: Array[Complex[F]] = []
    let c: Array[Complex[F]] = []
    let rand = Rand
    let x = Complex[F](F.from[F64](rand.i32().f64() * rand.real()),
                       F.from[F64](rand.i32().f64() * rand.real()))
    let y = Complex[F](F.from[F64](rand.i32().f64() * rand.real()),
                       F.from[F64](rand.i32().f64() * rand.real()))

    for i in Range(0, 1024) do
      let ca = Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real()))
      let cb = Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real()))
      let cc = (x * ca) + (y * cb)
      a.push(ca)
      b.push(cb)
      c.push(cc)
    end

    let a_fft = FFT[F].fourier(a)
    let b_fft = FFT[F].fourier(b)
    let c_fft = FFT[F].fourier(c)

    try
      for i in Range(0, 1024) do
        // To avoid impact of calculation errors, we approximate equality with relative tolerance
        h.assert_true(((x * a_fft(i)?) + (y * b_fft(i)?)).almost_eq(c_fft(i)?), "Linear operator")
      end
    else
      h.fail("Error while testing FFT linearity")
    end

    // Inverse FFT of FFT gives back original
    let reverse: Array[Complex[F]] = []
    let reverse_fft: Array[Complex[F]] = []
    for i in Range(0, 2048) do
      let cx = Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real()))
      reverse.push(cx)
      reverse_fft.push(cx)
    end

    let original = FFT[F].fourier(FFT[F].fourier(reverse), true)
    h.assert_array_almost_eq[Complex[F], F](reverse_fft, original, "Reverse FFT")

    // A more complex function
    let sinus: Array[Complex[F]] = []
    let sinus_fft: Array[Complex[F]] = []
    let pi2 = F.from[F64](F64.pi() * 2.0)
    let two = F.from[ISize](2)
    let five = F.from[ISize](5)
    let hundred = F.from[ISize](100)
    let range = F.from[ISize](256)
    for i in Range(0, 256) do
      let cx = Complex[F]((five * ((pi2 * F.from[USize](i)) / range).cos()) +
                          (two * (((pi2 * F.from[USize](i)) + hundred) / range).cos()),
                          ((pi2 * F.from[USize](i)) / range).sin())
      sinus.push(cx)

      if i == 1 then
        sinus_fft.push(Complex[F](F.from[F64](1004.71584293556111333601), F.from[F64](97.47620070205202580382)))
      elseif i == 255 then
        sinus_fft.push(Complex[F](F.from[F64](748.71584293556111333601), F.from[F64](-97.47620070205141473707)))
      else
        sinus_fft.push(Complex[F])
      end
    end

    let sfft = FFT[F].fourier(sinus)
    h.assert_array_almost_eq[Complex[F], F](sinus_fft, sfft, "Combined cosinus + sinus")

    // Compare with naive implementation
    let four = Array[Complex[F]](2048)
    let naiv = Array[Complex[F]](2048)
    for i in Range(0, 2048) do
      let cx = Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real()))
      four.push(cx)
      naiv.push(cx)
    end
    let orig = four.clone()
    let four_fft = FFT[F].fourier(four)
    let naiv_fft = _FFT[F].naive_fft(naiv)
    h.assert_array_almost_eq[Complex[F], F](four_fft, naiv_fft, "Naive implementation")

    // Now inverse
    let rev = FFT[F].fourier(four_fft, true)
    h.assert_array_almost_eq[Complex[F], F](rev, orig, "Inverse fourier")

    // Parseval's theorem: sum(|x[n]|^2) == (1/N) * sum(|X[k]|^2)
    // Energy is preserved under DFT (up to a 1/N normalization factor)
    let px: Array[Complex[F]] = []
    let parseval_n: USize = 256
    for _ in Range(0, parseval_n) do
      px.push(Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real())))
    end
    var energy_x: F = zero
    try
      for i in Range(0, parseval_n) do
        energy_x = energy_x + px(i)?.abs2()
      end
    end
    let px_fft = FFT[F].fourier(px)
    var energy_fft: F = zero
    try
      for i in Range(0, parseval_n) do
        energy_fft = energy_fft + px_fft(i)?.abs2()
      end
    end
    // Parseval: energy_x == energy_fft / N
    let parseval_lhs = energy_x
    let parseval_rhs = energy_fft / F.from[USize](parseval_n)
    h.assert_true(
      Complex[F](parseval_lhs).almost_eq(Complex[F](parseval_rhs)),
      "Parseval's theorem")


class iso _TestFFTUnsafe[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Tests on `fourier_unsafe`, which must produce the same results as `fourier`.
  """

  fun name(): String =>
    "fft/fourier_unsafe"

  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    // FFT of impulse is a constant
    let impulse: Array[Complex[F]] = []
    impulse.push(Complex[F](one))
    for _ in Range(1, 256) do
      impulse.push(Complex[F](zero))
    end
    let impulse_fft: Array[Complex[F]] = []
    for _ in Range(0, 256) do
      impulse_fft.push(Complex[F](one))
    end
    h.assert_array_eq[Complex[F]](impulse_fft,
      FFT[F].fourier_unsafe(impulse), "fourier_unsafe of impulse is constant")

    // FFT of constant is impulse
    let constant: Array[Complex[F]] = []
    for _ in Range(0, 512) do
      constant.push(Complex[F](one))
    end
    let constant_fft: Array[Complex[F]] = []
    constant_fft.push(Complex[F](F.from[ISize](512)))
    for _ in Range(1, 512) do
      constant_fft.push(Complex[F](zero))
    end
    h.assert_array_eq[Complex[F]](constant_fft,
      FFT[F].fourier_unsafe(constant), "fourier_unsafe of constant is impulse")

    // fourier_unsafe agrees with fourier on random input, and inverse round-trips
    let rand = Rand
    for n in Range(4, 12) do
      let pow2 = UnsignedComp.pow2(n)
      let a: Array[Complex[F]] = []
      let b: Array[Complex[F]] = []
      for _ in Range(0, pow2) do
        let cx = Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real()))
        a.push(cx)
        b.push(cx)
      end
      let orig = a.clone()
      let a_fft = FFT[F].fourier(a)
      let b_fft = FFT[F].fourier_unsafe(b)
      h.assert_array_almost_eq[Complex[F], F](a_fft, b_fft,
        "fourier_unsafe == fourier, n=" + n.string())

      // Inverse round-trip
      let rev = FFT[F].fourier_unsafe(b_fft, true)
      h.assert_array_almost_eq[Complex[F], F](orig, rev,
        "fourier_unsafe inverse, n=" + n.string())
    end


class iso _TestFFTComplex[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Tests on FFT with linearized complex.
  We only need to test the same results as TestFFT.
  """

  fun name(): String =>
    "fft/fourier_complex"
    
  fun apply(h: TestHelper) =>
    let rand = Rand
    for n in Range(4, 12) do
      try
        let pow2 = UnsignedComp.pow2(n)
        let a: Array[Complex[F]] = []
        let b: Array[F] = []

        for _ in Range(0, pow2) do
          let real = F.from[F64](rand.real() * I32.max_value().f64())
          let imag = F.from[F64](rand.real() * I32.max_value().f64())
          a.push(Complex[F](real, imag))
          b.push(real)
          b.push(imag)
        end
      	let orig = b.clone()

        let a_fft = _FFT[F].naive_fft(a)
        let b_fft = FFT[F].fourier_complex(b)

        for i in Range(0, pow2) do
          h.assert_true(a_fft(i)?.almost_eq(Complex[F](b_fft(2 * i)?, b_fft((2 * i) + 1)?)), "fourier_complex, i=" + i.string())
          if not a_fft(i)?.almost_eq(Complex[F](b_fft(2 * i)?, b_fft((2 * i) + 1)?)) then
            h.log("i=" + i.string() + ", a_fft=" + a_fft(i)?.string() + " / " + b_fft(2 * i)?.string() + ", " + b_fft((2 * i) + 1)?.string())
          end
        end

        // Now inverse
        let rev = FFT[F].fourier_complex(b_fft, true)

        for i in Range(0, pow2) do
          h.assert_true(((rev(i)? - orig(i)?) / orig(i)?).abs() < F.from[F64](1e-3), "Inverse fourier_complex, i=" + i.string())
          if not (((rev(i)? - orig(i)?) / orig(i)?).abs() < F.from[F64](1e-3)) then
            h.log("i=" + i.string() + ", rev=" + rev(i)?.string() + " / " + orig(i)?.string())
          end
	      end
      else
        h.fail("Equality naive fourier == fourier_complex")
      end
    end


class iso _TestFFTReal[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Tests on FFT for real series
  """

  fun name(): String =>
    "fft/fourier_real"
    
  fun apply(h: TestHelper) =>
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    // FFT of impulse is a constant
    let impulse = Array[F](256)
    impulse.push(one)
    for i in Range(1, 256) do
      impulse.push(zero)
    end

    let original_imp = impulse.clone()

    let impulse_fft = Array[F](256)
    impulse_fft.push(one) // Real part of FFT[0]
    impulse_fft.push(one) // Real part of FFT[127]
    for i in Range(1, 128) do
      impulse_fft.push(one)  // Real part
      impulse_fft.push(zero) // Imaginary part
    end

    let imp = FFT[F].fourier_real(impulse)
    h.assert_array_eq[F](impulse_fft, imp, "FFT of impulse is constant")
    let orig_imp = FFT[F].fourier_real(imp, true)
    h.assert_array_eq[F](original_imp, orig_imp, "Reverse FFT of constant is impulse")

    // FFT of constant is impulse
    let constant = Array[F](512)
    for i in Range(0, 512) do
      constant.push(one)
    end

    let original_const = constant.clone()

    let constant_fft = Array[F](512)
    constant_fft.push(F.from[ISize](512))
    for i in Range(1, 512) do
      constant_fft.push(zero)
    end

    let const = FFT[F].fourier_real(constant)
    h.assert_array_eq[F](constant_fft, const, "FFT of constant is impulse")
    let orig_const = FFT[F].fourier_real(const, true)
    h.assert_array_eq[F](original_const, orig_const, "Reverse FFT of impulse is constant")

    let rand = Rand
    for n in Range(4, 12) do
      try
        let pow2 = UnsignedComp.pow2(n)
        let a = Array[Complex[F]](pow2)
        let b = Array[F](pow2)
      	let orig = Array[F](pow2)

        for _ in Range(0, pow2) do
          let real = F.from[F64](rand.real() * I32.max_value().f64())
          a.push(Complex[F](real, F.from[F64](0.0)))
          b.push(real)
	        orig.push(real)
        end

        let a_fft = _FFT[F].naive_fft(a)
        let b_fft = FFT[F].fourier_real(b)

        h.assert_almost_eq[Complex[F], F](Complex[F](b_fft(0)?), Complex[F](a_fft(0)?.real()), "fourier_real(0)")
        h.assert_almost_eq[Complex[F], F](Complex[F](b_fft(1)?), Complex[F](a_fft(pow2 / 2)?.real()), "fourier_real(pow2 / 2)")
        for i in Range(1, pow2 / 2) do
      	  let c_i = Complex[F](b_fft(2 * i)?, b_fft((2 * i) + 1)?)
          h.assert_almost_eq[Complex[F], F](c_i, a_fft(i)?, "fourier_real, i=" + i.string())
          if not c_i.almost_eq(a_fft(i)?) then
            h.log("i=" + i.string() + ", a_fft=" + a_fft(i)?.string() + " / " + b_fft(2 * i)?.string() + ", " + b_fft((2 * i) + 1)?.string())
          end
        end

        // Inverse FFT
        let rev = FFT[F].fourier_real(b_fft, true)
        // We can't use assert_array_almost_eq here because F32 and F64 don't have the
        // trait Approximated.
        //h.assert_array_almost_eq[F, F](rev, orig, "Inverse fourier_real")
        for i in Range(0, pow2) do
          h.assert_true(((rev(i)? - orig(i)?) / orig(i)?).abs() < F.from[F64](1e-2), "Inverse fourier_real, i=" + i.string())
          if not (((rev(i)? - orig(i)?) / orig(i)?).abs() < F.from[F64](1e-2)) then
            h.log("i=" + i.string() + ", rev=" + rev(i)?.string() + " / " + orig(i)?.string() + " eps=" + ((rev(i)? - orig(i)?) / orig(i)?).abs().string())
          end
        end

      else
        h.fail("Equality naive fourier == fourier_real")
      end
    end


class iso _TestFFT2[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Tests on FFT using trigonometric table.
  We only need to test the same results as naive FFT..
  """

  fun name(): String =>
    "fft/fourier2"
    
  fun apply(h: TestHelper) =>
    let rand = Rand
    for n in Range(4, 12) do
      let pow2 = UnsignedComp.pow2(n)
      let a: Array[Complex[F]] = []
      let b: Array[Complex[F]] = []

      for _ in Range(0, pow2) do
        let real = F.from[F64](rand.real() * I32.max_value().f64())
        let imag = F.from[F64](rand.real() * I32.max_value().f64())
        a.push(Complex[F](real, imag))
        b.push(Complex[F](real, imag))
      end

      let orig = a.clone()
      let a_fft = _FFT[F].naive_fft(a)
      let b_fft = FFT[F].fourier2(b)

      h.assert_array_almost_eq[Complex[F], F](a_fft, b_fft, "fourier2 almost equal")

      // Now inverse
      let rev = FFT[F].fourier2(b_fft, true)
      h.assert_array_almost_eq[Complex[F], F](rev, orig, "Inverse fourier2")
    end


class iso _TestFFTConvolution[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Tests on circular convolution
  """

  fun name(): String =>
    "fft/convolution"
    
  fun apply(h: TestHelper) =>
    let rand = Rand
    for n in Range(4, 12) do
      let pow2 = UnsignedComp.pow2(n)
      let a: Array[Complex[F]] = []
      let b: Array[Complex[F]] = []

      for _ in Range(0, pow2) do
        let real_a = F.from[F64](rand.real() * I32.max_value().f64())
        let imag_a = F.from[F64](rand.real() * I32.max_value().f64())
        a.push(Complex[F](real_a, imag_a))
        let real_b = F.from[F64](rand.real() * I32.max_value().f64())
        let imag_b = F.from[F64](rand.real() * I32.max_value().f64())
        b.push(Complex[F](real_b, imag_b))
      end

      let naive_conv = _FFT[F].naive_convolve(a, b)
      let convolution = FFT[F].convolve(a, b)
      h.assert_array_almost_eq[Complex[F], F](convolution, naive_conv, "Almost equal convolutions")
    end

    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    // Convolution identity: convolve(a, impulse) == a
    // impulse[0]=1, rest=0 is the identity element for circular convolution
    let n_id: USize = 256
    let sig = Array[Complex[F]](n_id)
    let imp_id = Array[Complex[F]](n_id)
    imp_id.push(Complex[F](one))
    for _ in Range(1, n_id) do
      imp_id.push(Complex[F](zero))
    end
    for _ in Range(0, n_id) do
      sig.push(Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real())))
    end
    let sig_orig = sig.clone()
    // convolve clones its inputs internally so sig is not modified
    h.assert_array_almost_eq[Complex[F], F](sig_orig, FFT[F].convolve(sig, imp_id),
      "Convolution with impulse is identity")

    // Convolution commutativity: convolve(a, b) == convolve(b, a)
    let n_comm: USize = 128
    let ca = Array[Complex[F]](n_comm)
    let cb = Array[Complex[F]](n_comm)
    for _ in Range(0, n_comm) do
      ca.push(Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real())))
      cb.push(Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real())))
    end
    let ab = FFT[F].convolve(ca.clone(), cb.clone())
    let ba = FFT[F].convolve(cb, ca)
    h.assert_array_almost_eq[Complex[F], F](ab, ba, "Convolution is commutative")


class iso _TestFFTBluestein[F: (Float & FloatingPoint[F])] is UnitTest
  """
  Tests on Bluestein DFT.
  """

  fun name(): String =>
    "fft/bluestein"
    
  fun apply(h: TestHelper) =>
    let rand = Rand
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    for n in Range(64, 65) do
      // FFT of impulse is a constant
      let impulse = Array[Complex[F]](n)
      impulse.push(Complex[F](one))
      for i in Range(1, n) do
        impulse.push(Complex[F](zero))
      end

      let impulse_fft = Array[Complex[F]](n)
      for i in Range(0, n) do
        impulse_fft.push(Complex[F](one))
      end

      let imp = FFT[F].bluestein(impulse)
      h.assert_array_almost_eq[Complex[F], F](impulse_fft, imp, "Bluestein FFT of impulse is constant")

      // FFT of constant is impulse
      let constant = Array[Complex[F]](n)
      for i in Range(0, n) do
        constant.push(Complex[F](one))
      end

      let constant_fft = Array[Complex[F]](n)
      constant_fft.push(Complex[F](F.from[USize](n)))
      for i in Range(1, n) do
        constant_fft.push(Complex[F](zero))
      end

      let const = FFT[F].bluestein(constant)
      h.assert_array_almost_eq[Complex[F], F](constant_fft, const, "Bluestein FFT of constant is impulse")

      // FFT is linear
      let a = Array[Complex[F]](n)
      let b = Array[Complex[F]](n)
      let c = Array[Complex[F]](n)
      let a_orig = Array[Complex[F]](n)
      let x = Complex[F](F.from[F64](rand.i32().f64() * rand.real()),
                         F.from[F64](rand.i32().f64() * rand.real()))
      let y = Complex[F](F.from[F64](rand.i32().f64() * rand.real()),
                         F.from[F64](rand.i32().f64() * rand.real()))

      for i in Range(0, n) do
        let ca = Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real()))
        let cb = Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real()))
        let cc = (x * ca) + (y * cb)
        a.push(ca)
        b.push(cb)
        c.push(cc)
      	a_orig.push(ca)
      end

      let a_fft = FFT[F].bluestein(a)
      let b_fft = FFT[F].bluestein(b)
      let c_fft = FFT[F].bluestein(c)

      try
        for i in Range(0, n) do
          // To avoid impact of calculation errors, we approximate equality with relative tolerance
          h.assert_almost_eq[Complex[F], F]((x * a_fft(i)?) + (y * b_fft(i)?), c_fft(i)?, "Linear operator")
        end
      else
        h.fail("Error while testing Bluestein FFT linearity")
      end

      // Check that we have a DFT when n is a power of 2
      if (n and (n - 1)) == 0 then
        let a_orig_fft = FFT[F].fourier2(a_orig)
	      h.assert_array_almost_eq[Complex[F], F](a_fft, a_orig_fft, "Bluestein DFT")
      end
    end

    // Non-power-of-2 sizes: verify Bluestein against naive DFT and test inverse round-trip
    let non_pow2_sizes: Array[USize] = [3; 5; 10; 100; 127; 500]
    for n2 in non_pow2_sizes.values() do
      let input = Array[Complex[F]](n2)
      let input_orig = Array[Complex[F]](n2)
      for _ in Range(0, n2) do
        let cx = Complex[F](F.from[F64](rand.real()), F.from[F64](rand.real()))
        input.push(cx)
        input_orig.push(cx)
      end

      // naive_fft reads input without modifying it
      let naive = _FFT[F].naive_fft(input)
      // bluestein modifies input in-place and returns it
      let blue = FFT[F].bluestein(input)
      h.assert_array_almost_eq[Complex[F], F](naive, blue,
        "Bluestein vs naive DFT, n=" + n2.string())

      // Inverse round-trip
      let rev = FFT[F].bluestein(blue, true)
      h.assert_array_almost_eq[Complex[F], F](input_orig, rev,
        "Bluestein inverse, n=" + n2.string())
    end

