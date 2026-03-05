// Tests of FFT implementations
//
// Some of these tests compare the results of FFT class with those of a naive DFT
// implementations. The comparisons of both values, arrays of complex numbers,
// is done relative to the modulus of the complex numbers. You'll see that two
// complex are considered almost equals with a precision of 1e-2 or 1e-3. These
// high error accuracy numbers results of the tests being done with big random
// numbers being used to detect if overflows occur. In practice, relative errors
// are much lower when using these functions.

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
        h.assert_true(((x * a_fft(i)?) + (y * b_fft(i)?)).almost_eq(c_fft(i)?, F.from[F64](1e-6)), "Linear operator")
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
    // TODO DELETE
    let original = FFT[F].fourier(FFT[F].fourier(reverse), true)
    try
      for i in Range(0, 2048) do
        h.assert_true(reverse_fft(i)?.almost_eq(original(i)?, F.from[F64](1e-6)), "Reverse FFT")
      end
    else
      h.fail("Reverse FFT")
    end
    //h.assert_array_almost_eq[Complex[F], F](reverse_fft, original, "Reverse FFT")
    // TODO HERE
    h.assert_array_almost_eq[Complex[F], F](reverse_fft, original, F.from[F64](1e-6), F.from[F64](0.0), "Reverse FFT")

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

    // TODO DELETE
    try
      for i in Range(0, 256) do
        h.assert_true(sinus_fft(i)?.almost_eq(sfft(i)?, F.from[F64](1e-6)), "Combined cosinus + sinus, i=" + i.string())
      end
    else
      h.fail("Combined sinus + cosinus")
    end
    //h.assert_array_almost_eq[Complex[F], F](sinus_fft, sfft, "Combined cosinus + sinus")
    // TODO HERE
    h.assert_array_almost_eq[Complex[F], F](sinus_fft, sfft, F.from[F64](1e-6), F.from[F64](0.0), "Combined cosinus + sinus")

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
    // TODO DELETE
    try
      for i in Range(0, 2048) do
        // TODO: We need to lower accuracy to 0.05! Higher than 1e-2!
        h.assert_true(four_fft(i)?.almost_eq(naiv_fft(i)?, F.from[F64](1e-6)), "Naive implementation, i=" + i.string())
        if not four_fft(i)?.almost_eq(naiv_fft(i)?, F.from[F64](1e-6)) then
          h.log("i=" + i.string() + ", four_fft=" + four_fft(i)?.string() + " / " + naiv_fft(i)?.string())
        end
      end
    else
      h.fail("Naive implementation")
    end
    //h.assert_array_almost_eq[Complex[F], F](four_fft, naiv_fft, "Naive implementation")
    // TODO HERE
    h.assert_array_almost_eq[Complex[F], F](four_fft, naiv_fft, F.from[F64](1e-6), F.from[F64](0.0), "Naive implementation")

    // Now inverse
    let rev = FFT[F].fourier(four_fft, true)

    // TODO DELETE
    try
      for i in Range(0, 2048) do
        h.assert_true(rev(i)?.almost_eq(orig(i)?, F.from[F64](1e-6)), "Inverse fourier, i=" + i.string())
        if not rev(i)?.almost_eq(orig(i)?, F.from[F64](1e-6)) then
          h.log("i=" + i.string() + ", rev=" + rev(i)?.string() + " / " + orig(i)?.string())
        end
      end
    else
      h.fail("Inverse fourier")
    end
    //h.assert_array_almost_eq[Complex[F], F](rev, orig, "Inverse fourier")
    // TODO HERE
    h.assert_array_almost_eq[Complex[F], F](rev, orig, F.from[F64](1e-6), F.from[F64](0.0), "Inverse fourier")


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
          h.assert_true(a_fft(i)?.almost_eq(Complex[F](b_fft(2 * i)?, b_fft((2 * i) + 1)?), F.from[F64](1e-6)), "fourier_complex, i=" + i.string())
          if not a_fft(i)?.almost_eq(Complex[F](b_fft(2 * i)?, b_fft((2 * i) + 1)?), F.from[F64](1e-6)) then
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

        // TODO DELETE h.assert_true(((a_fft(0)?.real() - b_fft(0)?) / a_fft(0)?.real()).abs() < F.from[F64](1e-3), "fourier_real(0)")
        // TODO DELETE h.assert_true(((a_fft(pow2 / 2)?.real() - b_fft(1)?) / a_fft(pow2 / 2)?.real()).abs() < F.from[F64](1e-3), "fourier_real(pow2 / 2)")
        h.assert_almost_eq[F, F](b_fft(0)?, a_fft(0)?.real()?, F.from[F64](1e-6), F.from[F64](0.0), "fourier_real(0)")
        h.assert_almost_eq[F, F](b_fft(1)?, a_fft(pow2 / 2)?.real(), F.from[F64](1e-6), F.from[F64](0.0), "fourier_real(pow2 / 2)")
        for i in Range(1, pow2 / 2) do
	  let c_i = Complex[F](b_fft(2 * i)?, b_fft((2 * i) + 1)?)
          h.assert_almost_eq[Complex[F], F](c_i, a_fft(i)?, F.from[F64](1e-6), F.from[F64](0.0), "fourier_real, i=" + i.string())
          if not c_i.almost_eq(a_fft(i)?, F.from[F64](1e-6)) then
            h.log("i=" + i.string() + ", a_fft=" + a_fft(i)?.string() + " / " + b_fft(2 * i)?.string() + ", " + b_fft((2 * i) + 1)?.string())
          end
        end

	// Inverse FFT
	let rev = FFT[F].fourier_real(b_fft, true)

        // TODO DELETE
	for i in Range(0, pow2) do
          h.assert_true(((rev(i)? - orig(i)?) / orig(i)?).abs() < F.from[F64](1e-2), "Inverse fourier_real, i=" + i.string())
          if not (((rev(i)? - orig(i)?) / orig(i)?).abs() < F.from[F64](1e-2)) then
            h.log("i=" + i.string() + ", rev=" + rev(i)?.string() + " / " + orig(i)?.string() + " eps=" + ((rev(i)? - orig(i)?) / orig(i)?).abs().string())
          end
	end
	// TODO HERE
	h.assert_array_almost_eq[F, F](rev, orig, F.from[F64](1e-6), F.from[F64](0.0), "Inverse fourier_real")
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
      try
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

        // TODO DELETE
        for i in Range(0, pow2) do
          h.assert_true(a_fft(i)?.almost_eq(b_fft(i)?, F.from[F64](1e-6)), "Almost equal")
          if not a_fft(i)?.almost_eq(b_fft(i)?, F.from[F64](1e-6)) then
            h.log("i=" + i.string() + ", a_fft=" + a_fft(i)?.string() + " / b_fft=" + b_fft(i)?.string())
          end
        end
	// TODO THERE
	h.assert_array_almost_eq[Complex[F], F](a_fft, b_fft, F.from[F64](1e-6), F.from[F64](0.0), "Almost equal")

        // Now inverse
	let rev = FFT[F].fourier2(b_fft, true)

        // TODO DELETE
	for i in Range(0, pow2) do
          h.assert_true(rev(i)?.almost_eq(orig(i)?, F.from[F64](1e-6)), "Inverse fourier2, i=" + i.string())
          if not rev(i)?.almost_eq(orig(i)?, F.from[F64](1e-6)) then
            h.log("i=" + i.string() + ", rev=" + rev(i)?.string() + " / " + orig(i)?.string())
          end
	end
	// TODO HERE
	h.assert_array_almost_eq[Complex[F], F](rev, orig, F.from[F64](1e-6), F.from[F64](0.0), "Inverse fourier2")
      else
        h.fail("Equality naive fourier == fourier2")
      end
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
      try
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

        // TODO DELETE
        for i in Range(0, pow2) do
          h.assert_true(convolution(i)?.almost_eq(naive_conv(i)?, F.from[F64](1e-6)), "Almost equal convolutions")
          if not convolution(i)?.almost_eq(naive_conv(i)?, F.from[F64](1e-6)) then
            h.log("i=" + i.string() + ", convolution=" + convolution(i)?.string() + " / naive_conv=" + naive_conv(i)?.string())
          end
        end
	// TODO HERE
	h.assert_array_almost_eq[Complex[F], F](convolution, naive_conv, F.from[F64](1e-6), F.from[F64](0.0), "Almost equal convolutions")
      else
        h.fail("Equality naive convolution == convolution")
      end
    end


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
      h.assert_array_almost_eq[Complex[F], F](impulse_fft, imp, F.from[F64](1e-6), F.from[F64](0.0), "Bluestein FFT of impulse is constant")
      // TODO DELETE
      try
        for i in Range(0, n) do
          // To avoid impact of calculation errors, we approximate equality with relative tolerance
          h.assert_true(imp(i)?.almost_eq(impulse_fft(i)?, F.from[F64](1e-6)), "Bluestein FFT of impulse is constant")
          if not imp(i)?.almost_eq(impulse_fft(i)?, F.from[F64](1e-6)) then
            h.log("i=" + i.string() + ", imp=" + imp(i)?.string() + " / " + impulse_fft(i)?.string())
          end
        end
      else
        h.fail("Error while testing Bluestein FFT of impulse")
      end
      // TODO HERE

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
      h.assert_array_almost_eq[Complex[F], F](constant_fft, const, F.from[F64](1e-6), F.from[F64](0.0), "Bluestein FFT of constant is impulse")
      // TODO DELETE
      try
        for i in Range(0, n) do
          // To avoid impact of calculation errors, we approximate equality with relative tolerance
          h.assert_true(const(i)?.almost_eq(constant_fft(i)?, F.from[F64](1e-6)), "Bluestein FFT of constant is impulse")
          if not const(i)?.almost_eq(constant_fft(i)?, F.from[F64](1e-6)) then
            h.log("i=" + i.string() + ", const=" + const(i)?.string() + " / " + constant_fft(i)?.string())
          end
        end
      else
        h.fail("Error while testing Bluestein FFT of constant")
      end
      // TODO HERE

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
          h.assert_true(((x * a_fft(i)?) + (y * b_fft(i)?)).almost_eq(c_fft(i)?, F.from[F64](1e-6)), "Linear operator")
        end
      else
        h.fail("Error while testing Bluestein FFT linearity")
      end

      // Check that we have a DFT when n is a power of 2
      if (n and (n - 1)) == 0 then
        let a_orig_fft = FFT[F].fourier2(a_orig)
	// TODO DELETE
	try
	  for i in Range(0, n) do
	    h.assert_true(a_fft(i)?.almost_eq(a_orig_fft(i)?, F.from[F64](1e-6)), "Bluestein DFT of " + i.string())
	  end
	else
	  h.fail("Bluestein is not a DFT")
	end
	// TODO HERE
	h.assert_array_almost_eq[Complex[F], F](a_fft, a_orig_fft, F.from[F64](1e-6), F.from[F64](0.0), "Bluestein DFT")
      end
    end


