
use "../../mathx"
use "../../assertx"

use "random"
use "collections"

use "debug"

actor Main
  new create(env: Env) =>

  // Donne le résultat correct: 1, 2, 3, 4, 5, 6, 7, 8, 7, 6, 5, 4, 4, 4, 3, 2
    let a = Array[F64](16)
    a.push(1.0)
    a.push(-2.0)
    a.push(-2.0)
    a.push(4.0)
    a.push(4.0)
    a.push(4.0)
    a.push(4.0)
    a.push(-8.0)
    a.push(-8.0)
    a.push(-8.0)
    a.push(-8.0)
    a.push(-8.0)
    a.push(-8.0)
    a.push(-8.0)
    a.push(-8.0)
    a.push(16.0)
/* */
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(16.0)
    a.push(-32.0) // 32
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
    a.push(-32.0)
/* */
    let b = Array[Complex[F64]](16)
    b.push(Complex[F64](1.0))
    b.push(Complex[F64](-2.0))
    b.push(Complex[F64](-2.0))
    b.push(Complex[F64](4.0))
    b.push(Complex[F64](4.0))
    b.push(Complex[F64](4.0))
    b.push(Complex[F64](4.0))
    b.push(Complex[F64](-8.0))
    b.push(Complex[F64](-8.0))
    b.push(Complex[F64](-8.0))
    b.push(Complex[F64](-8.0))
    b.push(Complex[F64](-8.0))
    b.push(Complex[F64](-8.0))
    b.push(Complex[F64](-8.0))
    b.push(Complex[F64](-8.0))
    b.push(Complex[F64](16.0))
/* */
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](16.0))
    b.push(Complex[F64](-32.0)) // 32
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
    b.push(Complex[F64](-32.0))
/* */
    env.out.print("orig a =" + dump_complex(a))
    //FFT[F64].fourier_complex(a)
    FFT[F64].fourier_real(a) // WRONG!
    env.out.print("real a =" + dump_complex(a))
    FFT[F64].fourier_real(a, true)
    env.out.print("inv a  =" + dump_complex(a))

//    env.out.print("naiv b =" + dump(naive_fft[F64](b)))
    FFT[F64].fourier(b)
    try
      var res: String = "target ="
      for i in Range(0, b.size() / 2) do
        if i == 0 then
          res = res + "; " + b(i)?.real().string() + "; " + b(b.size() / 2)?.real().string()
        else
          res = res + "; " + b(i)?.real().string() + "; " + b(i)?.imag().string()
        end
      end
      env.out.print(res)
    end

    //env.out.print("orig b =" + dump(b))
    //FFT[F64].fourier_complex(b)
    //FFT[F64].fourier(b)
//    env.out.print("fft b  =" + dump(b))
//    env.out.print("ninv b =" + dump(naive_fft[F64](b, true)))

/*
    let a: Array[Complex[F64]] = []

    a.push(Complex[F64](1.0))
    a.push(Complex[F64](1.0))
    a.push(Complex[F64](1.0))
    a.push(Complex[F64](1.0))
    a.push(Complex[F64](0.0))
    a.push(Complex[F64](0.0))
    a.push(Complex[F64](0.0))
    a.push(Complex[F64](0.0))

    a.push(Complex[F64](1.0))
    a.push(Complex[F64](1.0))
    a.push(Complex[F64](1.0))
    a.push(Complex[F64](1.0))
    a.push(Complex[F64](0.0))
    a.push(Complex[F64](0.0))
    a.push(Complex[F64](0.0))
    a.push(Complex[F64](0.0))

    let b: Array[Complex[F64]] = []

    b.push(Complex[F64](1.0))
    b.push(Complex[F64](1.0))
    b.push(Complex[F64](1.0))
    b.push(Complex[F64](1.0))
    b.push(Complex[F64](0.0))
    b.push(Complex[F64](0.0))
    b.push(Complex[F64](0.0))
    b.push(Complex[F64](0.0))

    b.push(Complex[F64](1.0))
    b.push(Complex[F64](1.0))
    b.push(Complex[F64](1.0))
    b.push(Complex[F64](1.0))
    b.push(Complex[F64](0.0))
    b.push(Complex[F64](0.0))
    b.push(Complex[F64](0.0))
    b.push(Complex[F64](0.0))

    a.push(Complex[F64](0.0))
    a.push(Complex[F64](1.0))
    a.push(Complex[F64](2.0))
    a.push(Complex[F64](3.0))
    a.push(Complex[F64](4.0))
    a.push(Complex[F64](5.0))
    a.push(Complex[F64](6.0))
    a.push(Complex[F64](7.0))
    a.push(Complex[F64](8.0))
    a.push(Complex[F64](9.0))
    a.push(Complex[F64](10.0))
    a.push(Complex[F64](11.0))
    a.push(Complex[F64](12.0))
    a.push(Complex[F64](13.0))
    a.push(Complex[F64](14.0))
    a.push(Complex[F64](15.0))

    env.out.print("a=" + dump(a))
    //FFT.rearrange[Complex[F64]](FFT.rearrange[Complex[F64]](a))

    //env.out.print("b=" + dump(b))
    //FFT.rearrange_bis[Complex[F64]](FFT.rearrange_bis[Complex[F64]](b))

    let reference = FFT.naive_fft(a)
    env.out.print("Reference naive       =" + dump(reference))
    let inverse = FFT.naive_fft(reference, true)
    env.out.print("Inverse naive         =" + dump(inverse))

    FFT.fourier(a)
    env.out.print("fourrier(a)           =" + dump(a))
    FFT.fourier(a, true) // Inverse
    env.out.print("fourier-1(fourier(a)  =" + dump(a))

    FFT.fourier2(b)
    env.out.print("fourrier2(b)          =" + dump(b))
    FFT.fourier2(b, true) // Inverse
    env.out.print("fourier2-1(fourier2(b)=" + dump(b))

*/
/*
    let a': Array[F64] iso = []

    a'.push(1.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)
    a'.push(0.0)

    //let a'': Array[F64] iso = []
    let a'': Array[F64] = []
    a''.push(1.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)

    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)
    a''.push(0.0)

    a''.push(0.0)
    a''.push(1.0)
    a''.push(2.0)
    a''.push(3.0)
    a''.push(4.0)
    a''.push(5.0)
    a''.push(6.0)
    a''.push(7.0)
    a''.push(8.0)
    a''.push(9.0)
    a''.push(10.0)
    a''.push(11.0)
    a''.push(12.0)
    a''.push(13.0)
    a''.push(14.0)
    a''.push(15.0)

    //var b' = FFT.fourier_complex(consume a')
    //b' = FFT.fourier_complex(consume b', true)
    //env.out.print("Fourier=" + dump_complex(consume b'))

/*
    var c = FFT.fourier_real(consume a'')
    //c = FFT.fourier_real(consume c, true)
    env.out.print("Fourier_real=" + dump_complex(consume c))

    try
      //let m = MPInt("42572974527947842435")?
      //let n = MPInt("8724357928742796304")?
      let m = MPInt("425141")?
      let n = MPInt("58765")?

      //env.out.print("fast m*n=" + m.fast_mul(n).string())
      env.out.print("mul  m*n=" + (m * n).string())

      env.out.print("fast m*n=" + m.fast_mul(n).dump())
      env.out.print("mul  m*n=" + (m * n).dump())
    else
      env.out.print("Can't create m")
    end

*/
*/
  fun dump(a: Array[Complex[F64]]): String =>
    let size = a.size()
    var i: USize = 0
    var r: String = ""
    while i < size do
      r = r + "; " + try a(i)?.string() else "***" end
      i = i + 1
    end
    r

  fun dump_complex(a: Array[F64 val]): String =>
    let size = a.size()
    var i: USize = 0
    var r: String = ""
    while i < size do
      r = r + "; " + try a(i)?.string() else "***" end
      i = i + 1
    end
    r

  fun naive_fft[F: F64](a: Array[Complex[F]], inverse: Bool = false): Array[Complex[F]]^ =>
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

