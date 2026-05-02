// Calculation of some constants

use "collections"

use "../../mathx/gmp"

use @mpfr_snprintf[I32](buf: Pointer[U8] tag, buf_size: U64, fmt: Pointer[U8] tag, ...)
use @pony_ctx[Pointer[None]]()
use @pony_alloc[Pointer[U8]](ctx: Pointer[None], size: USize)


actor Main
  """
  Calculates various mathematical constants using multiple precision floating point numbers.
  """
  new create(env: Env) =>
    """
    Calculate a few *well-known* mathematical constants using big floats.
    """
    env.out.print("Calculate some well-known constants using formulas and" +
                  " infinite-precision floats.")
    env.out.print("")
    env.out.print("Pi is calculated with 10_000 decimals. e calculation uses 1000 decimals")
    env.out.print("")
    e(env)
    mpfr_pi(env)
    pi_native(env)
    pi(env)


  fun mpfr_pi(env: Env) =>
    """
    The Pi constant calculated by MPFR
    """
    env.out.print("Using MPFR")
    env.out.print("----------")

    let pi_const = MPFloat.pi(10000)
    env.out.print("Pi constant from MPFR =")
    env.out.print(pi_const.string())
    env.out.print("")


  fun pi(env: Env) =>
    """
    Calculate Pi = 3.14149... ($$\pi$$) using the `MPFloat` class.

    Translation of the mathematical formulas is almost immediate. From a
    performance point of view, this version allocates more memory than
    `pi_native`.

    This algorithm uses [Borwein quadratic](https://en.wikipedia.org/wiki/Borwein%27s_algorithm#Quadratic_convergence_(1984))
    formula:

    Use the serie:
      * x_0 = sqrt(2)
      * y_0 = 0
      * pi_0 = 2 + sqrt(2)
    
      * x_n+1 = (sqrt(x_n) + 1/sqrt(x_n)) / 2
      * y_n+1 = (1 + y_n) * sqrt(x_n) / (x_n + y_n)
      * pi_n+1 = pi_n * (1 + x_n+1) * y_n+1 / (1 + y_n+1)
    """
    var i: ULong = 0
    let precision: ULong = 10000 // binary digits for the calculations, and result
    let zero = MPFloat.from[F64](0.0, precision)
    let one = MPFloat.from[F64](1.0, precision)
    let two = MPFloat.from[F64](2.0, precision)
    let two_sqrt = two.sqrt()

    let x_0 = two_sqrt
    let y_0 = zero
    let pi_0 = two + two_sqrt

    var x_n = x_0
    var y_n = y_0
    var pi_n = pi_0
    var old_pi = pi_0
    repeat
      old_pi = pi_n
      let x_n1 = (x_n.sqrt() + (one / x_n.sqrt())) / two
      let y_n1 = ((one + y_n) * x_n.sqrt()) / (x_n + y_n)
      let pi_n1 = (pi_n * (one + x_n1) * y_n1) / (one + y_n1)

      // Loop to next n
      x_n = x_n1
      y_n = y_n1
      pi_n = pi_n1
      i = i + 1
      // We stop if i > precision: there's probably a bug
    until (old_pi == pi_n) or (i > precision) end

    // Show result
    env.out.print("Using MPFloat")
    env.out.print("-------------")
    env.out.print("pi (after " + i.string() + " loops) = ")
    env.out.print(pi_n.string())
    env.out.print("")

    // Compare with pi constant
    env.out.print("Compare result with MPFloat constant")
    let pi_const = MPFloat.pi(precision)
    if pi_const == pi_n then
      env.out.print("Calculation is exact")
    else
      env.out.print("Calculation is inexact")
      let delta = pi_n - pi_const
      //env.out.print("Delta = " + delta.string())
    end
    env.out.print("")


  fun pi_native(env: Env) =>
    """
    This version of the calculation of Pi = 3.14159 ($$\pi$$) uses also the
    Borwein's quadratic serie but implement it using the `MPF` class that is
    a raw binding to [MPF](https://www.mpfr.org/) library. This class limits
    memory use by reusing variables, but the developer has to manage memory
    allocations and deallocations. Compared with `pi` implementation, this
    version feels less *natural* to mathematicians and is much longer as
    each line of code is a step of calculation in a formula. It looks more
    like machine code...

    Use the serie:
      * x_0 = sqrt(2)
      * y_0 = 0
      * pi_0 = 2 + sqrt(2)
    
      * x_n+1 = (sqrt(x_n) + 1/sqrt(x_n)) / 2
      * y_n+1 = (1 + y_n) * sqrt(x_n) / (x_n + y_n)
      * pi_n+1 = pi_n * (1 + x_n+1) * y_n+1 / (1 + y_n+1)
    """
    let prec: ILong = 10000 // Precision of calculation
    var i: ILong = 0
    var x: SMPFr = SMPFr
    var y: SMPFr = SMPFr
    var p: SMPFr = SMPFr
    var d: SMPFr = SMPFr
    var t: SMPFr = SMPFr
    var s: SMPFr = SMPFr
    var u: SMPFr = SMPFr

    // x_0 = sqrt(2)
    MPF.init2(x, prec)
    MPF.set_d(x, 2.0, RoundingNearest)
    MPF.sqrt(x, x, RoundingNearest)

    // y_0 = 0
    MPF.init2(y, prec)
    MPF.set_d(y, 0.0, RoundingNearest)

    // pi_0 = 2 + sqrt(2)
    MPF.init2(p, prec)
    MPF.set_d(p, 2.0, RoundingNearest)
    MPF.add(p, p, x, RoundingNearest)

    // Temporaries for calculations
    MPF.init2(d, prec) // Previous pi value
    MPF.init2(t, prec)
    MPF.init2(s, prec)
    MPF.init2(u, prec)

    repeat
      // Keep previous value of pi to determine when to stop
      MPF.set(d, p, RoundingNearest)

      //   x_n+1 = (sqrt(x_n) + 1/sqrt(x_n)) / 2
      MPF.sqrt(t, x, RoundingNearest)             // t = sqrt(x_n)
      MPF.d_div(s, 1.0, t, RoundingNearest)       // s = 1/sqrt(x_n)
      MPF.add(s, t, s, RoundingNearest)           // s = sqrt(x_n) + 1/sqrt(x_n)
      MPF.div_d(s, s, 2.0, RoundingNearest)       // s = s / 2 --> x_n+1

      //   y_n+1 = (1 + y_n) * sqrt(x_n) / (x_n + y_n)
      MPF.add_d(u, y, 1.0, RoundingNearest)       // u = 1 + y_n
      MPF.mul(u, u, t, RoundingNearest)           // u = (1 + y_n) * sqrt(x_n)
      MPF.add(t, x, y, RoundingNearest)           // t = x_n + y_n
      MPF.div(t, u, t, RoundingNearest)           // t = u / t --> y_n+1

      //   pi_n+1 = pi_n * (1 + x_n+1) * y_n+1 / (1 + y_n+1)
      MPF.add_d(u, s, 1.0, RoundingNearest)       // u = 1 + x_n+1
      MPF.mul(p, p, u, RoundingNearest)           // p = pi_n * (1 + x_n+1)
      MPF.mul(p, p, t, RoundingNearest)           // p = pi_n * (1 + x_n+1) * y_n+1
      MPF.add_d(u, t, 1.0, RoundingNearest)       // u = 1 + y_n+1
      MPF.div(p, p, u, RoundingNearest)           // p = p / u --> pi_n+1
      
      // Now n -> n+1
      MPF.set(x, s, RoundingNearest)              // x_n+1 = s
      MPF.set(y, t, RoundingNearest)              // y_n+1 = t

      // Stop when new value has not changed compared with old value
      i = i + 1
      // We stop if i > prec: there's a bug
    until MPF.equal_p(p, d) or (i > prec) end

    env.out.print("Using MPF")
    env.out.print("---------")
    env.out.print("Pi (after " + i.string() + " loops) = ")
    env.out.print(_string(p))
    env.out.print("")

    // Compare result with constant
    env.out.print("Compare result with MPF constant")
    MPF.const_pi(d, RoundingNearest)
    if MPF.equal_p(d, p) then
      env.out.print("Calculation is exact")
    else
      env.out.print("Calculation is inexact")
      MPF.sub(d, p, d, RoundingNearest)
      //env.out.print("Delta = " + _string(d))
    end
    env.out.print("")

    // Clean memory and cache
    MPF.clear(x)
    MPF.clear(y)
    MPF.clear(p)
    MPF.clear(d)
    MPF.clear(t)
    MPF.clear(s)
    MPF.clear(u)
    MPF.free_cache()


  fun _string(x: SMPFr): String iso^ =>
    """
    Print the result using MPFR's formatting with full precision as the `MPF`
    does not have an output function yet..

    The format string is built dynamically from the number's precision (in bits),
    converted to decimal digits via ceil(prec * log2 / log10). Without this,
    a bare `"%.RNe"` format would have precision=0 (the `.` with no digit count
    before `R` means 0 significant digits), producing output like `"3e+00"`.
    """
    let prec = MPF.get_prec(x)
    let p = ((prec.f64() * F64(2.0).log()) / F64(10.0).log()).ceil().i32()
    let fmt = ("%." + p.string() + "RNe").cstring()
    let buf_size: I32 = @mpfr_snprintf(Pointer[U8], 0, fmt, x) + 1
    let result: String iso = recover iso
      let buffer: Pointer[U8] = @pony_alloc(@pony_ctx(), buf_size.usize())
      @mpfr_snprintf(buffer, buf_size.u64(), fmt, x)
      String.copy_cstring(buffer)
    end
    consume result


  fun e(env: Env) =>
    """
    Sample calculation of e, base of Neper Logarithm using the series:
    e = 1 + 1/1! + 1/2! + ... + 1/n!
    """
    // Sample calculation of e, base of Neper Logarithm
    //
    // e = 1 + 1/1! + 1/2! + ... + 1/1000!
    //
    let one = MPFloat.from[F64](1.0, 1000)
    var res = one
    var prev_fact = one
    var old_res = res

    env.out.print("Using MPFloat")
    env.out.print("-------------")

    for i in Range(1, 1001) do
      let fact = prev_fact * MPFloat.from[F64](F64.from[USize](i), 1000)
      res = res + (one / fact)
      prev_fact = fact

      // If value does not change, the sum has converged and there's no gain
      // to continue calculation
      if res == old_res then
        env.out.print("Stopping after " + i.string() + " iterations...")
        break
      end
      old_res = res
    end
    env.out.print("e = " + res.string())
    env.out.print("")


