// Arithmetic of complex numbers

use "debug"
use "collections"

use "../assertx"
use "../pony_testx"


class val Complex[F: (Float & FloatingPoint[F]) = F64]
  is (Equatable[Complex[F]] & Stringable & Approximated[Complex[F], F])
  """
  Complex numbers done right for numeric calculations. These functions try to
  avoid some overflows, underflows or loss of precision.

  Default precision for `Complex` is `F64`.

  Most of all operations are also available using unsafe arithmetic operators.
  """
  
  let _re: F
    """
    Real part
    """


  let _im: F
    """
    Imaginary part
    """

  
  new val create(re: F = F.from[ISize](0), im: F = F.from[ISize](0)) =>
    """
    Create a new complex number from real `re` and imaginary `im`parts. If the
    real or imaginary part are not specified, they default to 0.
    Create `re + i * im`.
    """
    _re = re
    _im = im


  new val j() =>
    """
    The `i` complex in literature that is known for `i * i == -1`.

    Name *j* was selected instead of *i* for the reasons:

    - `i` is frequently used as loop variable and to prevent name clash
    - Better visual distinction from other characters like 1 in source
    - Follows some physics conventions using j instead of i.
    """
    _re = F.from[ISize](0)
    _im = F.from[ISize](1)


  new val from_polar(r: F, theta: F) =>
    """
    Create a complex number from polar form `r * exp(i * theta)`.
    `re = r * cos(theta)`, `im = r * sin(theta)`.

    The inverse is `z.abs()` for `r` and `z.arg()` for `theta`.
    """
    _re = r * theta.cos()
    _im = r * theta.sin()


  fun string(): String iso^ =>
    """
    Convert the complex to a string for display with the cartesian form
    `re + im * i`. Real and imaginary parts are printed only when non-null.

    There's no special representation for a complex number with `inf` real
    or imaginary parts. When `NaN` real or imaginary, the string is `NaN`.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    let result = match (_re, _im)
    | if (_re == zero) and (_im == zero) => "0"
    | if (_re == zero) and (_im == one) => "i"
    | if (_re == zero) and (_im == -one) => "-i"
    | if (_re == zero) => _im.string() + "i"
    | if (_im == zero) => _re.string()
    | if (_im == one) => _re.string() + "+i"
    | if (_im == -one) => _re.string() + "-i"
    | if (_im < zero) => _re.string() + _im.string() + "i"
    | if (_im > zero) => _re.string() + "+" + _im.string() + "i"
    else
      "nan"
    end
    result.clone()


  fun real(): F =>
    """
    Get the real part of a complex number.
    """
    _re


  fun imag(): F =>
    """
    Get the imaginary part of a complex number.
    """
    _im


  fun is_real(tol: F = F.epsilon()): Bool =>
    """
    `true` when the complex has no imaginary part, with a tolerance `tol`
    (default to `F.epsilon`).
    """
    ifdef debug then
      try
        Assert(tol >= F.from[ISize](0), "Precision tolerance (" + tol.string() +
              ") must be positive", true)?
      end
    end
    _im.abs() <= tol


  fun is_imag(abs_tol: F = F.epsilon()): Bool =>
    """
    `true` when the complex has no real part, with an absolute tolerance `abs_tol`
    (default to `F.epsilon)`.
    """
    ifdef debug then
      try
        Assert(abs_tol >= F.from[ISize](0), "Precision tolerance (" +
              abs_tol.string() + ") must be positive", true)?
      end
    end
    _re.abs() <= abs_tol


  fun is_null(abs_tol: F = F.epsilon()): Bool =>
    """
    `true` when the complex is equal to 0.0 with an absolute tolerance `abs_tol`
    (default to `F.epsilon)`. The tolerance is applied individually on real and
    imaginary parts and not on the modulus.
    """
    ifdef debug then
      try
        Assert(abs_tol >= F.from[ISize](0), "Precision tolerance (" +
              abs_tol.string() + ") must be positive", true)?
      end
    end
    (_re.abs() <= abs_tol) and (_im.abs() <= abs_tol)


  fun finite(): Bool =>
    """
    Return `true` when real and imaginary parts of complex are finite.
    """
    _re.finite() and _im.finite()


  fun infinite(): Bool =>
    """
    Return `false` when real or imaginary part of complex is infinite.
    """
    _re.infinite() or _im.infinite()


  fun nan(): Bool =>
    """
    Return `true` when real or imaginary part of complex is NaN (Not a Number).
    """
    _re.nan() or _im.nan()


  fun conj(): Complex[F]^ =>
    """
    Create the conjugate of the complex, `re - i * im`.
    """
    Complex[F](_re, -_im)


  fun conj_unsafe(): Complex[F]^ =>
    """
    Create the conjugate of the complex, `re - i * im` using unsafe arithmetic.
    """
    Complex[F](_re, -~_im)


  fun eq(that: Complex[F]): Bool =>
    """
    Return `true` if `this` and `that` represent the same complex number. The
    two complex numbers are equal if
    `this.real() == that.real() and this.imag() == that.imag()`
    
    Because `Complex` can't be defined as a primitive, there can be multiple
    instances of the same complex number.
    """
    (_re == that._re) and (_im == that._im)


  fun eq_unsafe(that: Complex[F]): Bool =>
    """
    Return `true` if `this` and `that` represent the same complex number, using
    unsafe arithmetic. The two complex numbers are equal if
    `this.real() ==~ that.real() and this.imag() ==~ that.imag()`
    
    Because `Complex` can't be defined as a primitive, there can be multiple
    instances of the same complex number.
    """
    (_re ==~ that._re) and (_im ==~ that._im)


  fun almost_eq(that: box->Complex[F],
                rel_tol: F = F.epsilon().sqrt(),
                abs_tol: F = F.epsilon().sqrt())
               : Bool =>
    """
    Return `true` when `this` and `that` are almost equal with a relative
    precision tolerance of `rel_tol` and absolute tolerance `abs_tol`. Both
    default to `sqrt(F.epsilon())` (~1.49e-08 for `F64`, ~3.45e-04 for
    `F32`). All tolerances are measured using the complex modulus `abs()`.

    Two complex numbers `this` and `that` are considered almost equal when
    `(this - that).abs() <= max(rel_tol * max(this.abs(), that.abs()), abs_tol)`.

    When `this` and `that` moduli are large enough, then they are compared
    using the relative tolerance. One can write that `this` is within 5% of
    `that` with `this.almost_eq(that, 0.05)`. But when `this` and `that` are
    very small, using relative comparison makes no more sense, like when one
    value is 0.0. In that case, the absolute tolerance is used to check that
    both numbers are within that range.

    A complex number that is `NaN` is never equal to another number, according
    to IEEE 754 and so is not "almost equal" either. A complex number that is
    `infinite` is only almost equal to another similar infinite number.

    This function is mostly usefull when writing tests and that rounding and
    other types of errors accumulate, and one needs to accept results within
    predetermined tolerances.

    See Knuth TAOCP 4.2.2.A for approximation of real numbers.
    """
    let zero: F = F.from[F64](0.0)
    let one: F = F.from[F64](1.0)
    ifdef debug then
      try
        Assert(rel_tol >= zero, "Relative tolerance (" + rel_tol.string() +
              ") must be positive", true)?
        Assert(rel_tol <= one, "Relative tolerance (" + rel_tol.string() +
              ") should be lower than 1.0", true)?
        Assert(abs_tol >= zero, "Absolute tolerance (" + abs_tol.string() +
              ") must be positive", true)?
      end
    end

    if nan() or that.nan() then
      false
    elseif infinite() then
      if that.infinite() then
        if (one.copysign(_re) == one.copysign(that._re)) and
          (one.copysign(_im) == one.copysign(that._im)) then
          true
        else
          false
        end
      else
        false
      end
    else
      (this - that).abs() <= (rel_tol * this.abs().max(that.abs())).max(abs_tol)
    end


  fun add(that: box->Complex[F]): Complex[F]^ =>
    """
    Addition of two complex numbers. Result is
    `(this.rea1() + that.real()) + (j() * (this.imag() + that.imag()))`.
    """
    Complex[F](_re + that._re, _im + that._im)


  fun add_unsafe(that: box->Complex[F]): Complex[F]^ =>
    """
    Usafe addition of two complex numbers using usafe arithmetic. Result is
    `(this.rea1() +~ that.real()) + (j() * (this.imag() +~ that.imag()))`.
    """
    Complex[F](_re +~ that._re, _im +~ that._im)


  fun sub(that: box->Complex[F]): Complex[F]^ =>
    """
    Substraction of two complex numbers. Result is
    `(this.rea1() - that.real()) + (j() * (this.imag() - that.imag()))`.
    """
    Complex[F](_re - that._re, _im - that._im)


  fun sub_unsafe(that: box->Complex[F]): Complex[F]^ =>
    """
    Unsafe substraction of two complex numbers using usafe arithmetic. Result is
    `(this.rea1() -~ that.real()) + (j() * (this.imag() -~ that.imag()))`.
    """
    Complex[F](_re -~ that._re, _im -~ that._im)


  fun neg(): Complex[F]^ =>
    """
    Negation of complex.
    """
    Complex[F](-_re, -_im)


  fun neg_unsafe(): Complex[F]^ =>
    """
    Negation of complex, using unsafe arithmetic.
    """
    Complex[F](-~_re, -~_im)


  fun mul(that: box->Complex[F]): Complex[F]^ =>
    """
    Multiplication of two complex numbers. Result is
    `((this.real() * that.real()) - (this.imag() * that.imag())) + (j() * ((this.imag() * that.real()) + (this.real() * that.imag()))`
    but calculation is done in less operations.
    
    This does not report an error if the calculation overflows.
    """
    let ac = _re * that._re
    let bd = _im * that._im
    let abcd = (_re + _im) * (that._re + that._im)
    Complex[F](ac - bd, abcd - ac - bd)


  fun mul_unsafe(that: box->Complex[F]): Complex[F]^ =>
    """
    Unsafe multiplication of two complex numbers using arithmetic. Result is
    `((this.real() *~ that.real()) -~ (this.imag() *~ that.imag())) + (j() * ((this.imag() *~ that.real()) +~ (this.real() *~ that.imag()))`
    but actual calculation is done in less operations.

    This does not report an error if the calculation overflows.
    """
    let ac = _re *~ that._re
    let bd = _im *~ that._im
    let abcd = (_re +~ _im) *~ (that._re +~ that._im)
    Complex[F](ac -~ bd, abcd -~ ac -~ bd)


  fun scale(s: F): Complex[F]^ =>
    """
    Multiply by the real scalar `s`. Equivalent to `mul(Complex(s, 0))` but
    avoids the overhead of a full complex multiplication.
    """
    Complex[F](_re * s, _im * s)


  fun scale_unsafe(s: F): Complex[F]^ =>
    """
    Multiply by the real scalar `s` using unsafe arithmetic.
    """
    Complex[F](_re *~ s, _im *~ s)


  fun abs(): F =>
    """
    Calculate the modulus of complex number. That is `sqrt((real() * real()) + (imag() * imag()))`,
    trying not to overflow.
    """
    let a = _re.abs()
    let b = _im.abs() 
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    if a == zero then
      b
    elseif b == zero then
      a
    elseif a >= b then
      let tmp = b / a
      a * (one + (tmp * tmp)).sqrt()
    else
      let tmp = a / b
      b * (one + (tmp * tmp)).sqrt()
    end
    
    
  fun abs_unsafe(): F =>
    """
    Calculate the unsafe modulus of complex number using unsafe arithmentic.
    That is `sqrt((real() *~ real()) +~ (imag() *~ imag()))`,
    trying not to overflow.
    """
    let a = _re.abs()
    let b = _im.abs()
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    if a ==~ zero then
      b
    elseif b ==~ zero then
      a
    elseif a >=~ b then
      let tmp = b /~ a
      a *~ (one +~ (tmp *~ tmp)).sqrt_unsafe()
    else
      let tmp = a /~ b
      b *~ (one +~ (tmp *~ tmp)).sqrt_unsafe()
    end


  fun abs2(): F =>
    """
    Returns the squared modulus `re² + im²`, avoiding the `sqrt` in `abs()`.
    Satisfies `abs2() == abs() * abs()`.

    Useful when comparing magnitudes (avoids sqrt) or in performance-critical
    inner loops where only relative magnitudes matter.
    """
    (_re * _re) + (_im * _im)


  fun abs2_unsafe(): F =>
    """
    Returns the squared modulus using unsafe arithmetic.
    """
    (_re *~ _re) +~ (_im *~ _im)
    
    
  fun invert(): Complex[F]^ =>
    """
    Inverse of complex, `1 / this`.
    
    This does not report an error if the calculation overflows.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let inf = F.from[F64](1.0 / 0.0)

    if (_re == zero) and (_im == zero) then
      Complex[F](inf, inf)
    elseif _im == zero then
      Complex[F](one / _re, zero)
    elseif _re == zero then
      Complex[F](zero, -one / _im)
    elseif _re.abs() >= _im.abs() then
      let dc = _im / _re
      let denom = _re + (_im * dc)
      Complex[F](one / denom, -dc / denom)
    else
      let cd = _re / _im
      let denom = (_re * cd) + _im
      Complex[F](cd / denom, -_re / denom)
    end
    
  
  fun invert_unsafe(): Complex[F]^ =>
    """
    Inverse of complex, `1 / this`, using unsafe arithmetic.
    
    This does not report an error if the calculation overflows.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let inf = F.from[F64](1.0 / 0.0)

    if (_re ==~ zero) and (_im ==~ zero) then
      Complex[F](inf, inf)
    elseif _im ==~ zero then
      Complex[F](one /~ _re, zero)
    elseif _re ==~ zero then
      Complex[F](zero, -one /~ _im)
    elseif _re.abs() >=~ _im.abs() then
      let dc = _im /~ _re
      let denom = _re +~ (_im *~ dc)
      Complex[F](one /~ denom, -~dc /~ denom)
    else
      let cd = _re /~ _im
      let denom = (_re *~ cd) +~ _im
      Complex[F](cd /~ denom, -~_re /~ denom)
    end
    
  
  fun div(that: box->Complex[F]): Complex[F]^ =>
    """
    Division of two complex numbers `this / that`.
    
    This does not report an error if the calculation overflows but an infinite
    complex.
    """
    let zero = F.from[ISize](0)
    let inf = F.from[F64](1.0 / 0.0)

    if (that._re == zero) and (that._im == zero) then
      Complex[F](inf.copysign(_re), inf.copysign(_im))
    elseif that._im == zero then
      Complex[F](_re / that._re, _im / that._re)
    elseif that._re == zero then
      Complex[F](_im / that._im, -_re / that._im)
    elseif that._re.abs() >= that._im.abs() then
      let dc = that._im / that._re
      let denom = that._re + (that._im * dc)
      Complex[F]((_re + (_im * dc)) / denom, (_im - (_re * dc)) / denom)
    else
      let cd = that._re / that._im
      let denom = (that._re * cd) + that._im
      Complex[F](((_re * cd) + _im) / denom, ((_im * cd) - _re) / denom)
    end
    
  
  fun div_unsafe(that: box->Complex[F]): Complex[F]^ =>
    """
    Unsafe division of two complex numbers `this /~ that`.

    It uses usafe arithmetic operations.    
    This does not report an error if the calculation overflows but returns an
    infinite complex.
    """
    let zero = F.from[ISize](0)
    let inf = F.from[F64](1.0 / 0.0)

    if (that._re ==~ zero) and (that._im ==~ zero) then
      Complex[F](inf.copysign(_re), inf.copysign(_im))
    elseif that._im ==~ zero then
      Complex[F](_re /~ that._re, _im /~ that._re)
    elseif that._re.abs() >=~ that._im.abs() then
      let dc = that._im /~ that._re
      let denom = that._re +~ (that._im *~ dc)
      Complex[F]((_re +~ (_im *~ dc)) /~ denom, (_im -~ (_re *~ dc)) /~ denom)
    else
      let cd = that._re /~ that._im
      let denom = (that._re *~ cd) +~ that._im
      Complex[F](((_re *~ cd) +~ _im) /~ denom, ((_im *~ cd) -~ _re) /~ denom)
    end
    
  
  fun sqrt(): Complex[F]^ =>
    """
    Square root of the complex.
    """
    let zero = F.from[ISize](0)
    if (_re == zero) and (_im == zero) then
      return Complex[F](zero, zero)
    end

    let c = _re.abs()
    let d = _im.abs()
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    
    let w =
      if c >= d then
        let dc = d / c
        c.sqrt() * ((one + (one + (dc * dc)).sqrt()) / two).sqrt()
      else
        let cd = c / d
        d.sqrt() * ((cd + (one + (cd * cd)).sqrt()) / two).sqrt()
      end
    
    if _re >= zero then
      Complex[F](w, _im / (two * w))
    elseif (_re < zero) and (_im >= zero) then
      Complex[F](d / (two * w), w)
    else
      Complex[F](d / (two * w), -w)
    end
    
  
  fun sqrt_unsafe(): Complex[F]^ =>
    """
    Unsafe square root of the complex, using unsafe arithmetic operations.
    """
    let zero = F.from[ISize](0)
    if (_re ==~ zero) and (_im ==~ zero) then
      return Complex[F](zero, zero)
    end

    let c = _re.abs()
    let d = _im.abs()
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    
    let w =
      if c >=~ d then
        let dc = d /~ c
        c.sqrt() *~ ((one +~ (one +~ (dc *~ dc)).sqrt()) /~ two).sqrt()
      else
        let cd = c /~ d
        d.sqrt() *~ ((cd +~ (one +~ (cd *~ cd)).sqrt()) /~ two).sqrt()
      end
    
    if _re >=~ zero then
      Complex[F](w, _im /~ (two *~ w))
    elseif (_re <~ zero) and (_im >=~ zero) then
      Complex[F](d /~ (two *~ w), w)
    else
      Complex[F](d /~ (two *~ w), -~w)
    end
    
  
  fun arg(): F =>
    """
    The argument (angle) of the complex number, in (-π, π].

    `arg(re + i*im) = atan2(im, re)`
    """
    _im.atan2(_re)
    
    
  fun exp(): Complex[F]^ =>
    """
    Exponential of complex.

    `exp(re + i*im) = e^re * (cos(im) + i*sin(im))`
    """
    let e_re = _re.exp()
    Complex[F](e_re * _im.cos(), e_re * _im.sin())
    
    
  fun exp_unsafe(): Complex[F]^ =>
    """
    Exponential of complex using unsafe arithmetic.

    `exp(re + i*im) = e^re *~ (cos(im) + i*sin(im))`
    """
    let e_re = _re.exp()
    Complex[F](e_re *~ _im.cos(), e_re *~ _im.sin())
 

  fun cos(): Complex[F]^ =>
    """
    Cosine of complex.

    `z.cos() = ((i * z).exp() + (-i * z).exp())) / 2`
    """
    let jz = j() * Complex[F](_re, _im)
    let two = Complex[F](F.from[ISize](2), F.from[ISize](0))
    (jz.exp() + ((-jz).exp())) / two


  fun sin(): Complex[F]^ =>
    """
    Sine of complex.

    `z.sin() = ((i * z).exp() - (-i * z).exp()) / (2 * i)`
    """
    let jz = j() * Complex[F](_re, _im)
    let twoi = Complex[F](F.from[ISize](0), F.from[ISize](2))
    (jz.exp() - ((-jz).exp())) / twoi


  fun tan(): Complex[F]^ =>
    """
    Tangent of complex.

    `z.tan() = z.sin() / z.cos()`
    """
    sin() / cos()


  fun asin(): Complex[F]^ =>
    """
    Inverse sine of complex.

    `asin(z) = -i * log(iz + sqrt(1 - z²))`
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, zero)
    let neg_i = Complex[F](zero, -one)
    // Explicit val copy of this to allow arithmetic (fun receiver is box)
    let z = Complex[F](_re, _im)
    // iz = i * (re + i*im) = -im + i*re
    let iz = Complex[F](-_im, _re)
    let z2 = z * z
    let w = (c_one - z2).sqrt()
    neg_i * (iz + w).log()


  fun acos(): Complex[F]^ =>
    """
    Inverse cosine of complex.

    `acos(z) = -i * log(z + i * sqrt(1 - z²))`
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, zero)
    let neg_i = Complex[F](zero, -one)
    let z = Complex[F](_re, _im)
    let z2 = z * z
    let w = (c_one - z2).sqrt()
    // i * w = (-w._im + i*w._re)
    neg_i * (z + Complex[F](-w._im, w._re)).log()


  fun atan(): Complex[F]^ =>
    """
    Inverse tangent of complex.

    `atan(z) = (i/2) * log((i + z) / (i - z))`
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let half_i = Complex[F](zero, one / F.from[ISize](2))
    let i = Complex[F](zero, one)
    let z = Complex[F](_re, _im)
    half_i * ((i + z) / (i - z)).log()


  fun cosh(): Complex[F]^ =>
    """
    Hyperbolic cosine of complex.

    `z.cosh() = (e^z + e^-z) / 2 = (i * z).cos()`
    """
    let two = Complex[F](F.from[ISize](2), F.from[ISize](0))
    (this.exp() + ((-this).exp())) / two


  fun sinh(): Complex[F]^ =>
    """
    Hyperbolic sine of complex.

    `z.sinh() = (e^z - e^-z) / 2 = (-i * (i * z).sin())`
    """
    let two = Complex[F](F.from[ISize](2), F.from[ISize](0))
    (this.exp() - ((-this).exp())) / two


  fun tanh(): Complex[F]^ =>
    """
    Hyperbolic tangent of complex.

    `z.tanh() = z.sinh() / z.cosh()`
    """
    sinh() / cosh()


  fun asinh(): Complex[F]^ =>
    """
    Inverse hyperbolic sine of complex.

    `asinh(z) = log(z + sqrt(z² + 1))`
    """
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, F.from[ISize](0))
    let z = Complex[F](_re, _im)
    let z2 = z * z
    (z + (z2 + c_one).sqrt()).log()


  fun acosh(): Complex[F]^ =>
    """
    Inverse hyperbolic cosine of complex.

    `acosh(z) = log(z + sqrt(z² - 1))`
    """
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, F.from[ISize](0))
    let z = Complex[F](_re, _im)
    let z2 = z * z
    (z + (z2 - c_one).sqrt()).log()


  fun atanh(): Complex[F]^ =>
    """
    Inverse hyperbolic tangent of complex.

    `atanh(z) = (1/2) * log((1 + z) / (1 - z))`
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, zero)
    let half = Complex[F](one / F.from[ISize](2), zero)
    let z = Complex[F](_re, _im)
    half * ((c_one + z) / (c_one - z)).log()


  fun powi(n: I32): Complex[F]^ =>
    """
    Integer power of complex.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    // Exponent special values
    if n == 0 then
      return Complex[F](one, zero)
    elseif n == 1 then
      return Complex[F](_re, _im)
    elseif n == -1 then
      return this.invert()
    end

    // Treat real and imaginary as special cases
    if _im == zero then
      return Complex[F](_re.powi(n), zero)
    elseif _re == zero then
      match n % 4
      | -3 => return Complex[F](zero, _im.powi(n))
      | -2 => return Complex[F](-_im.powi(n), zero)
      | -1 => return Complex[F](zero, -_im.powi(n))
      |  0 => return Complex[F](_im.powi(n), zero)
      |  1 => return Complex[F](zero, _im.powi(n))
      |  2 => return Complex[F](-_im.powi(n), zero)
      |  3 => return Complex[F](zero, -_im.powi(n))
      else
        // Can't happen
        return Complex[F](zero, zero)
      end
    end

    // General case
    var i = if n > 0 then n else -n end
    var result = Complex[F](one, zero)
    var x: Complex[F] = Complex[F](_re, _im)
    while i > 0 do
      if (i % 2) == 1 then
        result = result * x
      end
      x = x * x
      i = i / 2
    end

    if n > 0 then
      result
    else
      result.invert()
    end
    
    
  fun powi_unsafe(n: I32): Complex[F]^ =>
    """
    Integer power of complex using unsafe arithmetic.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    // Exponent special values
    if n == 0 then
      return Complex[F](one, zero)
    elseif n == 1 then
      return Complex[F](_re, _im)
    elseif n == -1 then
      return this.invert()
    end

    // Treat real and imaginary as special cases
    if _im ==~ zero then
      return Complex[F](_re.powi(n), zero)
    elseif _re ==~ zero then
      match n % 4
      | -3 => return Complex[F](zero, _im.powi(n))
      | -2 => return Complex[F](-_im.powi(n), zero)
      | -1 => return Complex[F](zero, -_im.powi(n))
      |  0 => return Complex[F](_im.powi(n), zero)
      |  1 => return Complex[F](zero, _im.powi(n))
      |  2 => return Complex[F](-_im.powi(n), zero)
      |  3 => return Complex[F](zero, -_im.powi(n))
      else
        // Can't happen
        return Complex[F](zero, zero)
      end
    end

    // General case
    var i = if n > 0 then n else -n end
    var result: Complex[F] = Complex[F](one, zero)
    var x: Complex[F] = Complex[F](_re, _im)
    while i > 0 do
      if (i %~ 2) == 1 then
        result = result *~ x
      end
      x = x *~ x
      i = i /~ 2
    end

    if n > 0 then
      result
    else
      result.invert_unsafe()
    end


  fun log(): Complex[F]^ =>
    """
    Natural logarithm of complex (principal value).

    `log(z) = ln|z| + i * arg(z)`
    """
    Complex[F](abs().log(), arg())


  fun log2(): Complex[F]^ =>
    """
    Logarithm base 2 of complex: `log(this) / log(2)`.

    `log2(z).real() = log2(|z|)`, `log2(z).imag() = arg(z) / ln(2)`
    """
    let ln2 = F.from[ISize](2).log()
    Complex[F](abs().log() / ln2, arg() / ln2)


  fun log10(): Complex[F]^ =>
    """
    Logarithm base 10 of complex: `log(this) / log(10)`.

    `log10(z).real() = log10(|z|)`, `log10(z).imag() = arg(z) / ln(10)`
    """
    let ln10 = F.from[ISize](10).log()
    Complex[F](abs().log() / ln10, arg() / ln10)


  fun pow(that: box->Complex[F]): Complex[F]^ =>
    """
    Complex power: `this^that = exp(that * log(this))`.
    """
    (that * log()).exp()


  fun powf(r: F): Complex[F]^ =>
    """
    Power with a real exponent: `this^r = |this|^r * exp(i * r * arg(this))`.

    More efficient than `pow(Complex(r, 0))`. Consistent with polar form:
    the result has modulus `|this|^r` and argument `r * arg(this)`.
    """
    let rho = abs().pow(r)
    let theta = arg() * r
    Complex[F](rho * theta.cos(), rho * theta.sin())


  fun dotp(that: box->Complex[F]): F =>
    """
    Scalar product of `this` and `that`,
    `(this.real() * that.real()) + (this.imag() * that.imag())`
    """
    (_re * that._re) + (_im * that._im)
  
  
  fun dotp_unsafe(that: box->Complex[F]): F =>
    """
    Unsafe scalar product of `this` and `that`,
    `this.real() * that.real() + this.imag() * that.imag()`
    """
    (_re *~ that._re) +~ (_im *~ that._im)

