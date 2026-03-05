use "collections"
use "debug"

class FLimits[F: FloatingPoint[F] val]
  """
  Display numerical limits of the current platform affecting floating-point
  calculations.

  This is somewhat equivalent to what is provided by C-language `limits.h`
  or C++ class `std::numeric_limits`.

  From [MACHAR](https://web.archive.org/web/20060904152057/http://orion.math.iastate.edu/burkardt/c_src/machar/machar.c)
  and [Numerical Recipes](http://numerical.recipes/), adapted for pony compiler.

  The name of variables and fields have been kept like in the original
  algorithm, but functions have been renamed to more descriptive ones.
  """

  let _ibeta: ISize
    """
    The radix in which floating-points numbers are represented.
    """

  var _it: ISize
    """
    Number of `ibeta`-base digits in the floating point mantissa
    """

  var _irnd: ISize
    """
    Type of rounding done by the compiler/processor.
    """

  var _negeps: ISize
    """
    `negeps` is the exponent of the smallest power of `ibeta` that, subtracted
    from 1.0, gives something different from 1.0.
    """

  let _epsneg: F
    """
    `epsneg` is `ibeta^negeps`, another way of defining floating-point precision.
    """

  var _machep: ISize
    """
    `machep` is the exponent of the smallest (most negative) power of `ibeta`
    that, added to 1.0, gives something different from 1.0.
    """

  let _eps: F
    """
    `eps` is the floating-point number `ibeta^machep`, loosely referred to as
    the “floating-point precision.”
    """

  var _ngrd: ISize
    """
    `ngrd` is the number of “guard digits” used when truncating the product of
    two mantissas to fit the representation.
    """

  var _iexp: ISize
    """
    `iexp` is the number of bits in the exponent (including its sign or bias).
    """

  var _minexp: ISize
    """
    `minexp` is the smallest (most negative) power of `ibeta` consistent with
    there being no leading zeros in the mantissa.
    """

  var _xmin: F
    """
    `xmin` is the `F` floating-point number `ibeta^minexp`, generally the smallest
    (in magnitude) useable `F` floating value.
    """

  var _maxexp: ISize
    """
    `maxexp` is the smallest (positive) power of `ibeta` that causes overflow.
    """

  var _xmax: F
    """
    `xmax` is `(1 - epsneg) * ibeta^maxexp`, generally the largest (in magnitude)
    useable `F` floating value.
    """


  new create() =>
    """
    Calculate machine limits.
    """
    let one = F.from[ISize](1)
    let two = one + one
    let zero = one - one

    var temp: F = F.from[ISize](0)
    var temp1: F = F.from[ISize](0)
    var itemp: ISize = 0

    // Calculate beta and ibeta by method of M. Malcolm
    var a: F = one
    repeat
      a = a + a
      temp = a + one
      temp1 = temp - a
    until (temp1 - one) != zero end

    var b: F = one
    repeat
      b = b + b
      temp = a + b
      itemp = (temp - a).isize()
    until itemp != 0 end

    _ibeta = itemp
    let beta: F = F.from[ISize](_ibeta)

    // Determine it and irnd
    _it = 0
    b = one
    repeat
      _it = _it + 1
      b = b * beta
      temp = b + one
      temp1 = temp - b
    until (temp1 - one) != zero end

    // Floating-point addition chops by default
    _irnd = 0
    let betah = beta / two
    temp = a + betah
    if (temp - a) != zero then
      // Floating-point addition rounds, but not in the IEEE style
      _irnd = 1
    end

    let tempa = a + beta
    temp = tempa + betah
    if (_irnd == 0) and ((temp - tempa) != zero) then
      // Floating-point addition rounds in the IEEE style
      _irnd = 2
    end

    // Determine negeps and epsneg
    _negeps = _it + 3
    let betain = one / beta
    a = one
    for _ in Range[ISize](1, _negeps + 1) do
      a = a * betain
    end
    b = a
    while true do
      temp = one - a
      if (temp - one) != zero then
        break
      end
      a = a * beta
      _negeps = _negeps - 1
    end
    _negeps = -_negeps
    _epsneg = a

    // Determine machep and eps
    _machep = -_it - 3
    a = b
    while true do
      temp = one + a
      if (temp - one) != zero then
        break
      end
      a = a * beta
      _machep = _machep + 1
    end
    _eps = a
    
    // Determine ngrd
    _ngrd = 0
    temp = one + _eps
    if (_irnd == 0) and (((temp * one) - one) != zero) then
      _ngrd = 1
    end

    // Determin iexp
    var i: ISize = 0
    var k: ISize = 1
    var y: F = betain
    var z: F = betain
    var iz: ISize = 0
    let t: F = one + _eps
    var nxres: ISize = 0
    var mx: ISize = 0

    // Loop to determine largest i such that (1 / beta) ^ (2 ^i) does not
    // underflow.
    // Loop until underflow
    while true do
      y = z
      z = y * y

      // Check for underflow
      a = z * one
      temp = z * t
      if ((a + a) == zero) or (z.abs() > y) then
        break
      end
      temp1 = temp * betain
      if (temp1 * beta) == z then
        break
      end
      i = i + 1
      k = k + k
    end

    // Determine k such that (1 / beta) ^ k does not underflow
    // First set k = 2 ^ i
    if _ibeta != 10 then
      _iexp = i + 1
      mx = k + k
    else
      // For decimal hardware
      _iexp = 2
      iz = _ibeta
      while k >= iz do
        iz = iz * _ibeta
        _iexp = _iexp + 1
      end
      mx = (iz + iz) - 1
    end

    // Determine minexp and xmin
    _xmin = y

    // Loop until underflow
    while true do
      _xmin = y
      y = y * betain
      a = y * one
      temp = y * t
      if ((a + a) != zero) and (y.abs() < _xmin) then
        k = k + 1
        temp1 = temp * betain
        if ((temp1 * beta) == y) and (temp != y) then
          nxres = 3
          _xmin = y
          break
        end
      else
        break
      end
    end
    _minexp = -k

    // Determine maxexp and xmax
    if (mx <= ((k + k) - 3)) and (_ibeta != 10) then
      mx = mx + mx
      _iexp = _iexp + 1
    end
    _maxexp = mx + _minexp

    // Adjust irnd to reflect partial underflow
    // _irnd == 3 if floating-point addition chops, and there is partial
    //            underflow
    //          4 if floating-point addition rounds, but not in the IEEE
    //            style, and there is partial underflow
    //          5 if floating-point addition rounds in the IEEE style,
    //            and there is partial underflow
    _irnd = _irnd + nxres

    // Adjust for IEEE-style machines
    if (_irnd >= 2) then
      _maxexp = _maxexp - 2
    end

    // Adjust for machines with implicit leading bit in binary mantissa, and
    // machines with radix point at extreme right of mantissa
    i = _maxexp + _minexp
    if (_ibeta == 2) and (i == 0) then
      _maxexp = _maxexp - 1
    end

    if i > 20 then
      _maxexp = _maxexp - 1
    end

    if a != y then
      _maxexp = _maxexp - 2
    end

    _xmax = one - _epsneg
    if ((_xmax * one) != _xmax) then
      _xmax = one - (beta * _epsneg)
    end
    _xmax = _xmax / (_xmin * beta * beta * beta)
    i = _maxexp + _minexp + 3
    for _ in Range[ISize](1, i + 1) do
      if _ibeta == 2 then
        _xmax = _xmax + _xmax
      else
        _xmax = _xmax * beta
      end
    end


  fun radix(): ISize =>
    """
    Return the radix for the floating-point representation of type `F`.
    Usually 2 but it can be different on specialized or legacy hardware.
    """
    _ibeta


  fun digit(): ISize =>
    """
    Return the number of `radix`-based digits in the floating-point `F` mantissa.
    """
    _it


  fun round_style(): ISize =>
    """
    Return a code describing the type of rounding done for addition.
    * 0, if floating-point addition chops.
    * 1, if floating-point addition rounds, but not in the IEEE style.
    * 2, if floating-point addition rounds in the IEEE style.
    * 3, if floating-point addition chops, and there is partial underflow.
    * 4, if floating-point addition rounds, but not in the IEEE style, and 
      there is partial underflow.
    * 5, if floating-point addition rounds in the IEEE style, and there is 
      partial underflow.
    """
    _irnd


  fun guard_digits(): ISize =>
    """
    Return the number of "guard digits" when trucating the product of two
    mantissages to fit the representation. It is

    * 0, if floating-point arithmetic rounds, or if it truncates and only 
      `digit` base `radix` digits participate in the post-normalization shift
      of the floating-point mantissa in multiplication;
    * 1, if floating-point arithmetic truncates and more than `digit` base
      `radix` digits participate in the post-normalization shift of the
      floating-point mantissa in multiplication.
    """
    _ngrd


  fun machep(): ISize =>
    """
    `machep` is the largest negative integer such that `1.0 + radix^machep != 1.0`
    except that `machep` is bounded below by `- (digit + 3)`.
    """
    _machep


  fun negeps(): ISize =>
    """
    `negeps` is the largest negative integer sunch that `1.0 - radix^negeps != 1.0`
    except that `negeps` is bounded below by `- (digit + 3)`.
    """
    _negeps


  fun exponent(): ISize =>
    """
    `exponent` is the number of bits (decimal places if `radix` = 10)
     reserved for the representation of the exponent (including the bias or
     sign) of a `F` floating-point number.
    """
    _iexp


  fun min_exponent(): ISize =>
    """
    `min_exponent` is the largest in magnitude negative integer such that
    `radix^min_exponent` is positive and normalized.
    """
    _minexp


  fun max_exponent(): ISize =>
    """
    `max_exponent` is the smallest positive power of `beta` that overflows.
    """
    _maxexp


  fun epsilon(): F =>
    """
    `epsilon` is the "floating-point precision", that is the smallest positive
    floating-point number such that `1.0 + epsilon != 1.0`.

    In particular, if either `radix == 2 or round_style == 0`, then
    `epsilon = radix^machep`. Otherwise, `epsilon = (radix^machep) / 2`.
    """
    _eps


  fun negative_epsilon(): F =>
    """
    `negative_epsilon` is a small positive floating-point number such that
    `1.0 - negative_epsilon != 1.0`.

    In particular, if `radix == 2 or round_style == 0`, then
    `negative_epsilon = radix^negeps`. Otherwise, `negative_epsilon = radix^negeps / 2`.  
    Because `negeps` is bounded below by `- (digit + 3)`, `negative_epsilon`
    might not be the smallest number that can alter 1.0 by subtraction.
    """
    _epsneg


  fun min_value(): F =>
    """
    `min_value` is the smallest non-vanishing normalized `F` floating-point power
    of the radix: `min_value = radix^min_exponent`.

    That's the smallest (in magnitude) useable floating value.
    """
    _xmin


  fun max_value(): F =>
    """
    `max_value` is the largest finite `F` floating-point number. In particular,
    `max_value = (1.0 - epsneg) * radix^max_exponent`.

    On some machines, the computed value of `xmax` will be only the second, 
    or perhaps third, largest number, being too small by 1 or 2 units in 
    the last digit of the mantissa.

    That's the largest (in magnitude) useable `F` floating value.
    """
    _xmax
