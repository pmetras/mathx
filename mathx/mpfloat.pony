// Multi-precision floating point numbers

use "debug"
use "../assertx"


class MPFloat
  """
  MPFloat represent real numbers with arbitrary precision.

  https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic

  Contrarily to fixed precision floating point numbers, MPFloat supports infinite
  precision limited only by available memory and the needs of the user. Digits
  are encoded in base 256 in a `Array[U8]`.

  The multiplication of two `MPFloat` uses a fast-multiplication algorithm using
  real Fast Fourier Transform (https://en.wikipedia.org/wiki/Sch%C3%B6nhage%E2%80%93Strassen_algorithm),
  that allows floats with around 10e6 decimals. If more decimals were required,
  one should used Number Therotic Transforms
  [NTT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform_over_a_ring#Number-theoretic_transform)
  instead. If the system detects loss of accuracy because of too many base-256
  digits being used in multiplication-derived operation, assertion will raise.

  TODO: In this version of MPFloat, results of operations create new instances
  of `MPFloat` and keeps the accuracy of operands. I must see if memory allocation
  can be optimized to limit memory allocation, for instance with pre-allocating
  the resulting `MPFloat` and using it for all operations.
  """

  let _digits: Array[U8]
    """
    The arbitrary-precision number encoded in base 256.
    """


  var _exponent: U64 = 0
    """
    The exponent in base 256.
    """


  new create(size: USize = 0) =>
    """
    Create a new arbitrary-precision number with a capacity of `size` 256-base-
    encoded digits (by default 0, giving 0 10-base digits of accuracy).
    """
    _digits = Array[U8].init(0, size)


  new from_string(value: String, size: USize = 0) ? =>
    """
    Create a new `MPFloat` initialized from the string `value`. An error is
    raised if the `value` does not represent a valid floating point number.
    The accuracy of the float is set by the maximum of the accuracy of the
    specified `value` or from the `size` parameter.
    """
    // TODO
    _digits = Array[U8](size)
    error


  new pi(size: USize = 64) =>
    """
    The pi constant, calculated with the specified `size` accuracy.

    The value of pi is calculated by the following serie, with pi being the
    limit of `pi_n`:

    ```
    x_0 = sqrt(2)
    pi_0 = 2 + sqrt(2)
    y_0 = sqrt(sqrt(2))

    x_n+1 = (sqrt(x_n) + 1/sqrt(x_n)) / 2
    pi_n+1 = pi_n * ((x_n+1 + 1) / (y_n + 1))
    y_n+1 = ((y_n * sqrt(x_n+1)) + 1/sqrt(x_n+1)) / (y_n + 1)
    """
    var i: USize = 0
    var pi_n = MPFloat(size + 2)

    try
      // two = 2
      let two = MPFloat(size + 2)
      two._digits.update(0, 2)?
      i = 1
      while i < (size + 2) do
        two._digits.update(i, 0)?
        i = i + 1
      end

      // x_0 = sqrt(2)
      var x_n = two.sqrt()

      // pi_0 = 2 + sqrt(2)
      pi_n = (two + x_n) << 1

      // y_0 = sqrt(sqrt(2))
      var y_n = x_n.sqrt()

      // Now recurrence
      while true do
        // x_n+1 = (sqrt(x_n) + 1/sqrt(x_n)) / 2
        let s = x_n.sqrt()
        (let t, _) = (s + s.inv())._short_div(2)
        x_n = t << 1

        // y_n+1 = ((y_n * sqrt(x_n+1)) + 1/sqrt(x_n+1)) / (y_n + 1)
        let u = x_n.sqrt() // Now x_n is in fact x_n+1
        let v = (((y_n * u) << 1) + u.inv()) << 1

        // y_n = y_n + 1
        y_n._digits.update(0, y_n._digits(0)? + 1)?
        (y_n, _) = v / y_n
        y_n = y_n << 1

        // pi_n+1 = pi_n * ((x_n+1 + 1) / (y_n + 1))
        // x_n+1 = x_n+1 + 1
        x_n._digits.update(0, x_n._digits(0)? + 1)?
        (let w, _) = x_n / y_n
        pi_n = (pi_n * (w << 1)) << 1

        // Check convergence
        let m = w._digits(0)? - 1
        i = 1
        while i < (size + 1) do
          if w._digits(i)? != m then
            break
          end
          i = i + 1
        end

        // We've reached desired accuracy
        if i == (size + 1) then
          break
        end
      end
    else
      Debug("MPFloat.pi: Index out of bounds [0..~" + (size + 2).string() +
            ")! i=" + i.string())
      // TODO: Should exit
    end
    _digits = pi_n._digits


  fun _size(): USize =>
    """
    Gives the size of the internal representation of the floating point number.
    """
    _digits.size()


  fun _lowb(a: U16): U8 =>
    """
    Get the lower byte of `a`.
    """
    a.u8()


  fun _highb(a: U16): U8 =>
    """
    Get the high byte of `a`.
    """
    a.shr(8).u8()


  fun _addc(a: U8, b: U8, c: U16): (U8, U16) =>
    """
    Addition of `a + b` with carry over `c` on a single byte. Return the result
    and the new carry over.
    """
    let r = a.u16() + b.u16() + c
    (r.u8(), r.shr(8))


  fun add(that: MPFloat): MPFloat =>
    """
    Calculate `this + that`. The result is a new allocated `MPFloat` whose size
    is the minimum size of `this` and `that` plus one, to store it without lost
    of accuracy.
    """
    // The size of the result is the size of the shorter float + 1 digit
    // for the carry.
    let s = _size().min(that._size())
    let size = s + 1
    let res = MPFloat(size)
    var i = size - 1

    try
      var carry: U16 = 0
      repeat
        i = i - 1
        (let sum, carry) = _addc(_digits(i)?, that._digits(i)?, carry)
        res._digits.update(i + 1, sum)?
      until i == 0 end
      res._digits.update(0, _lowb(carry))?
    else
      Debug("MPFloat.add: Index out of bounds [0.." + (s + 1).string() +
            ")! i=" + i.string())
      // TODO: We should exit!
    end
    res


  fun _subc(a: U8, b: U8, c: U16): (U8, U16) =>
    """
    Substraction `a - b` with carry `c`.
    """
    let r = ((a.max_value().u16() + a.u16()) - b.u16()) + c.shr(8)
    (r.u8(), r.shr(8))


  fun sub(that: MPFloat): (MPFloat, I8) =>
    """
    Calculate `this - that`, considering `this` and `that` as unsigned
    256-base numbers. If the result is negative, the second part of the
    tuple is `-1`.
    """
    // The size of the result must be as long as the shorter float.
    let s = _size().min(that._size())
    let size = s
    let res = MPFloat(size)
    var i = size
    var negat: I8 = 1

    try
      let max = _digits(0)?.max_value()
      var carry: U16 = max.u16() + 1
      repeat
        i = i - 1
        (let substract, carry) = _subc(_digits(i)?, that._digits(i)?, carry)
        res._digits.update(i, substract)?
      until i == 0 end
      negat = _highb(carry).i8() - 1
    else
      Debug("MPFloat.sub: index out of bounds [0.." + s.string() +
            ")! i=" + i.string())
      // TODO: We should exit!
    end

    // Result and overflow
    (res, negat)


  fun _short_add(b: U8): MPFloat =>
    """
    Short addition of the single word `b`, added to the least significant digit
    of `this`. To ensure that the result does not require two digits before the
    radix point, one may right-shift to get a 0 first digit and keep track of
    the shift separately.
    """
    // Result is as large as `this`
    let size = _size()
    let res = MPFloat(size)
    var i = size

    try
      var carry: U16 = b.u16().shl(8)
      repeat
        i = i - 1
        carry = _digits(i)?.u16() + _highb(carry).u16()
        if (i + 1) < size then
          res._digits.update(i + 1, _lowb(carry))?
        end
      until i == 0 end
      res._digits.update(0, _highb(carry))?
    else
      Debug("MPFloat._short_add: Index out of bounds [0.." + size.string() +
            ")! i=" + i.string())
      // TODO: We should exit!
    end
    res


  fun _short_mul(b: U8): MPFloat =>
    """
    Short multiplication of `this` by the single word `b`. To ensure that the
    result does not require two digits before the radix point, one may
    right-shift to get a 0 first digit and keep track of the shift separately.
    """
    // Result is as large as `this`
    let size = _size()
    let res = MPFloat
    var i = size

    try
      var carry: U16 = 0
      repeat
        i = i - 1
        carry = (_digits(i)?.u16() * b.u16()) + _highb(carry).u16()
        if i < (size - 1) then
          res._digits.update(i + 1, _lowb(carry))?
        end
      until i == 0 end
      res._digits.update(0, _highb(carry))?
    else
      Debug("MPFloat._short_mul: Index out of bounds [0.." + size.string() +
            ")! i=" + i.string())
      // TODO: We should exit!
    end
    res


  fun _short_div(b: U8): (MPFloat, U8) =>
    """
    Short division of `this` by `b`. The second part of the result tuple is the
    remainder.
    """
    // Result is as large as `this`
    let size = _size()
    let res = MPFloat(size)
    var i: USize = 0
    var remain: U16 = 0

    try
      while i < size do
        (let quotient, remain) = (remain.shl(8) + _digits(i)?.u16()).divrem(b.u16())
        res._digits.update(i, _lowb(quotient))?
        i = i + 1
      end
    else
      Debug("MPFloat._short_div: Index out of bounds [0.." + size.string() +
            ")! i=" + i.string())
      // TODO: We should exit!
    end
    (res, remain.u8())


  fun neg(): MPFloat =>
    """
    One-complement negate the current float.
    """
    let size = _size()
    let res = MPFloat(size)
    var i = size

    try
      let max = _digits(0)?.max_value().u16()
      var carry: U16 = max + 1
      repeat
        i = i - 1
        carry = (max - _digits(i)?.u16()) + _highb(carry).u16()
        res._digits.update(i, _lowb(carry))?
      until i == 0 end
    else
      Debug("MPFloat.neg: Index out of bounds [0.." + size.string() +
            ")! i=" + i.string())
      // TODO: We should exit!
    end
    res


  fun clone(): MPFloat =>
    """
    Create of copy of the current float.
    """
    let size = _size()
    let res = MPFloat(size)
    _digits.copy_to(res._digits, 0, 0, size)
    res


  fun shl(n: USize = 1): MPFloat =>
    """
    Left shift current float by `n` digits. The final elements in the underlying
    digits array are set to 0. By default, it left-shifts by 1 digit.
    """
    let size = _size()
    let s = size - n
    let res = MPFloat(size)
    var i: USize = 0

    try
      while i < s do
        res._digits.update(i, _digits(i + n)?)?
        i = i + 1
      end
      while i < size do
        res._digits.update(i, 0)?
        i = i + 1
      end
    else
      Debug("MPFloat.shl: Index out of bounds [0.." + size.string() +
            ")! i=" + i.string())
      // TODO: We should exit!
    end
    res


  fun mul(that: MPFloat): MPFloat =>
    """
    Fast multiplication of `this` and `that` using FFT. The result has a size
    of `this.size() + that.size() - 1`.
    """
    let this_size = _size()
    let that_size = that._size()
    let max_size = this_size.max(that_size)
    let size = (this_size + that_size) - 1
    // The minimal power of 2 for the Fourier transforms
    let pow2 = (max_size + 1).next_pow2()

    let res = MPFloat(size)
    var i: USize = 0

    try
      // Convert to arrays of F64
      let a = Array[F64].init(0.0, pow2)
      i = 0
      while i < this_size do
        a.update(i, _digits(i)?.f64())?
        i = i + 1
      end

      let b = Array[F64].init(0.0, pow2)
      i = 0
      while i < that_size do
        b.update(i, that._digits(i)?.f64())?
        i = i + 1
      end

      // Convolution
      FFT.fourier_real(a)
      FFT.fourier_real(b)

      // Multiply in dual space with complex multiplication
      b.update(0, b(0)? * a(0)?)?
      b.update(1, b(1)? * a(1)?)?
      i = 2
      while i < max_size do
        let temp = b(i)?
        b.update(i, (temp * a(i)?) - (b(i + 1)? * a(i + 1)?))?
        b.update(i + 1, (temp * a(i + 1)?) + (b(i + 1)? * a(i)?))?
        i = i + 2
      end

      // Inverse FFT
      FFT.fourier_real(b, true)

      // Convert back to U8
      var carry: F64 = 0.0
      i = pow2
      let s2 = (max_size / 2).f64()
      let base: F64 = 256.0
      repeat
        i = i - 1
        let temp = (b(i)? / s2) + carry + 0.5
        carry = (temp / base).u64().f64()
        b.update(i, temp - (carry * base))?

        ifdef debug then
          Assert(carry < base, "We can't have a carry higher than 256! carry=" +
                carry.string(), true)?
        end
      until i == 0 end

      // Create result
      res._digits.update(0, carry.u8())?
      i = 1
      let sp = size.min(pow2)
      while i < sp do
        res._digits.update(i, b(i - 1)?.u8())?
        i = i + 1
      end
    else
      Debug("MPFloat.mul: Index out of bounds [0.." + size.string() +
            ")! i=" + i.string())
      // TODO: Should exit
    end
    res


  fun inv(): MPFloat =>
    """
    Calculate `1 / this`.

    This is done with [Newton approximation](https://en.wikipedia.org/wiki/Multiplicative_inverse#Algorithms)
    of serie `x_n+1 = x_n * (2 - this * x_n)`.

    `this._digits` are interpreted as a radix-256 number with the radix point
    after `this._digits(0)`. The result is set to the most signiﬁcant digits
    of its reciprocal, with the radix point after `res._digits(0)`.
    """
    // We require at least twice digits of accuracy
    let size = (2 * _size()) + 1
    var res = MPFloat(size)
    var i: USize = 0

    try
      // Get the F64 corresponding to the most significant part of `this`
      // (4 bytes are enough)
      i = if size > 4 then size.min(4) else size end
      var ftemp: F64 = 0.0
      repeat
        i = i - 1
        ftemp = (ftemp / 256.0) + _digits(i)?.f64()
      until i == 0 end

      ftemp = 1.0 / ftemp

      // Now initialise `res` with the F64-estimate of `1/this`
      i = 0
      while i < size do
        res._digits.update(i, ftemp.u8())?
        ftemp = 256.0 * (ftemp - res._digits(i)?.f64())
        i = i + 1
      end

      // Apply Newton method
      while true do
        // temp = 2 - res * this
        // TODO: See if clone() can be removed by optimizing MPFloat capabilty
        var temp: MPFloat = res * this.clone()
        temp = -(temp << 1)
        temp._digits.update(0, temp._digits(0)? + 2)?
        // res = temp * res
        res = temp * res
        res = res << 1

        // If fractional part of `temp` is not 0, `temp` has not converged to 1
        i = 1
        while i < (size - 1) do
          if temp._digits(i)? != 0 then
            break
          end
          i = i + 1
        end

        // We've reached the required accuracy
        if i == (size - 1) then
          break
        end
      end
    else
      Debug("MPFloat.inv: Index out of bounds [0.." + size.string() +
            ")! i=" + i.string())
      // TODO: We should exit
    end
    res


  fun div(that: MPFloat): (MPFloat, MPFloat) =>
    """
    Calculate `this / that` and returns the result of the division and its
    remainder.

    Accuracy of divisor `that` must be less or equal to accuracy of dividend
    `this`. The quotient has an accuracy equal to the accuracy of dividend minus
    accuracy of divisor, while remainder has accuracy of divisor.

    Note that this function should be named `divrem` to remain consistent with
    usage of Pony integer arithmetic, but has be name `div` to allow usage of
    alias `/`.
    """
    let this_size = _size()
    let that_size = that._size()
    let off = this_size - that_size

    var res = MPFloat(0) // MPFloat(off + 1)
    var remain = MPFloat(0) // MPFloat(that_size)
    var i: USize = 0
    try
      ifdef debug then
        Assert(that._size() > _size(), "Divisor longer (" + that._size().string() +
              ") than dividend (" + _size().string() + ")", true)?
      end

      // temp = 1/that
      var temp = that.inv()
      // TODO: See how to remove clone()
      temp = temp * this.clone()
      temp = temp._short_add(1)
      res = temp << 2
      // TODO: Remove clone()
      remain = res * this.clone()
      remain = remain << 1
      // TODO: Remove clone
      (remain, let negat) = remain - this.clone()

      ifdef debug then
        Assert(negat == 0, "Accuracy too small! negat=" + negat.string(), true)?
      end

      if off > 0 then
        i = 0
        while i < that_size do
          remain._digits.update(i, remain._digits(i + off)?)?
          i = i + 1
        end
      end
    else
      Debug("MPFloat.div: Index out of bounds [0.." + that_size.string() +
            ")! i=" + i.string())
      // TODO: We should exit
    end
    (res, remain)


  fun sqrt(): MPFloat =>
    """
    Calculate the square root of the current float using
    [Newton method](https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Iterative_methods_for_reciprocal_square_roots).

    We use the serie `x_n+1 = x_n * (3 - this * x_n^2) / 2` that converges
    quadratically to `1 / sqrt(this)`, and then take the inverse.
    """
    let size = _size()

    var res1 = MPFloat(size)
    var i: USize = 0
    try
      // Get the F64 corresponding to the most significant part of `this`
      // (4 bytes are enough)
      i = if size > 4 then size.min(4) else size end
      var ftemp: F64 = 0.0
      repeat
        i = i - 1
        ftemp = (ftemp / 256.0) + _digits(i)?.f64()
      until i == 0 end

      ftemp = 1.0 / ftemp.sqrt()

      // Now initialise `res1` with the F64-estimate of `1/this.sqrt()`
      i = 0
      while i < size do
        res1._digits.update(i, ftemp.u8())?
        ftemp = 256.0 * (ftemp - res1._digits(i)?.f64())
        i = i + 1
      end

      // Apply Newton method
      while true do
        res1 = (res1 * res1) << 1
        // TODO: remove clone()
        res1 = res1 * this.clone()
        res1 = res1 << 2
        res1._digits.update(0, res1._digits(0)? + 3)?
        (res1, _) = res1._short_div(2)

        i = 1
        while i < (size - 1) do
          if res1._digits(i)? != 0 then
            res1 = (res1 * res1) << 1
            break
          end
          i = i + 1
        end

        if i == (size - 1) then
          break
        end
      end
    else
      Debug("MPFloat.sqrt: Index of of bounds [0.." + size.string() +
            ")! i=" + i.string())
      // TODO: We should exit
    end
    // TODO: remove clone()
    let res = (res1 * this.clone()) << 1
    res
    

  fun string(): String ref =>
    """
    Return the decimal representation of the float.

    TODO: Manage sign
    """
    let size = _size()
    // The case of the non-initialized MPFloat
    if size == 0 then
      return "0.0".clone()
    end

    var res = String(3 * size)
    var i: USize = 0
    try
      res.add(_digits(0)?.string())
      res.add(".")

      var temp = this << 1
      i = 0
      while i < size do
        temp = temp._short_mul(10)
        res.add(temp._digits(0)?.string())

        if (i % 3) == 0 then
          res.add("_")
        end

        temp = temp << 1
        i = i + 1
      end
    else
      Debug("MPFloat.string: Index out of bounds [0.." + size.string() +
            ")! i=" + i.string())
      // TODO: We should exit
    end
    res


  
