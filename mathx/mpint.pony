// Multi-precision integers

use "debug"
use "../assertx"
use "../bitsx"
use "format"
use "collections"


class val MPInt is (SignedInteger[MPInt, MPInt] & UnsignedInteger[MPInt] & Comparable[MPInt] & Stringable)
  """
  A multiple-precision integer that can be used for arbitrary precision
  calculation on integers.

  MPInt is a `class val`, meaning it is globally immutable. Internal 
  calculations use local mutable arrays within `recover` blocks for efficiency.
  https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic

  This pure Pony implementation must not be used for cryptographic applications.

  A `MPInt` is represented by an array of `U16`, with the least important digit
  in index `0`. The digits are encoded in base `65536`. Because multiplication
  and division use a Fast Fourier Transform algorithm using `F64` values having
  a 53 bits mantissa, this limits the size of `MPInt` (and `MPFloat`) to around
  1e6 base-65536 digits. This should be enough for most uses of these classes.
  If you need more precision, use the GMP/MPFR implementations of `MPInt` and
  `MPFloat`.
  """

  let _digits: Array[U16] val
    """
    The arbitrary-precision number encoded in base 65536 (see [`_base`](#_base)).
    The least important digit is in index `0` (Little-Endian order).

    Digits are encoded in reverse of traditional occidental practice, with
    the least important first. So number 1234 is encoded as `4321`, with
    `_digits(0)? == 4, _digits(1)? == 3, _digits(2)? == 2, _digits(3) == 1`.
    This is the natural way for many integer operations, where having a carry
    needs only to `push` it to the array.

    Why select an array of `U16` instead of `U64` or `U32` that would be best
    fitted for modern CPU architectures? The size of a **digit** is related to
    the size of **flint**, representing integers as floats... I want to use
    a fast multiplication algorithm. Mutliplication of arbitrary-long integers
    involves lots of sum of products:
    `(_digits(i)? * _digits(j - k)?) + (_digits(i + 1) * _digits(j - k - 1)?) + ...` 
    These are convolutions, that can be calculated efficiently by using FFT
    (Fast Fourier Transforms). So the multiplication is transformed to a dual
    space with FFT and converted back with inverse FFT. But the result must be
    an integer, while the FFT is calculated using floats. In order not to
    lose precision (and to be exact), we should use a range of floats that can
    represent integers (**flint**). `F64` in Pony is using iEEE-754 format and
    has a mantissa of 53 bits. Without taking into account the summation, a
    product of integers must fit into 53 bits, meaning that the integers must
    be less than 26 bits. We can use `U16` or `U8`... I've selected `U16`.
    If the implementation and algorithms are changed, one can consider using a
    larger **digit**, particularly if using cyclic convolutions (see Knuth TAOCP
    4.6.4.E59).

    We assume that `_digits` contains at least one **digit** item. The zero
    values is coded with `_digits(0) == 0`.
    """


  let _negative: Bool
    """
    The sign of the integer, `true` when less than `0`.
    """


  fun tag _base(): U32 =>
    """
    The base used to encode the 'digits' of `MPInt` numbers. The digits are
    coded as `U16` giving a `2^16` base.
    """
    65536


  fun tag _base_bits(): USize =>
    """
    The number of bits in a `MPInt` 'digit'.
    """
    16


  new val _create(negative: Bool, digits: Array[U16] val) =>
    """
    Private constructor for sharing digits arrays.
    """
    _digits = digits
    // -0 is positive... This assumes that digits is normalized.
    _negative = if (digits.size() == 1) and (try digits(0)? == 0 else false end) then
      false
    else
      negative
    end


  new val from_string(s: String = "", base: U32 = 10) ? =>
    """
    Create a new arbitrary-precision integer interpreting the string `s` as
    a number in base `base` (default decimal). An error is raised if `s` is
    not a valid number. An empty `s` string (the default) create a `MPInt`
    of value 0.

    The base can be in the range `[2..36]`. In case of a base larger than `10`,
    the case of the digits is not important (i.e. `a` and `A` represent the
    same digit).

    Because `MPInt` string representation can have a large number of digits, an
    exponential format representation is supported. The exponent is always in
    decimal and must be positive. For base larger than 10, the exponent
    character is `'@`'. For instance, you must write `123@987` instead of
    `123e987` if you are using base 16. For bases lower than 10, one can use
    `'e'`, `'E'` or `'@'` to indicate the exponent part of the integer. The
    value of the exponent is an `ILong`. Note that even though the value of
    the exponent is a decimal number, the exponent is in `base` base. For
    instance, `1@3` in hexadecimal represent the integer with the value
    `1 * 16 * 16 * 16` in decimal.

    The string format uses the following grammar:

    ```BNF
    mpint ::= sign? (digit | separator)+ exponent?
    sign ::= '+' | '-'
    digit ::= '0' | '1' | ... | '9' | 'a' | ... | 'z' | 'A' | ... | 'Z'
    separator ::= '_'
    exponent ::= [ 'e' | 'E' | '@' ] '+'? digit+
    ```

    The whole string `s` must represent an integer number. Leading spaces
    are accepted but not trailing spaces. Note that the syntax of the `"mpint"`
    is more flexible that the one accepted by Pony's parser, for instance
    accepting `+e86` as a valid integer or when using consecutive separators.
    Don't abuse it!
    """
    ifdef debug then
      Assert((2 <= base) and (base <= 36), "The base must be in [2..36]", true)?
    end

    if s.size() == 0 then
      _negative = false
      _digits = [0]
      return
    end

    var i: USize = 0
    let size = s.size()
    
    // Skip leading spaces
    while (i < size) and (s(i)? == ' ') do
      i = i + 1
    end

    if i >= size then
      error
    end

    var is_neg = false
    if s(i)? == '-' then
      is_neg = true
      i = i + 1
    elseif s(i)? == '+' then
      i = i + 1
    end

    if (i < size) and (s(i)? == ' ') then
      error
    end // No space after sign

    // Check for exponent
    var iexp: ISize = -1
    if base > 10 then
      iexp = try s.find("@", i.isize())? else iexp end
    else
      iexp = try s.find("e", i.isize())? else iexp end
      if iexp < 0 then
        iexp = try s.find("E", i.isize())? else iexp end
      end
      if iexp < 0 then
        iexp = try s.find("@", i.isize())? else iexp end
      end
    end

    let tmp_str: String val = if iexp >= 0 then
      let ie = iexp.usize()
      let exp_val_str = s.substring(iexp + 1)
      if exp_val_str.size() == 0 then
        error
      end
      
      var exp: ILong = 0
      var seen_digit = false
      var j: USize = 0
      if try exp_val_str(0)? == '+' else false end then
        j = 1
      end
      if j >= exp_val_str.size() then
        error
      end

      while j < exp_val_str.size() do
        let ec = exp_val_str(j)?
        if (ec >= '0') and (ec <= '9') then
          exp = (exp * 10) + (ec - '0').ilong()
          seen_digit = true
        elseif ec == '_' then
          None
        elseif ec == ' ' then
          error // Trailing space in exponent
        else
          error
        end
        j = j + 1
      end

      if not seen_digit then
        error
      end

      recover val
        let pads = String.create(exp.usize())
        var k: ILong = 0
        while k < exp do
          pads.push('0')
          k = k + 1
        end
        s.substring(i.isize(), iexp) + pads
      end
    else
      s.substring(i.isize())
    end

    if tmp_str.size() == 0 then
      error
    end

    _digits = recover val
      let d: Array[U16] ref = [0]
      var has_real_digit = false
      for dc in tmp_str.values() do
        if dc == '_' then
          continue
        end
        if dc == ' ' then
          error
        end // No spaces inside number
        let digit: U32 = if (dc >= '0') and (dc <= '9') then
          (dc - '0').u32()
        elseif (dc >= 'a') and (dc <= 'z') then
          (dc - 'a').u32() + 10
        elseif (dc >= 'A') and (dc <= 'Z') then
          (dc - 'A').u32() + 10
        else
          error
        end
        if digit >= base then
          error
        end
        _short_mul(d, base)
        _short_add(d, digit)
        has_real_digit = true
      end
      if not has_real_digit then
        error
      end
      _normalize(d)
      let d2 = Array[U16](d.size())
      for x in d.values() do
        d2.push(x)
      end
      d2
    end
    _negative = is_neg


  new val from_ilong(n: ILong = 0) =>
    """
    Create a new arbitrary-precision integer with value `n` (default 0).
    """
    _negative = (n < 0)
    var q: ULong = n.abs().ulong()
    let base: ULong = _base().ulong()
    _digits = recover val
      let d: Array[U16] iso = Array[U16]
      while q >= base do
        (q, let r) = q.divrem(base)
        d.push(r.u16())
      end
      d.push(q.u16())
      consume d
    end


  new val _from_array(negative: Bool, digits: Array[U16] val) =>
    """
    Create a new `MPInt` from its sign `negative` and its internal array
    representation `digits`.
    """
    _negative = negative
    _digits = digits


  new val create(value: MPInt) =>
    """
    The default copy constructor, sharing internally the same numbers.
    Required by Real trait.
    """
    _negative = value._negative
    _digits = value._digits


  new val from[A: (Number & Real[A] val)](a: A) =>
    """
    Create an `MPInt` from another number.

    CAUTION: This constructor converts `a` to a `U128` internally and
    overflow can occur, for instance when constructing the `MPInt` from a
    large `F64`.
    """
    _negative = (a.f64() < 0)
    let base = U128.from[U32](_base())
    _digits = recover val
      var q: U128 = if _negative then
        (-a).u128()
      else
        a.u128()
      end
      let d: Array[U16] iso = Array[U16]
      while q >= base do
        (q, let r) = q.divrem(base)
        d.push(r.u16())
      end
      d.push(q.u16())
      consume d
    end


  new \do_not_use\ val min_value() =>
    """
    There's no minimal `MPInt` value, so we return the `0` value.
    Required by Real trait.

    DO NOT USE THIS METHOD.
    """
    _negative = false
    _digits = [0]


  new \do_not_use\ val max_value() =>
    """
    There's no maximal `MPInt` value, so we return the `0` value.
    Required by Real trait.

    DO NOT USE THIS METHOD.
    """
    _negative = false
    _digits = [0]


  //- Status ------------------------------------------------------------------

  fun bitwidth(): MPInt =>
    """
    Returns the number of bits used in the `MPInt` representation.
    Contrarily to a fixed-width integer like `I32`, this value is variable
    and depends on the content of the `MPInt`. 

    As an `MPInt` object is not represented by a array of bits, like
    traditional binary integer representation, the size returned by this
    function is an approximate of the real memory size. The sign of the
    `MPInt` is assumed to be 1 bit, and each *digit* is supposed to be
    `base_bit` bits.

    CAUTION, for **very** large `MPInt`, this method could overflow without
    raising an error as calculations are done with `USize` values.
    """
    if is_zero() then
      return MPInt.from_ilong(0)
    end
    let sign: USize = 1
    let n = _digits.size() - 1
    let last = try _digits(n)? else 0 end
    let bits = _base_bits() - last.clz().usize()
    MPInt.from_ilong(((n * _base_bits()) + bits + sign).ilong())


  fun bytewidth(): USize =>
    """
    Returns the number of bytes used in the `MPInt` representation.
    Contrarily to a fixed-width integer like `I32`, this value is variable
    and depends on the content of the `MPInt`. 

    As an `MPInt` object is not represented by a array of bits, like
    traditional binary integer representation, the size returned by this
    function is an approximate of the real memory size. The sign of the
    `MPInt` is assumed to be 1 bit, and each *digit* is supposed to be
    `base_bit` bits.

    CAUTION, for **very** large `MPInt`, this method could overflow without
    raising an error as calculations are done with `USize` values.
    """
    ((bitwidth().ilong().usize() + 7) / 8)


  fun popcount(): MPInt =>
    """
    Returns the number of bits set to 1 in the absolute value of the `MPInt`.
    """
    var c: USize = 0
    for d in _digits.values() do
      c = c + d.popcount().usize()
    end
    MPInt.from_ilong(c.ilong())


  fun \do_not_use\ clz(): MPInt =>
    """
    Leading zeros count. For arbitrary precision, returns 0.

    DO NOT USE THIS METHOD.
    """
    MPInt.from_ilong(0)


  fun ctz(): MPInt =>
    """
    Trailing zeros count in the absolute value.
    """
    if is_zero() then
      return MPInt.from_ilong(0)
    end
    var c: USize = 0
    for d in _digits.values() do
      if d == 0 then
        c = c + _base_bits()
      else
        c = c + d.ctz().usize()
        break
      end
    end
    MPInt.from_ilong(c.ilong())


  fun \do_not_use\ clz_unsafe(): MPInt =>
    """
    Unsafe leading zeros count. Always returns 0.

    DO NOT USE THIS METHOD.
    """
    clz()


  fun ctz_unsafe(): MPInt =>
    """
    Unsafe trailing zeros count.
    """
    ctz()


  fun \do_not_use\ rotl(that: MPInt val): MPInt =>
    """
    Rotate left. Not supported for arbitrary precision, returns a copy.

    DO NOT USE THIS METHOD. NOT IMPLEMENTED.
    """
    _create(_negative, _digits)


  fun \do_not_use\ rotr(that: MPInt val): MPInt =>
    """
    Rotate right. Not supported for arbitrary precision, returns a copy.

    DO NOT USE THIS METHOD. NOT IMPLEMENTED.
    """
    _create(_negative, _digits)


  fun bit_shl(that: MPInt val): MPInt =>
    """
    Bit-level shift left: return `this * 2^that`, preserving sign.
    If `that` is negative or `this` is zero, return a copy unchanged.

    The shift is decomposed into a *digit*-level shift (`that / 16`) and a
    sub-digit bit shift (`that % 16`), with carry propagation between *digits*.
    """
    if that._negative or is_zero() then
      return _create(_negative, _digits)
    end

    let n = that.usize()
    let digit_shift = n / _base_bits()
    let bit_shift = n % _base_bits()
    let d = recover val
      let res = Array[U16](_digits.size() + digit_shift + 1)
      // Zero-fill the low digits
      for _ in Range(0, digit_shift) do
        res.push(0)
      end

      // Shift each digit left, propagating carry
      var carry: U32 = 0
      for digit in _digits.values() do
        let v = (digit.u32() << bit_shift.u32()) or carry
        res.push(v.u16())
        carry = v >> 16
      end
      if carry > 0 then
        res.push(carry.u16())
      end
      res
    end
    _create(_negative, d)


  fun bit_shr(that: MPInt val): MPInt =>
    """
    Bit-level shift right: return `this / 2^that` (truncating), preserving sign.
    If `that` is negative or `this` is zero, return a copy unchanged.
    If `that >= bitwidth()`, return 0.

    The shift is decomposed into a digit-level shift (`that / 16`) and a
    sub-digit bit shift (`that % 16`), reading overlapping bits from adjacent digits.
    """
    if that._negative or is_zero() then
      return _create(_negative, _digits)
    end

    let n = that.usize()
    let digit_shift = n / _base_bits()
    let bit_shift = n % _base_bits()
    if digit_shift >= _digits.size() then
      return MPInt.from_ilong(0)
    end

    let new_size = _digits.size() - digit_shift
    let d = recover val
      let res = Array[U16](new_size)
      for i in Range(0, new_size) do
        let lo: U16 = try _digits(i + digit_shift)? >> bit_shift.u16() else 0 end
        let hi: U16 = if (bit_shift > 0) and ((i + digit_shift + 1) < _digits.size()) then
          try (_digits(i + digit_shift + 1)? << (16 - bit_shift).u16()) else 0 end
        else
          0
        end
        res.push(lo or hi)
      end
      _normalize(res)
      res
    end
    _create(_negative, d)


  fun shl(that: MPInt val): MPInt =>
    """
    Bitwise shift left. Delegates to `bit_shl`.

    See [`bit_shl`](#bit_shl)
    """
    bit_shl(that)


  fun shr(that: MPInt val): MPInt =>
    """
    Bitwise shift right. Delegates to `bit_shr`.

    See [`bit_shr`](#bit_shr)
    """
    bit_shr(that)


  fun shl_unsafe(that: MPInt val): MPInt =>
    """
    Bitwise shift left (no overflow check; it can't happen).
    Delegates to `bit_shl`.

    See [`bit_shl`](#bit_shl)
    """
    bit_shl(that)


  fun shr_unsafe(that: MPInt val): MPInt =>
    """
    Bitwise shift right (no overflow check; it can't happen).
    Delegates to `bit_shr`.

    See [`bit_shr`](#bit_shr)
    """
    bit_shr(that)


  //- Arithmetic --------------------------------------------------------------

  fun neg(): MPInt =>
    """
    Return a new `MPInt` with value `-this`.
    """
    _create(not _negative, _digits)


  fun add(that: MPInt): MPInt =>
    """
    Return a new `MPInt` with the value `this + that`.
    """
    if _negative == that._negative then
      let d = recover val
        let res = Array[U16](_digits.size().max(that._digits.size()) + 1)
        for x in _digits.values() do
          res.push(x)
        end
        _add_arrays(res, that._digits)
        res
      end
      _create(_negative, d)
    else
      if not this._uabs_lt(that) then
        let d = recover val
          let res = Array[U16](_digits.size())
          for x in _digits.values() do
            res.push(x)
          end
          _sub_arrays(res, that._digits)
          res
        end
        _create(_negative, d)
      else
        let d = recover val
          let res = Array[U16](that._digits.size())
          for x in that._digits.values() do
            res.push(x)
          end
          _sub_arrays(res, _digits)
          res
        end
        _create(that._negative, d)
      end
    end


  fun sub(that: MPInt): MPInt =>
    """
    Return a new `MPInt` whose value is equal to the subtraction `this - that`.
    """
    if _negative != that._negative then
      let d = recover val
        let res = Array[U16](_digits.size().max(that._digits.size()) + 1)
        for x in _digits.values() do
          res.push(x)
        end
        _add_arrays(res, that._digits)
        res
      end
      _create(_negative, d)
    else
      if not this._uabs_lt(that) then
        let d = recover val
          let res = Array[U16](_digits.size())
          for x in _digits.values() do
            res.push(x)
          end
          _sub_arrays(res, that._digits)
          res
        end
        _create(_negative, d)
      else
        let d = recover val
          let res = Array[U16](that._digits.size())
          for x in that._digits.values() do
            res.push(x)
          end
          _sub_arrays(res, _digits)
          res
        end
        _create(not _negative, d)
      end
    end


  fun mul(that: MPInt): MPInt =>
    """
    Multiply `this` by `that` and return a new `MPInt`.
    """
    if is_zero() or that.is_zero() then
      return MPInt.from_ilong(0)
    end
    if is_one() then
      return _create(that._negative, that._digits)
    end
    if is_minus_one() then
      return _create(not that._negative, that._digits)
    end
    if that.is_one() then
      return _create(_negative, _digits)
    end
    if that.is_minus_one() then
      return neg()
    end

    let d = recover val
      let res = Array[U16].init(0, _digits.size() + that._digits.size())
      var i: USize = 0
      while i < _digits.size() do
        var carry: U32 = 0
        var j: USize = 0
        try
          while j < that._digits.size() do
            let m = (_digits(i)?.u32() * that._digits(j)?.u32()) +
                    res(i + j)?.u32() + carry
            carry = m >> _base_bits().u32()
            res.update(i + j, m.u16())?
            j = j + 1
          end
          res.update(i + that._digits.size(), carry.u16())?
        end
        i = i + 1
      end
      _normalize(res)
      let d2 = Array[U16](res.size())
      for x in res.values() do
        d2.push(x)
      end
      d2
    end
    _create(_negative xor that._negative, d)


  fun fast_mul(that: MPInt): MPInt =>
    """
    Fast multiplication of `this` and `that` using FFT. The result has a size
    of `this.size() + that.size()`.

    See Knuth TAOCP 4.3.3.C
    Why is it fast? See https://en.wikipedia.org/wiki/Sch%C3%B6nhage%E2%80%93Strassen_algorithm
    """
    if is_zero() or that.is_zero() then
      return MPInt.from_ilong(0)
    end
    if is_one() then
      return _create(that._negative, that._digits)
    end
    if is_minus_one() then
      return _create(not that._negative, that._digits)
    end
    if that.is_one() then
      return _create(_negative, _digits)
    end
    if that.is_minus_one() then
      return neg()
    end

    let pow2 = (_digits.size() + that._digits.size()).next_pow2()
    let d = recover val
      let res = Array[U16].init(0, _digits.size() + that._digits.size())
      try
        var a = Array[F64].init(0, pow2)
        var i: USize = 0
        for x in _digits.values() do
          a.update(i, x.f64())?
          i = i + 1
        end
        
        var b = Array[F64].init(0, pow2)
        i = 0
        for x in that._digits.values() do
          b.update(i, x.f64())?
          i = i + 1
        end

        FFT.fourier_real(a)
        FFT.fourier_real(b)
        
        // Complex multiplication in packed format:
        // F(0) and F(N/2) are real and stored in b(0) and b(1)
        b.update(0, b(0)? * a(0)?)?
        b.update(1, b(1)? * a(1)?)?
        i = 2
        while i < pow2 do
          let ar = a(i)?
          let ai = a(i+1)?
          let br = b(i)?
          let bi = b(i+1)?
          b.update(i, (ar * br) - (ai * bi))?
          b.update(i+1, (ar * bi) + (ai * br))?
          i = i + 2
        end
        FFT.fourier_real(b, true)

        var carry: F64 = 0
        var k: USize = 0
        let base: F64 = _base().f64()
        while k < res.size() do
          let temp = b(k)? + carry + 0.5
          let carry' = (temp / base).trunc()
          res.update(k, (temp - (carry' * base)).u16())?
          carry = carry'
          k = k + 1
        end
      end
      _normalize(res)
      let d2 = Array[U16](res.size())
      for x in res.values() do
        d2.push(x)
      end
      d2
    end
    _create(_negative xor that._negative, d)


  fun karatsuba_mul(that: MPInt): MPInt =>
    """
    Multiply `this` by `that`using the
    [Karatsuba algorithm](https://en.wikipedia.org/wiki/Karatsuba_algorithm).
    This algorith uses a recursive divide and conquer approach that requires
    creating new temporary `MPInt`. It reduces the multiplication of two n-
    digits numbers to 3 multiplications of n/2-digits, and by recursion
    providing at the limit O(n^log2(3)) operations.

    The Karatsuba multiplication triggers only when the numbers are more
    than 2048-bits long (128 base-digits). When the size is less than that
    limit, classical `mul` multiplication is used.
    """
    if (_digits.size() <= 128) or (that._digits.size() <= 128) then
      mul(that)
    else
      let size = _digits.size().max(that._digits.size())
      let half = size / 2
      let this_low = _create(_negative, recover val 
        let d_low = Array[U16](half)
        for k in Range(0, half) do
          try d_low.push(_digits(k)?) end
        end
        d_low
      end)
      let this_high = _create(_negative, recover val
        let d_high = Array[U16](_digits.size() - half)
        for k in Range(half, _digits.size()) do
          try d_high.push(_digits(k)?) end
        end
        d_high
      end)
      let that_low = _create(that._negative, recover val
        let d_low_that = Array[U16](half)
        for k in Range(0, half) do
          try d_low_that.push(that._digits(k)?) end
        end
        d_low_that
      end)
      let that_high = _create(that._negative, recover val
        let d_high_that = Array[U16](that._digits.size() - half)
        for k in Range(half, that._digits.size()) do
          try d_high_that.push(that._digits(k)?) end
        end
        d_high_that
      end)

      let z2 = this_high.karatsuba_mul(that_high)
      let z0 = this_low.karatsuba_mul(that_low)
      let z1 = (this_low + this_high).karatsuba_mul(that_low + that_high) - z2 - z0
      z2.digit_shl(2 * half) + z1.digit_shl(half) + z0
    end


  //- Division ----------------------------------------------------------------

  fun divrem(that: MPInt): (MPInt, MPInt) =>
    """
    Calculates the quotient and the remainder where `this = quotient * that + remainder`.
    The result is `(this / that, this % that)`.
    """
    // Division by 0
    if that.is_zero() then
      return (MPInt.from_ilong(0), MPInt.from_ilong(0))
    end
    if _uabs_lt(that) then
      return (MPInt.from_ilong(0), _create(_negative, _digits))
    end
    
    // Division by a single digit
    if that._digits.size() == 1 then
      let v0 = try that._digits(0)? else 1 end
      let q_digits = recover val
        let d_q = Array[U16](_digits.size())
        for x in _digits.values() do
          d_q.push(x)
        end
        _short_div(d_q, v0.u32())
        let d2_q = Array[U16](d_q.size())
        for x in d_q.values() do
          d2_q.push(x)
        end
        d2_q
      end
      let r_val = recover val
        let d_r = Array[U16](_digits.size())
        for x in _digits.values() do
          d_r.push(x)
        end
        let rem_val = _short_div(d_r, v0.u32())
        [rem_val.u16()]
      end
      return (_create(_negative xor that._negative, q_digits), _create(_negative, r_val))
    end
    
    _divrem_alg_d(that)


  fun _divrem_alg_d(that: MPInt): (MPInt, MPInt) =>
    """
    Implementation of Knuth's Algorithm D for multi-precision division.
    TAOCP Vol 2, 4.3.1.
    """
    // Constants
    let base2: U32 = _base() / 2
    let base: U64 = _base().u64()

    // Step D1: Normalize
    let v_high = try that._digits(that._digits.size() - 1)? else 0 end
    var s_shift: U32 = 0
    var v_h = v_high.u32()
    while v_h < base2 do
      v_h = v_h << 1
      s_shift = s_shift + 1
    end
    let d_val: U32 = U32(1) << s_shift

    let u_ref: Array[U16] iso = recover iso
      let a: Array[U16] = Array[U16].create(_digits.size() + 1)
      for x in _digits.values() do a.push(x) end
      _short_mul(a, d_val)
      if a.size() == _digits.size() then
        a.push(0)
      end
      a
    end

    let v_ref: Array[U16] val = recover val
      let a: Array[U16] = Array[U16].create(that._digits.size())
      for x in that._digits.values() do a.push(x) end
      _short_mul(a, d_val)
      a
    end

    let n_val: USize = v_ref.size()
    let m_val: USize = u_ref.size() - n_val - 1
    let q_ref: Array[U16] iso = recover iso Array[U16].init(0, m_val + 1) end

    let v_n_1: U64 = try v_ref(n_val - 1)?.u64() else 0 end
    let v_n_2: U64 = try v_ref(n_val - 2)?.u64() else 0 end

    let shift = _base_bits().u64()
    let mask = _base().i64() - 1

    // Step D2: Loop on j
    var j_idx: USize = m_val
    try
      while true do
        // Step D3: Calculate q_hat
        let u_j_n: U64 = u_ref(j_idx + n_val)?.u64()
        let u_j_n_1: U64 = u_ref((j_idx + n_val) - 1)?.u64()
        let u_j_n_2: U64 = u_ref((j_idx + n_val) - 2)?.u64()

        let u_2_digits = (u_j_n << shift) + u_j_n_1
        var q_hat: U64 = u_2_digits / v_n_1
        var r_hat: U64 = u_2_digits % v_n_1

        while true do
          if (q_hat >= base) or ((q_hat * v_n_2) > ((r_hat << shift) + u_j_n_2)) then
            q_hat = q_hat - 1
            r_hat = r_hat + v_n_1
            if r_hat < base then
              continue
            end
          end
          break
        end

        // Step D4: Multiply and subtract
        var borrow: I64 = 0
        for i_sub in Range(0, n_val) do
          let p = q_hat * v_ref(i_sub)?.u64()
          let sub_res: I64 = u_ref(j_idx + i_sub)?.i64() - (p.i64() and mask) - borrow
          u_ref.update(j_idx + i_sub, (sub_res and mask).u16())?
          borrow = (p.i64() >> shift) - (sub_res >> shift)
        end
        
        let sub_res_last: I64 = u_ref(j_idx + n_val)?.i64() - borrow
        u_ref.update(j_idx + n_val, (sub_res_last and mask).u16())?
        
        if sub_res_last < 0 then
          // Step D5: Negative case - add back
          q_hat = q_hat - 1
          var carry: U32 = 0
          for i_add in Range(0, n_val) do
            let sum: U32 = u_ref(j_idx + i_add)?.u32() + v_ref(i_add)?.u32() + carry
            u_ref.update(j_idx + i_add, sum.u16())?
            carry = sum >> _base_bits().u32()
          end
          u_ref.update(j_idx + n_val, (u_ref(j_idx + n_val)?.u32() + carry).u16())?
        end

        q_ref.update(j_idx, q_hat.u16())?

        if j_idx == 0 then break end
        j_idx = j_idx - 1
      end
    end

    let q_digits: Array[U16] val = recover val
      let a: Array[U16] ref = consume q_ref
      _normalize(a)
      a
    end

    // Step D7: Unnormalize
    let r_digits: Array[U16] val = recover val
      let a: Array[U16] ref = Array[U16].init(0, n_val)
      try
        var k_rem: USize = 0
        while k_rem < n_val do
          a.update(k_rem, u_ref(k_rem)?)?
          k_rem = k_rem + 1
        end
      end
      _short_div(a, d_val)
      _normalize(a)
      a
    end

    (_create(_negative xor that._negative, q_digits), _create(_negative, r_digits))


  fun div(that: MPInt): MPInt =>
    """
    Integer division. Result is truncated quotiont: `this / that`.

    See [`rem`](#rem)
    """
    divrem(that)._1


  fun rem(that: MPInt): MPInt =>
    """
    Remainder of integer division. Result is truncated remainder of division
    `this / that`, that is `this % that`

    See [`div`](#div)
    """
    divrem(that)._2


  fun mod(that: MPInt): MPInt =>
    """
    Modulo operation gining the floored remainder of division of `this` by `that`:
    `this.mod(that)` follows the sign of the divisor (`that`) and is defined such
    that `this = fld(quotient) * that + mod(remainder)`.

    See [`MPInt.fld`](#fld)
    """
    let r = rem(that)
    if (not r.is_zero()) and (_negative xor that._negative) then
      r + that
    else
      r
    end


  fun fld(y: MPInt): MPInt =>
    """
    Floored integer division of `this` by `that`. The equation `this = fld(quotient) * that + mod(remainder)`
    holds where `this = quotient * that + remainder`.
    """
    (let q, let r) = divrem(y)
    if (not r.is_zero()) and (_negative xor y._negative) then
      q - MPInt.from_ilong(1)
    else
      q
    end


  fun pow(n: MPInt val): MPInt =>
    """
    Return `this` raised to the power `n` using binary (square-and-multiply)
    exponentiation. `O(log n)` multiplications.

    - `x^0 = 1` for any `x` (including 0).
    - `x^n` for negative `n` returns 0 (integer semantics, no fractions).
    """
    let one = MPInt.from_ilong(1)
    if n.is_zero() then
      return one
    end
    if n._negative then
      return MPInt.from_ilong(0)
    end

    var base: MPInt = _create(_negative, _digits)
    var exp: MPInt  = n
    var result: MPInt = one
    while not exp.is_zero() do
      if exp.bit_get(0) then
        result = result * base
      end
      base = base * base
      exp  = exp.bit_shr(one)
    end
    result


  fun isqrt(): MPInt =>
    """
    Return the integer square root, i.e. `floor(sqrt(this))`.
    Uses Newton's method, which converges quadratically without any
    floating-point arithmetic.

    Returns 0 for negative or zero input.
    """
    if _negative or is_zero() then
      return MPInt.from_ilong(0)
    end
    let one = MPInt.from_ilong(1)
    if is_one() then
      return one
    end

    // Initial guess: 2^(ceil(bitwidth/2)), always >= sqrt(this).
    let half_bits = ((bitwidth().usize() + 1) / 2)
    var x: MPInt  = one.bit_shl(MPInt.from_ilong(half_bits.ilong()))
    var x1: MPInt = (x + (this / x)).bit_shr(one)
    while x1 < x do
      x  = x1
      x1 = (x + (this / x)).bit_shr(one)
    end
    x


  fun gcd(that: MPInt): MPInt =>
    """
    Return the greatest common divisor of `|this|` and `|that|` using the
    Euclidean algorithm.

    `gcd(0, x) = x`, `gcd(x, 0) = x`, `gcd(0, 0) = 0`.
    """
    var a: MPInt = abs()
    var b: MPInt = that.abs()
    while not b.is_zero() do
      let t = b
      b = a.rem(b)
      a = t
    end
    a


  fun pow_mod(exp: MPInt val, m: MPInt val): MPInt =>
    """
    Return `(this ^ exp) mod m` using square-and-multiply with modular
    reduction at each step. Much faster than `pow(exp).mod(m)` for large
    exponents since intermediate values stay bounded by `m^2`.

    - `pow_mod(0, m) = 1` (for m > 1).
    - Returns 0 for negative `exp` (integer semantics).
    - Returns 0 for `m = 1` since every integer is congruent to 0 mod 1.
    """
    let zero = MPInt.from_ilong(0)
    let one = MPInt.from_ilong(1)
    if m.is_one() then
      return zero
    end
    if exp.is_zero() then
      return one
    end
    if exp._negative then
      return zero
    end

    var base: MPInt   = abs().mod(m)
    var e: MPInt      = exp
    var result: MPInt = one
    while not e.is_zero() do
      if e.bit_get(0) then
        result = (result * base).mod(m)
      end
      base = (base * base).mod(m)
      e = e.bit_shr(one)
    end
    result


  fun add_partial(that: MPInt val): MPInt ? =>
    """
    Partial addition. `MPInt` can't overflow or underflow and no errors are raised.

    See [`add`](#add)
    """
    if false then error end
    add(that)


  fun sub_partial(that: MPInt val): MPInt ? =>
    """
    Partial subtraction. `MPInt` can't overflow or underflow and no errors are raised.

    See [`sub`](#sub)
    """
    if false then error end
    sub(that)


  fun mul_partial(that: MPInt val): MPInt ? =>
    """
    Partial multiplication. `MPInt` can't overflow or underflow and no errors are raised.

    See [`mul`](#mul)
    """
    if false then error end
    mul(that)


  fun div_partial(that: MPInt): MPInt ? =>
    """
    Partial division. It raises an error only when dividing by zero (`that == 0`).

    See [`div`](#div)
    """
    if that.is_zero() then
      error
    else
      div(that)
    end


  fun rem_partial(that: MPInt): MPInt ? =>
    """
    Partial remainder. It raises an error only when dividing by zero (`that == 0`).

    see [`rem`](#rem)
    """
    if that.is_zero() then
      error
    else
      rem(that)
    end


  fun divrem_partial(that: MPInt): (MPInt, MPInt) ? =>
    """
    Partial divrem. It raises an error only when dividing by zero (`that == 0`).

    See [`divrem`][#divrem]
    """
    if that.is_zero() then
      error
    else
      divrem(that)
    end


  fun fld_partial(that: MPInt): MPInt ? =>
    """
    Partial floored division. It raises an error only when dividing by zero (`that == 0`).

    See [`fld`](#fld)
    """
    if that.is_zero() then
      error
    else
      fld(that)
    end


  fun mod_partial(that: MPInt): MPInt ? =>
    """
    Partial modulo. It raises an error only when dividing by zero (`that == 0`).

    See [`mod`](#mod)
    """
    if that.is_zero() then
      error
    else
      mod(that)
    end


  fun add_unsafe(that: MPInt): MPInt =>
    """
    Unsafe addition.

    See [`add`](#add)
    """
    add(that)


  fun sub_unsafe(that: MPInt): MPInt =>
    """
    Unsafe subtraction.

    See [`sub`](#sub)
    """
    sub(that)


  fun mul_unsafe(that: MPInt): MPInt =>
    """
    Unsafe multiplication.

    See [`mul`](#mul)
    """
    mul(that)


  fun div_unsafe(that: MPInt): MPInt =>
    """
    Unsafe division.

    See [`div`](#div)
    """
    div(that)


  fun rem_unsafe(that: MPInt): MPInt =>
    """
    Unsafe remainder.

    See [`rem`](#rem)
    """
    rem(that)


  fun divrem_unsafe(that: MPInt): (MPInt, MPInt) =>
    """
    Unsafe divrem.

    See [`divrem`](#divrem)
    """
    divrem(that)


  fun fld_unsafe(that: MPInt): MPInt =>
    """
    Unsafe floored division.

    See [`fld`](#fld)
    """
    fld(that)


  fun mod_unsafe(that: MPInt): MPInt =>
    """
    Unsafe modulo.

    See [`mod`](#mod)
    """
    mod(that)


  fun addc(that: MPInt): (MPInt, Bool) =>
    """
    Addition with carry. No overflow/underflow are possible and the boolean
    is always `false`.

    See [`add`](#add)
    """
    (add(that), false)


  fun subc(that: MPInt): (MPInt, Bool) =>
    """
    Subtraction with borrow. No overflow/underflow are possible and the boolean
    is always `false`.

    See [`sub`](#sub)
    """
    (sub(that), false)


  fun mulc(that: MPInt): (MPInt, Bool) =>
    """
    Multiplication with carry. No overflow/underflow are possible and the boolean
    is always `false`.

    See [`mul`](#mul)
    """
    (mul(that), false)


  fun divc(that: MPInt): (MPInt, Bool) =>
    """
    Division with overflow check. The boolean is set when `that` is 0.

    See [`div`](#div)
    """
    (div(that), that.is_zero())


  fun remc(that: MPInt): (MPInt, Bool) =>
    """
    Remainder with overflow check. The boolean is set when `that` is 0.

    See [`rem`](#rem)
    """
    (rem(that), that.is_zero())


  fun fldc(that: MPInt): (MPInt, Bool) =>
    """
    Floored division with overflow check. The boolean is set when `that` is 0.

    See [`fld`](#fld)
    """
    (fld(that), that.is_zero())


  fun modc(that: MPInt): (MPInt, Bool) =>
    """
    Modulo with overflow check. The boolean is set when `that` is 0.

    See [`mod`](#mod)
    """
    (mod(that), that.is_zero())


  //- Bitwise -----------------------------------------------------------------

  fun op_and(that: MPInt): MPInt =>
    """
    Bitwise AND using two's complement semantics for arbitrary precision.
    """
    _bitwise_op(that, {(a: BitMap ref, b: BitMap ref): BitMap ref => a.and_in_place(b)})


  fun op_or(that: MPInt): MPInt =>
    """
    Bitwise OR using two's complement semantics for arbitrary precision.
    """
    _bitwise_op(that, {(a: BitMap ref, b: BitMap ref): BitMap ref => a.or_in_place(b)})


  fun op_xor(that: MPInt): MPInt =>
    """
    Bitwise XOR using two's complement semantics for arbitrary precision.
    """
    _bitwise_op(that, {(a: BitMap ref, b: BitMap ref): BitMap ref => a.xor_in_place(b)})


  fun _bitwise_op(that: MPInt,
      op: {(BitMap ref, BitMap ref): BitMap ref} ref): MPInt =>
    """
    Convert both operands to two's complement `BitMap`s, apply `op` in place on
    the first, then convert the result back to sign-magnitude `MPInt`.
    """
    // Size the BitMaps to hold both operands plus one sign bit.
    let n_bits = ((_digits.size().max(that._digits.size()) + 1) * _base_bits()) + 1
    let a = BitMap.from_u16_array(_digits, n_bits)
    if _negative then
      a.not_in_place()
      a.increment()
    end
    let b = BitMap.from_u16_array(that._digits, n_bits)
    if that._negative then
      b.not_in_place()
      b.increment()
    end

    op(a, b)

    // Two's complement sign bit determines result sign. And convert back
    // to MPInt.
    let negative = a(n_bits - 1)
    if negative then
      a.not_in_place()
      a.increment()
    end
    let n_digits = (n_bits + (_base_bits() - 1)) / _base_bits()
    let raw = a.to_u16_array(n_digits)
    let d = recover val
      let res = Array[U16](raw.size())
      for v in raw.values() do res.push(v) end
      _normalize(res)
      res
    end
    _create(negative, d)


  fun op_not(): MPInt =>
    """
    Bitwise NOT using two's complement identity: `~x = -(x + 1)`.
    """
    neg() - MPInt.from_ilong(1)


  fun bit_get(n: USize): Bool =>
    """
    Return the value of bit `n` in the absolute value of this `MPInt`.
    Bit 0 is the least significant. Returns `false` for bits above the
    most significant digit.
    """
    let digit_idx = n / _base_bits()
    let bit_pos = (n % _base_bits()).u16()
    if digit_idx >= _digits.size() then
      false
    else
      try ((_digits(digit_idx)? >> bit_pos) and 1) == 1 else false end
    end


  fun bit_set(n: USize): MPInt =>
    """
    Return a new `MPInt` with bit `n` of the absolute value set to 1.
    The sign is preserved.
    """
    let digit_idx = n / _base_bits()
    let bit_pos = (n % _base_bits()).u16()
    let d = recover val
      let res = Array[U16](_digits.size().max(digit_idx + 1))
      for k in Range(0, _digits.size()) do
        try res.push(_digits(k)?) end
      end
      while res.size() <= digit_idx do res.push(0) end
      try res(digit_idx)? = res(digit_idx)? or (U16(1) << bit_pos) end
      res
    end
    _create(_negative, d)


  fun bit_clear(n: USize): MPInt =>
    """
    Return a new `MPInt` with bit `n` of the absolute value cleared to 0.
    The sign is preserved. If the bit is already 0, returns a copy.
    """
    if not bit_get(n) then
      return _create(_negative, _digits)
    end
    let digit_idx = n / _base_bits()
    let bit_pos = (n % _base_bits()).u16()
    let d = recover val
      let res = Array[U16](_digits.size())
      for k in Range(0, _digits.size()) do
        try res.push(_digits(k)?) end
      end
      try res(digit_idx)? = res(digit_idx)? and (not (U16(1) << bit_pos)) end
      _normalize(res)
      res
    end
    _create(_negative, d)


  fun bit_flip(n: USize): MPInt =>
    """
    Return a new `MPInt` with bit `n` of the absolute value toggled.
    The sign is preserved.
    """
    if bit_get(n) then bit_clear(n) else bit_set(n) end


  fun \do_not_use\ bit_reverse(): MPInt =>
    """
    Reverse bits.

    DO NOT USE THIS METHOD. NOT IMPLEMENTED.
    """
    _create(_negative, _digits) // TODO


  fun \do_not_use\ bswap(): MPInt =>
    """
    Swap bytes.

    DO NOT USE THIS METHOD. NOT IMPLEMENTED.
    """
    _create(_negative, _digits) // TODO


  //- Conversions -------------------------------------------------------------
  // All conversions are silent: no error or panic is raised on overflow or
  // precision loss. Use `bitwidth()` to verify the value fits before converting.

  fun i8(): I8 =>
    """
    Convert to `I8`. Delegates to `ilong().i8()`, so two silent truncations
    can occur: `ilong()` clips to the low 64 bits of `|this|`, then `.i8()`
    further clips to 8 bits. Safe only when `-128 ≤ this ≤ 127`.
    """
    ilong().i8()


  fun i16(): I16 =>
    """
    Convert to `I16`. Delegates to `ilong().i16()`, so two silent truncations
    can occur: `ilong()` clips to the low 64 bits of `|this|`, then `.i16()`
    further clips to 16 bits. Safe only when `-32768 ≤ this ≤ 32767`.
    """
    ilong().i16()


  fun i32(): I32 =>
    """
    Convert to `I32`. Delegates to `ilong().i32()`, so two silent truncations
    can occur: `ilong()` clips to the low 64 bits of `|this|`, then `.i32()`
    further clips to 32 bits. Safe only when `-2^31 ≤ this ≤ 2^31 - 1`.
    """
    ilong().i32()


  fun i64(): I64 =>
    """
    Convert to `I64`. Delegates to `ilong().i64()`. On 64-bit platforms
    `ILong` is `I64`, so `.i64()` is a no-op cast and the only data loss is
    the `ilong()` truncation: values outside `[-2^63, 2^63 - 1]` silently
    return the low 64 bits of `|this|`, negated if negative.
    """
    ilong().i64()


  fun i128(): I128 =>
    """
    Convert to `I128` by reconstructing the value from the little-endian
    base-65536 *digit* array using I128 accumulation.

    **Overflow behaviour (silent, no error).**
    Safe for `|this| ≤ I128.max_value()` (= 2^127 - 1). When the `MPInt`
    value exceeds that range, the `weight` variable wraps to 0 at *digit* index
    8 (`65536^8 = 2^128 ≡ 0`), higher *digits* contribute nothing, and the
    intermediate additions wrap silently. The net effect is that the function
    returns the low 128 bits of `|this|` reinterpreted as `I128`, then negated
    if negative.
    """
    var res: I128 = 0
    var weight: I128 = 1
    for d in _digits.values() do
      res = res +~ (d.i128() *~ weight)
      weight = weight *~ _base().i128()
    end
    if _negative then
      -res
    else
      res
    end


  fun ilong(): ILong =>
    """
    Convert to `ILong` by reconstructing the value from the little-endian
    base-65536 digit array: `sum(digit[i] * 65536^i)`, negated when negative.

    **Overflow behaviour (silent, no error).**
    `ILong` is a signed 64-bit integer on 64-bit platforms (LP64 / C `long`).
    When `|this| > ILong.max_value()` (= 2^63 - 1 on 64-bit), the
    intermediate arithmetic wraps in the following way:

    - The per-digit `weight` value (`65536^i`) overflows to 0 starting at
      `i = 4` (`65536^4 = 2^64 ≡ 0` under 64-bit wrapping).
    - All digits at position 4 and beyond therefore contribute 0 to `res`.
    - The additions and multiplications on `res` itself also wrap silently.

    The net effect is that the function returns the value of the **lowest 64
    bits** of `|this|` reinterpreted as `ILong`, then negated if negative —
    i.e. `this mod 2^64` truncated to `ILong`, without any indication that
    data was lost. Check [`bitwidth`](#bitwidth) before converting if the range matters.
    """
    var res: ILong = 0
    var weight: ILong = 1
    for d in _digits.values() do
      res = res +~ (d.ilong() *~ weight)
      weight = weight *~ _base().ilong()
    end
    if _negative then
      -res
    else
      res
    end


  fun isize(): ISize =>
    """
    Convert to `ISize`. Delegates to `ilong().isize()`. On 64-bit platforms
    `ISize` equals `ILong`, so the only loss is the `ilong()` truncation.
    On 32-bit platforms an additional 32-bit truncation occurs afterward.
    Safe when `ISize.min_value() ≤ this ≤ ISize.max_value()`.
    """
    ilong().isize()


  fun u8(): U8 =>
    """
    Convert to `U8`. Delegates to `ilong().u8()`, so two silent truncations
    can occur: `ilong()` clips to the low 64 bits of `|this|`, then `.u8()`
    further clips to 8 bits (and reinterprets the sign bit).
    Safe only when `0 ≤ this ≤ 255`.
    """
    ilong().u8()


  fun u16(): U16 =>
    """
    Convert to `U16`. Delegates to `ilong().u16()`, so two silent truncations
    can occur: `ilong()` clips to the low 64 bits of `|this|`, then `.u16()`
    further clips to 16 bits. Negative values yield the two's-complement
    reinterpretation. Safe only when `0 ≤ this ≤ 65535`.
    """
    ilong().u16()


  fun u32(): U32 =>
    """
    Convert to `U32`. Delegates to `ilong().u32()`, so two silent truncations
    can occur: `ilong()` clips to the low 64 bits of `|this|`, then `.u32()`
    further clips to 32 bits. Negative values yield the two's-complement
    reinterpretation. Safe only when `0 ≤ this ≤ 2^32 - 1`.
    """
    ilong().u32()


  fun u64(): U64 =>
    """
    Convert to `U64` by accumulating the little-endian base-65536 digits
    using wrapping `U64` arithmetic.

    Each digit contributes `digit × 65536^i`. Because `65536^4 = 2^64 ≡ 0
    (mod 2^64)`, digits at index ≥ 4 contribute 0; only the low 64 bits of
    the magnitude are retained.

    Safe for `0 ≤ this ≤ U64.max_value()` (= 2^64 − 1). Larger positive
    values silently truncate to the low 64 bits. Negative values yield the
    two's-complement `U64` representation of the magnitude.
    """
    var res: U64 = 0
    var weight: U64 = 1
    for d in _digits.values() do
      res = res +~ (d.u64() *~ weight)
      weight = weight *~ _base().u64()
    end
    if _negative then -res else res end


  fun u128(): U128 =>
    """
    Convert to `U128` by accumulating the little-endian base-65536 digits
    using wrapping `U128` arithmetic.

    Each digit contributes `digit × 65536^i`. Because `65536^8 = 2^128 ≡ 0
    (mod 2^128)`, digits at index ≥ 8 contribute 0; only the low 128 bits of
    the magnitude are retained.

    Safe for `0 ≤ this ≤ U128.max_value()` (= 2^128 − 1). Larger positive
    values silently truncate to the low 128 bits. Negative values yield the
    two's-complement `U128` representation of the magnitude.
    """
    var res: U128 = 0
    var weight: U128 = 1
    for d in _digits.values() do
      res = res +~ (d.u128() *~ weight)
      weight = weight *~ _base().u128()
    end
    if _negative then -res else res end


  fun ulong(): ULong =>
    """
    Convert to `ULong`. Delegates to `ilong().ulong()`. On 64-bit platforms
    `ULong` equals `U64`. The `ilong()` truncation applies first (low 64
    bits of `|this|`); negative values then yield the two's-complement
    `ULong` representation. Safe for `0 ≤ this ≤ ULong.max_value()`.
    """
    ilong().ulong()


  fun usize(): USize =>
    """
    Convert to `USize`. Delegates to `ilong().usize()`. On 64-bit platforms
    `USize` equals `U64`. The `ilong()` truncation applies first; negative
    values yield the two's-complement `USize`. Safe for `0 ≤ this ≤ USize.max_value()`.
    """
    ilong().usize()


  fun f32(): F32 =>
    """
    Convert to `F32`. Delegates to `f64().f32()`.

    `F32` has a 24-bit mantissa, so integers are represented exactly only for
    `|this| ≤ 2^24` (= 16 777 216). Larger values are rounded to the nearest
    representable single-precision float. Values exceeding `F32.max_value()`
    (≈ 3.4 x 10^38) become `+Inf` or `-Inf`. No error is raised.
    """
    f64().f32()


  fun f64(): F64 =>
    """
    Convert to `F64` by reconstructing the value from the little-endian
    base-65536 digit array using `F64` accumulation.

    **Precision loss (silent, no error).**
    `F64` has a 53-bit mantissa (IEEE 754 double), so integers are represented
    exactly only for `|this| ≤ 2^53` (≈ 9 x 10^15, roughly three U16 digits).
    Beyond that, the result is rounded to the nearest representable double.
    Values exceeding `F64.max_value()` (≈ 1.8 x 10^308) become `+Inf` or
    `-Inf`. No error is raised.
    """
    var res: F64 = 0
    var weight: F64 = 1.0
    for d in _digits.values() do
      res = res + (d.f64() * weight)
      weight = weight * _base().f64()
    end
    if _negative then
      -res
    else
      res
    end


  fun i8_unsafe(): I8 =>
    """
    Identical to `i8()`. For `MPInt`, the `_unsafe` suffix signals that the
    caller asserts the value fits in `I8`, but no runtime check is removed —
    `MPInt` conversions never panic; they silently truncate.
    """
    i8()


  fun i16_unsafe(): I16 =>
    """
    Identical to `i16()`. See `i16` for overflow behaviour. The caller asserts
    `-32768 ≤ this ≤ 32767`.
    """
    i16()


  fun i32_unsafe(): I32 =>
    """
    Identical to `i32()`. See `i32` for overflow behaviour. The caller asserts
    `-2^31 ≤ this ≤ 2^31 - 1`.
    """
    i32()


  fun i64_unsafe(): I64 =>
    """
    Identical to `i64()`. See `i64` for overflow behaviour. The caller asserts
    `-2^63 ≤ this ≤ 2^63 - 1`.
    """
    i64()


  fun i128_unsafe(): I128 =>
    """
    Identical to `i128()`. See `i128` for overflow behaviour. The caller
    asserts `|this| ≤ 2^127 - 1`.
    """
    i128()


  fun ilong_unsafe(): ILong =>
    """
    Identical to `ilong()`. See `ilong` for overflow behaviour. The caller
    asserts `ILong.min_value() ≤ this ≤ ILong.max_value()`.
    """
    ilong()


  fun isize_unsafe(): ISize =>
    """
    Identical to `isize()`. See `isize` for overflow behaviour. The caller
    asserts `ISize.min_value() ≤ this ≤ ISize.max_value()`.
    """
    isize()


  fun u8_unsafe(): U8 =>
    """
    Identical to `u8()`. See `u8` for overflow behaviour. The caller asserts
    `0 ≤ this ≤ 255`.
    """
    u8()


  fun u16_unsafe(): U16 =>
    """
    Identical to `u16()`. See `u16` for overflow behaviour. The caller asserts
    `0 ≤ this ≤ 65535`.
    """
    u16()


  fun u32_unsafe(): U32 =>
    """
    Identical to `u32()`. See `u32` for overflow behaviour. The caller asserts
    `0 ≤ this ≤ 2^32 - 1`.
    """
    u32()


  fun u64_unsafe(): U64 =>
    """
    Identical to `u64()`. See `u64` for overflow behaviour. The caller asserts
    `0 ≤ this ≤ 2^64 - 1`.
    """
    u64()


  fun u128_unsafe(): U128 =>
    """
    Identical to `u128()`. See `u128` for overflow behaviour. The caller
    asserts `0 ≤ this ≤ U128.max_value()` (= 2^128 − 1).
    """
    u128()


  fun ulong_unsafe(): ULong =>
    """
    Identical to `ulong()`. See `ulong` for overflow behaviour. The caller
    asserts `0 ≤ this ≤ ULong.max_value()`.
    """
    ulong()


  fun usize_unsafe(): USize =>
    """
    Identical to `usize()`. See `usize` for overflow behaviour. The caller
    asserts `0 ≤ this ≤ USize.max_value()`.
    """
    usize()


  fun f32_unsafe(): F32 =>
    """
    Identical to `f32()`. See `f32` for precision-loss behaviour. The caller
    asserts `|this| ≤ 2^24` for an exact result.
    """
    f32()


  fun f64_unsafe(): F64 =>
    """
    Identical to `f64()`. See `f64` for precision-loss behaviour. The caller
    asserts `|this| ≤ 2^53` for an exact result.
    """
    f64()


  //- String ------------------------------------------------------------------

  fun string(): String iso^ =>
    """
    Return the decimal representation of the `MPInt`. No tentatives to
    shorten the string representation using the exponential format (see
    [`from_string`](#from_string)) are done. This can result in long strings 
    but make it easier to know the number of digits from string size.
    """
    if is_zero() then
      return "0".clone()
    end
    var result: String iso = String
    
    // We do all calculations on a local mutable array
    let digits_ref: Array[U16] ref = recover
      let d_tmp = Array[U16](_digits.size())
      for x in _digits.values() do
        d_tmp.push(x)
      end
      d_tmp
    end

    while true do
      if _is_zero(digits_ref) then
        break
      end
      let r = _short_div(digits_ref, 10)
      result.push(r.u8() + '0')
    end
    if _negative then
      result.push('-')
    end
    result.reverse_in_place()
    consume result


  //- Comparisons -------------------------------------------------------------

  fun is_zero(): Bool =>
    """
    Return `true` when the current `MPInt` is zero.
    """
    if _digits.size() == 0 then
      return true
    end
    for d in _digits.values() do
      if d != 0 then
        return false
      end
    end
    true


  fun is_one(): Bool =>
    """
    Return `true` when the current `MPInt` is one.
    """
    (not _negative) and (_digits.size() == 1) and (try _digits(0)? == 1 else false end)


  fun is_minus_one(): Bool =>
    """
    Return `true` when the current `MPInt` is minus one.
    """
    _negative and (_digits.size() == 1) and (try _digits(0)? == 1 else false end)


  fun is_negative(): Bool =>
    """
    Return `true` when the current `MPInt` is strictly negative (< 0).
    """
    _negative and not is_zero()


  fun is_positive(): Bool =>
    """
    Return `true` when the current `MPInt` is strictly positive (> 0).
    """
    (not _negative) and (not is_zero())


  fun compare(that: box->MPInt): Compare =>
    """
    Three-way comparison. Returns `Less`, `Equal`, or `Greater`.

    Satisfies `Comparable[MPInt]`. All comparison operators (`<`, `<=`, `>`,
    `>=`) are derived from this method by the trait's default implementations.
    """
    if _negative != that._negative then
      if is_zero() and that.is_zero() then return Equal end
      return if _negative then Less else Greater end
    end
    let u_lt = _uabs_lt(that)
    let u_eq = _uabs_eq(that)
    if _negative then
      if u_lt then Greater elseif u_eq then Equal else Less end
    else
      if u_lt then Less elseif u_eq then Equal else Greater end
    end


  fun eq(that: box->MPInt): Bool =>
    """
    `this` and `that` are equals when they have the same sign and the same
    digits (they have the same value).
    """
    if _negative != that._negative then
      if is_zero() and that.is_zero() then
        return true
      else
        return false
      end
    end
    _uabs_eq(that)


  fun lt(that: box->MPInt): Bool =>
    """
    Return `true` when `this < that`.
    """
    compare(that) is Less


  fun le(that: box->MPInt): Bool =>
    """
    Return `true` when `this <= that`.
    """
    not (compare(that) is Greater)


  fun ge(that: box->MPInt): Bool =>
    """
    Return `true` when `this >= that`.
    """
    not (compare(that) is Less)


  fun gt(that: box->MPInt): Bool =>
    """
    Return `true` when `this > that`.
    """
    compare(that) is Greater


  fun ne(that: box->MPInt): Bool =>
    """
    Return `true` when `this != that`.
    """
    not (this == that)


  fun abs(): MPInt =>
    """
    Return a new `MPInt` that is the absolute value of `this`.
    """
    if _negative then
      neg()
    else
      _create(false, _digits)
    end


  fun min(that: MPInt): MPInt =>
    """
    Return the minimum of `this` and `that`.
    """
    if this < that then
      _create(_negative, _digits)
    else
      _create(that._negative, that._digits)
    end


  fun max(that: MPInt): MPInt =>
    """
    Return the maximum of `this` and `that`.
    """
    if this > that then
      _create(_negative, _digits)
    else
      _create(that._negative, that._digits)
    end


  fun hash(): USize =>
    """
    Return a hash of the MPInt.
    """
    var h: USize = if _negative then 1 else 0 end
    for d in _digits.values() do
      h = (h * 31) + d.usize()
    end
    h


  fun abs_eq(that: MPInt): Bool =>
    """
    Absolute value equality. This method is equivalent to `this.abs().eq(that.abs())`
    or `this.abs() == that.abs()` but does not require allocation of 2 new `MPInt`
    objects when calculating the absolute values.
    """
    _uabs_eq(that)


  fun abs_ne(that: MPInt): Bool =>
    """
    Absolute value inequality. This method is equivalent to `this.abs().ne(that.abs())`
    or `this.abs() != that.abs()` but does not require allocation of 2 new `MPInt`
    objects when calculating the absolute values.
    """
    not _uabs_eq(that)


  fun abs_lt(that: MPInt): Bool =>
    """
    Absolute value less than. This method is equivalent to `this.abs().lt(that.abs())`
    or `this.abs() < that.abs()` but does not require allocation of 2 new `MPInt`
    objects when calculating the absolute values.
    """
    _uabs_lt(that)


  fun abs_le(that: MPInt): Bool =>
    """
    Absolute value less or equal. This method is equivalent to `this.abs().le(that.abs())`
    or `this.abs() <= that.abs()` but does not require allocation of 2 new `MPInt`
    objects when calculating the absolute values.
    """
    _uabs_lt(that) or _uabs_eq(that)


  fun abs_ge(that: MPInt): Bool =>
    """
    Absolute value greater or equal. This method is equivalent to `this.abs().ge(that.abs())`
    or `this.abs() >= that.abs()` but does not require allocation of 2 new `MPInt`
    objects when calculating the absolute values.
    """
    not _uabs_lt(that)


  fun abs_gt(that: MPInt): Bool =>
    """
    Absolute value greater than. This method is equivalent to `this.abs().gt(that.abs())`
    or `this.abs() > that.abs()` but does not require allocation of 2 new `MPInt`
    objects when calculating the absolute values.
    """
    not (_uabs_lt(that) or _uabs_eq(that))


  //- Helpers -----------------------------------------------------------------

  fun _uabs_eq(that: box->MPInt): Bool =>
    """
    Equality of absolute values. Comparison is done `MPInt` *digit* by *digit*
    taking into account possible non-significant `0`.
    """
    if _digits.size() != that._digits.size() then
      let s1 = _digits.size()
      let s2 = that._digits.size()
      let m = s1.max(s2)
      try
        for k in Range(0, m) do
          let d1 = if k < s1 then _digits(k)? else 0 end
          let d2 = if k < s2 then that._digits(k)? else 0 end
          if d1 != d2 then
            return false
          end
        end
      end
      return true
    end
    try
      for k in Range(0, _digits.size()) do
        if _digits(k)? != that._digits(k)? then
          return false
        end
      end
    end
    true


  fun _uabs_lt(that: box->MPInt): Bool =>
    """
    Unsigned absolute value less than. We must take into account the possible
    non-significant leading `0`.
    """
    if _digits.size() < that._digits.size() then
      for k in Range(_digits.size(), that._digits.size()) do
        try
          if that._digits(k)? != 0 then
            return true
          end
        end
      end
    elseif _digits.size() > that._digits.size() then
      for k in Range(that._digits.size(), _digits.size()) do
        try
          if _digits(k)? != 0 then
            return false
          end
        end
      end
    end
    
    var idx = _digits.size().min(that._digits.size())
    while idx > 0 do
      idx = idx - 1
      try
        if _digits(idx)? < that._digits(idx)? then
          return true
        end
        if _digits(idx)? > that._digits(idx)? then
          return false
        end
      end
    end
    false


  fun tag _is_zero(d: Array[U16] box): Bool =>
    """
    Static zero check for array. Check that the *digits* array contains only 0 digits.
    """
    for x in d.values() do
      if x != 0 then
        return false
      end
    end
    true


  fun tag _normalize(d: Array[U16] ref) =>
    """
    Static normalization for array. Remove the leading non-significant
    0 *digits* in the array that is used by `MPInt` numbers.
    """
    while d.size() > 1 do
      try
        if d(d.size() - 1)? == 0 then
          d.pop()?
        else
          break
        end
      else
        break
      end
    end


  fun tag _short_add(d: Array[U16] ref, b: U32) =>
    """
    Short addition of the single word `b`, added to the least significant *digit*
    of `this`.

    This operation can increaze the size of the `_digits` array by 1.
    """
    var carry = b
    var k: USize = 0
    while (k < d.size()) and (carry > 0) do
      try
        let s = d(k)?.u32() + carry
        d.update(k, s.u16())?
        carry = s >> _base_bits().u32()
      end
      k = k + 1
    end
    if carry > 0 then
      d.push(carry.u16())
    end


  fun tag _short_mul(d: Array[U16] ref, b: U32) =>
    """
    Short multiplication of `this` by the single word `b`.

    This operation can increaze the size of the `_digits` array by 1.
    """
    var carry: U32 = 0
    var k: USize = 0
    while k < d.size() do
      try
        let m = (d(k)?.u32() * b) + carry
        d.update(k, m.u16())?
        carry = m >> _base_bits().u32()
      end
      k = k + 1
    end
    if carry > 0 then
      d.push(carry.u16())
    end


  fun tag _short_div(d: Array[U16] ref, b: U32): U32 =>
    """
    Short division of `this` by the single word `b`, returning the remainder.

    This operation can reduce the size of the `_digits` array.
    """
    var rem_val: U32 = 0
    var k = d.size()
    while k > 0 do
      k = k - 1
      try
        let cur = d(k)?.u32() + (rem_val * _base())
        let q = cur / b
        rem_val = cur % b
        d.update(k, q.u16())?
      end
    end
    _normalize(d)
    rem_val


  fun tag _add_arrays(u: Array[U16] ref, v: Array[U16] box) =>
    """
    Static array addition, adding the content of `v` to `u` and propagating
    the carry over *digits*. In the end, `u <-- u + v`.

    The size of `u` can increase of 1 *digit*.
    """
    var carry: U32 = 0
    var k: USize = 0
    while k < v.size() do
      try
        let s = (if k < u.size() then u(k)? else 0 end).u32() + v(k)?.u32() + carry
        if k < u.size() then
          u.update(k, s.u16())?
        else
          u.push(s.u16())
        end
        carry = s >> _base_bits().u32()
      end
      k = k + 1
    end
    while (k < u.size()) and (carry > 0) do
      try
        let s = u(k)?.u32() + carry
        u.update(k, s.u16())?
        carry = s >> _base_bits().u32()
      end
      k = k + 1
    end
    if carry > 0 then
      u.push(carry.u16())
    end


  fun tag _sub_arrays(u: Array[U16] ref, v: Array[U16] box) =>
    """
    Static array subtraction, subtracting the content of `v` to the *digits*
    of `u`, propagating borrow between *digits*. In the end, `u <-- u - v`.
    The final result is normalized, removing non-significant leading 0.
    """
    var borrow: U32 = 0
    var k: USize = 0
    while k < v.size() do
      try
        let v_val = v(k)?.u32() + borrow
        let u_val = u(k)?.u32()
        if u_val < v_val then
          u.update(k, ((u_val + _base()) - v_val).u16())?
          borrow = 1
        else
          u.update(k, (u_val - v_val).u16())?
          borrow = 0
        end
      end
      k = k + 1
    end
    while (k < u.size()) and (borrow > 0) do
      try
        let u_val = u(k)?.u32()
        if u_val < borrow then
          u.update(k, ((u_val + _base()) - borrow).u16())?
          borrow = 1
        else
          u.update(k, (u_val - borrow).u16())?
          borrow = 0
        end
      end
      k = k + 1
    end
    _normalize(u)


  fun digit_shl(n: USize): MPInt =>
    """
    Shift the *digits* left for `n` *digit* positions. Remember that the *digits*
    of an `MPInt` are in base `base` and are stored in Litle-Endian order internally. So a
    call `digit_shl(n)` multiplies the number by `base^n`.
    """
    if is_zero() then
      return _create(_negative, _digits)
    end
    let d = recover val
      let res = Array[U16].init(0, n + _digits.size())
      for k in Range(0, _digits.size()) do
        try res.update(k + n, _digits(k)?)? end
      end
      res
    end
    _create(_negative, d)


  fun digit_shr(n: USize): MPInt =>
    """
    Shift *digits* right for `n` *digit* positions. Remember that the *digits* of an
    `MPInt` are in base `base` and are stored in Litle-Endian order internally. So a
    call `digit_shr(n)` divides the number by `base^n)`
    """
    if n >= _digits.size() then
      return MPInt.from_ilong(0)
    end
    let d = recover val
      let res = Array[U16].init(0, _digits.size() - n)
      for k in Range(n, _digits.size()) do
        try res.update(k - n, _digits(k)?)? end
      end
      res
    end
    _create(_negative, d)


  fun dump(): String =>
    """
    Hexadecimal debug dump in digit groups used for debug.
    """
    var s = if _negative then "-" else "+" end
    for d in _digits.values() do
      s = s + "_" + Format.int[U16](d, FormatHexBare where width = _base_bits() / 4, fill = '0')
    end
    s
