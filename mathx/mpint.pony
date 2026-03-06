// Multi-precision integers

use "debug"

use "../assertx"

use "format"
use "collections"


class val MPInt is Comparable[MPInt val]
  """
  A multiple-precision integer that can be used for arbitrary precision
  calculation on integers.

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

  let _base: U32 = 65536
    """
    The base used to encode the 'digits' of `MPInt` numbers. The digits are
    coded as `U16` giving a `2^16` base.
    """


  let _digits: Array[U16]
    """
    The arbitrary-precision number encoded in base 65536. The least important
    digit is in index `0`.

    Digits are encoded in reverse of traditional occidental practice, with
    the least important first. So number 1234 is encoded as `4321`, with
    `_digits(0)? == 4, _digits(1)? == 3, _digits(2)? == 2, _digits(3) == 1`.
    This is the natural way for many integer operations, where having a carry
    needs only to `push` it to the array.

    Why select an array of `U16` instead of `U64` or `U32` that woulb be best
    fitted for modern CPU architectures? The size of a **digit** is related to
    the size of **flint**, representing integers as floats... I want to use
    a fast multiplication algorithm. Mutliplication of arbitrary-long integers
    involves lots of sum of products:
    `(_digits(i)? * _digits(j - k)?) + (_digits(i + 1) * _digits(j - k - 1)?) + ...` 
    These are convolutions, that can be calculated efficitiently by using FFT
    (Fast Fourier Transforms). So the multiplication is transformed to a dual
    space with FFT and converted back with inverse FFT. But the result must be
    and integer, while the FFT is calculated using floats. In order not to
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


  var _negative: Bool = false
    """
    The sign of the integer, `true` when less than `0`.
    """


  new val create(s: String = "", base: U32 = 10) ? =>
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
      Assert((2 <= base) and (base <= 36), "The base (" + base.string() +
            ") must be in the range [2..36]", true)?
    end

    let guess = s.size() / 3
    _digits = Array[U16](guess)

    // Fast exit: empty string is 0
    if s == "" then
      _digits.push(0)
      return
    end
//Debug("Fast exit")

    var i: USize = 0
    var c = s(i)?
    let size = s.size()
    var tmp = String

    // Leading spaces are accepted
    while c == ' ' do
      i = i + 1
      if i >= size then
        Debug("MPInt.create: Can't create MPInt from empty string value")
        error
      else
        c = s(i)?
      end
    end
//Debug("Leading space " + i.string() + "/" + c.string())

    // Optional sign
    if c == '-' then
      _negative = true
      i = i + 1
      if i >= size then
        Debug("MPInt.create: Can't create MPInt from minus sign only")
        error
      else
        c = s(i)?
      end
    elseif c == '+' then
      i = i + 1
      if i >= size then
        Debug("MPInt.create: Can't create MPInt from plus sign only")
        error
      else
        c = s(i)?
      end
    end
//Debug("Sign " + i.string() + "/" + c.string())

    // Check if exponent is present
    var imantissa = i.isize()
    var iexp: ISize = -1
    if base > 10 then
      iexp = try s.find("@")? else iexp end
    else
      iexp = try s.find("e")? else iexp end
      if iexp < 0 then
        iexp = try s.find("E")? else iexp end
      end
      if iexp < 0 then
        iexp = try s.find("@")? else iexp end
      end
    end
//Debug("Exponent present =" + (iexp >= 0).string() + ", iexp=" + iexp.string())

    if iexp >= 0 then
      // Read exponent
      i = iexp.usize() + 1
      if i >= size then
        Debug("MPInt.create: Missing exponent value")
        error
      else
        c = s(i)?
      end

      // Optional plus sign
      if c == '+' then
        i = i + 1
        if i >= size then
          Debug("MPInt.create: Missing exponent value")
          error
        else
          c = s(i)?
        end
      end
//Debug("Sign exp " + i.string() + "/" + c.string())

      var exp: ILong = 0
      var seen: Bool = false
      while true do
        if c == '_' then
          None
        elseif ('0' <= c) and (c <= '9') then
//Debug("Good exp digit " + c.string())
          exp = (10 * exp) + (c - '0').ilong()
          seen = true
        elseif c == ' ' then
          Debug("MPInt.create: Trailing spaces are not accepted")
          error
        else
          Debug("MPInt.create: The exponent part of the MPInt must be in " +
                "decimal. The character '" + c.string() + "' is invalid")
          error
        end

        i = i + 1
        if i >= size then
          break
        else
          c = s(i)?
        end
      end
      if not seen then
        Debug("MPInt.create: Missing exponent value")
        error
      end
//Debug("Exp=" + exp.string() + " " + i.string() + "/" + c.string())

      // Process exponent now
      tmp = String.create(exp.usize())
      while exp > 0 do
        exp = exp - 1
        tmp.push('0')
      end

      tmp = s.substring(imantissa, iexp) + tmp
    else
      tmp = s.substring(imantissa)
    end

//Debug("Mantissa = '" + tmp + "'")
    let tmp_size = tmp.size()
    i = 0
    while i < tmp_size do
      // We read batches of 3 digits because 36^3 (36 being the maximum base)
      // can fit in a U16 without overflow (see Knuth TAOCP 2.4.4.E).
      var n: U32 = 0
      var j: USize = 0
      var power: U32 = 1
      while (j < 3) and (i < tmp_size) do
        c = tmp(i)?
        // Separator
        if c == '_' then
          if i >= tmp_size then
            break
          else
            i = i + 1
            c = tmp(i)?
            continue
          end
        elseif ('0' <= c) and (c <= '9') then
          n = (n * base) + (c - '0').u32()
        elseif ('a' <= c) and (c <= 'z') then
          n = (n * base) + ((c - 'a') + 10).u32()
        elseif ('A' <= c) and (c <= 'Z') then
          n = (n * base) + ((c - 'A') + 10).u32()
        else
          Debug("MPInt.create: The character '" + c.string() + "' is not " +
                "accepted as a digit in base " + base.string())
          error
        end
        j = j + 1
        power = power * base
        i = i + 1
      end

      // Evaluate ((u_n*b + u_n-1)*b + ... u_1)*b + u_0
      // base^3
      _short_mul(power)
      _short_add(n)
    end


  new val from_ilong(n: ILong = 0) =>
    """
    Create a new arbitrary-precision integer with value `n` (default 0).
    """
    _negative = (n < 0)
    var q = n.abs()
    let base = _base.ulong()
    _digits = Array[U16](n.bytewidth() / 2)
    while q >= base do
      (q, let r) = q.divrem(base)
      _digits.push(r.u16())
    end
    _digits.push(q.u16())


  new iso _from_array(negative: Bool, digits: Array[U16] iso) =>
    """
    Create a new `MPInt` from its sign `negative` and its internal array
    representation `digits`.
    """
    if digits.size() == 0 then
      digits.push(0)
    end
    _negative = negative
    _digits = consume digits


  //- Clone -------------------------------------------------------------------

  fun clone(): MPInt iso^ =>
    """
    Return a copy of the current `MPInt`.
    """
    let size = _digits.size()
    let digits: Array[U16] iso = Array[U16](size)
    var i: USize = 0
    try
      while i < size do
        digits.push(_digits(i)?)
        i = i + 1
      end
    end
    _from_array(_negative, consume digits)


  //- String representations --------------------------------------------------

  fun string(): String iso^ =>
    """
    Return the decimal representation of the `MPInt`. No tentatives to
    shorten the string representation using the exponential format (see
    [`create`](#create)) are done. This can result in long strings but
    make it easier to know the number of digits from string size.
    """
    let size = _digits.size()

    // The case of the non-initialized MPInt
    try
      Assert(size > 0, "Non-initialized `_digits` array", true)?
    end
    if size == 0 then
      return "0".clone()
    end

    // Preallocate the string. A U16 can have 5 digits in decimal representation.
    var result: String iso = String(5 * size)
    // We get a copy because the _short_div has side effects
    let tmp = clone()
    while not tmp.zero() do
      // TODO: Check if we can use group-digit algorithm, dividing by 10000
      // instead of 10, for faster decimal conversion, without requiring
      // bursting memory allocation.
      let r = tmp._short_div(10)
      let c = r.u8() + '0'
      result.push(c)
    end
    if result == "" then
      result.push('0')
    end

    if _negative then
      result.push('-')
    end

    // Now make it human readable...
    result.reverse_in_place()
    consume result


  fun dump(): String =>
    """
    Dump `MPInt` content for debug.
    """
    var i: USize = 0
    let size = _digits.size()
    var result: String = if _negative then "-" else "+" end
    try
      while i < size do
        if i != 0 then
          result = result + "_"
        end
        result = result + Format.int[U16](_digits(i)?, FormatHexBare where width = 4, fill = '0')
        i = i + 1
      else
        result = result + "0"
      end
    else
      Fail("Index (" + i.string() + ") out of bounds [0.." + size.string() + ")")
    end
    result


  //- Short arithmetic --------------------------------------------------------

  fun _low(a: U32): U16 =>
    """
    Get the lower word of `a`.
    """
    a.u16()


  fun _high(a: U32): U16 =>
    """
    Get the high word of `a`.
    """
    a.shr(16).u16()


  fun ref _short_add(b: U32) =>
    """
    Short addition of the single word `b`, added to the least significant digit
    of `this`.

    This operation can increaze the size of the `_digits` array by 1.
    """
    let size = _digits.size()
    var i: USize = 0
    try
      var carry: U32 = b
      while (i < size) and (carry != 0) do
        let tmp = _digits(i)?.u32() + carry
        _digits.update(i, _low(tmp))?
        carry = _high(tmp).u32()
        i = i + 1
      end
      if ((i == size) and (carry != 0)) or (i == 0) then
        _digits.push(_low(carry))
      end
    else
      Fail("Index (" + i.string() + ") out of bounds [0.." + size.string() + ")")
    end


  fun ref _short_mul(b: U32) =>
    """
    Short multiplication of `this` by the single word `b`.

    This operation can increaze the size of the `_digits` array by 1.
    """
    let size = _digits.size()
    var i: USize = 0

    try
      var carry: U32 = 0
      while i < size do
        let tmp = (_digits(i)?.u32() * b) + carry
        _digits.update(i, _low(tmp))?
        carry = _high(tmp).u32()
        i = i + 1
      end
      if (i == size) and (carry != 0) then
        _digits.push(_low(carry))
      end
    else
      Fail("Index (" + i.string() + ") out of bounds [0.." + size.string() + ")")
    end


  fun ref _short_div(b: U32): U32 =>
    """
    Short division of `this` by the single word `b`, returning the remainder.

    This operation can reduce the size of the `_digits` array.
    """
    try
      Assert(b != 0, "Division by 0", true)?
    end

    let size = _digits.size()
    var i: USize = size - 1

    var quotient: U32 = 0
    var remaind: U32 = 0
    var leading: Bool = true
    try
      while i > 0 do
        (quotient, remaind) = (_digits(i)?.u32() + (remaind * _base)).divrem(b)
        if leading and (quotient == 0) then
          // Remove leading 0s from _digits array
          _digits.pop()?
        else
          _digits.update(i, quotient.u16())?
          leading = false
        end
        i = i - 1
      end
      if i == 0 then
        (quotient, remaind) = (_digits(i)?.u32() + (remaind * _base)).divrem(b)
        _digits.update(0, quotient.u16())?
      end
    else
      Fail("Index (" + i.string() + ") out of bounds [0.." + size.string() + ")")
    end
    remaind


  fun ref _normalize() =>
    """
    **Normalize** the `MPInt`, that is remove eventual leading `0`s.
    """
    let size = _digits.size()
    var i: USize = size - 1
    try
      while i > 0 do
        if _digits(i)? == 0 then
          _digits.pop()?
        else
          return
        end
        i = i - 1
      end
    else
      Fail("Index (" + i.string() + ") out of bounds [0.." + size.string() + ")")
    end


  //- Comparisons -------------------------------------------------------------

  fun zero(): Bool =>
    """
    Return `true` when the current `MPInt` is zero.
    """
    let size = _digits.size()
    var i: USize = size
    try
      while i > 0 do
        if _digits(i - 1)? != 0 then
          return false
        end
        i = i - 1
      end
      if i == 0 then
        if _digits(i)? != 0 then
          return false
        end
      end
    else
      Fail("Index (" + i.string() + ") out of bounds [0.." + size.string() + ")")
      return false
    end
    true


  fun eq(that: MPInt): Bool =>
    """
    `this` and `that` are equals when they have the same sign and the same
    digits.

    The internal representation can be different as one can have leading 0s.
    """
    if that._negative != _negative then
      if that.zero() and zero() then
        return true
      else
        return false
      end
    end
    abs_eq(that)


  fun lt(that: MPInt): Bool =>
    """
    Return `true` when `this < that`.
    """
    // Obvious if they have different signs
    if _negative and not that._negative then
      return true
    elseif not _negative and that._negative then
      return false
    end

    // Now, this and that have the same sign
    var i: USize = 0
    let this_size = _digits.size()
    let that_size = that._digits.size()
    try
      if this_size <= that_size then
        i = that_size
        while i > this_size do
          if that._digits(i - 1)? != 0 then
            return not _negative
          end
          i = i - 1
        end
      else
        i = this_size
        while i > that_size do
          if _digits(i - 1)? != 0 then
            return _negative
          end
          i = i - 1
        end
      end

      i = this_size.min(that_size)
      while i > 0 do
        if _digits(i - 1)? == that._digits(i - 1)? then
          i = i - 1
        elseif _negative then
          return _digits(i - 1)? > that._digits(i - 1)?
        else
          return _digits(i - 1)? < that._digits(i - 1)?
        end
      end
      if _negative then
        return _digits(0)? > that._digits(0)?
      else
        return _digits(0)? < that._digits(0)?
      end
    else
      Fail("Index (" + i.string() + ") not in range [0.." + this_size.string() +
           ") or " + "[0.." + that_size.string() + ")")
      false
    end


  //- Absolute value ----------------------------------------------------------
  // Operations on absolute value are usefull when doing in-place operations.
  // They can also be used when one wants to simulate operations on unsigned
  // arbitrary-precision integers.

  fun abs(): MPInt iso^ =>
    """
    Return a new `MPInt` that is the absolute value of `this`.
    """
    let result = clone()
    if _negative then
      result._negative = false
    end
    result


  fun abs_eq(that: MPInt): Bool =>
    """
    `this` and `that` are equals in absolute value when they have the same
    digits.

    The internal representation can be different as we can have leading 0s.
    """
    var i: USize = 0
    let this_size = _digits.size()
    let that_size = that._digits.size()
    try
      if this_size <= that_size then
        i = that_size
        while i > this_size do
          if that._digits(i - 1)? != 0 then
            return false
          end
          i = i - 1
        end
      else
        i = this_size
        while i > that_size do
          if _digits(i - 1)? != 0 then
            return false
          end
          i = i - 1
        end
      end

      i = this_size.min(that_size)
      while i > 0 do
        if _digits(i - 1)? != that._digits(i - 1)? then
          return false
        end
        i = i - 1
      end
      return _digits(0)? == that._digits(0)?
    else
      Fail("Index (" + i.string() + ") not in range [0.." + this_size.string() +
           ") or " + "[0.." + that_size.string() + ")")
      false
    end


  fun abs_ne(that: MPInt): Bool =>
    """
    Return `true` when `this` and `that` are not equal in absolute value:
    `this.abs() != that.abs()`.
    """
    not abs_eq(that)


  fun abs_lt(that: MPInt): Bool =>
    """
    Compare the absolute values of `this` and `that`. Return `true` when
    `this.abs() < that.abs()`.
    """
    var i: USize = 0
    let this_size = _digits.size()
    let that_size = that._digits.size()
    try
      if this_size <= that_size then
        i = that_size
        while i > this_size do
          if that._digits(i - 1)? != 0 then
            return true
          end
          i = i - 1
        end
      else
        i = this_size
        while i > that_size do
          if _digits(i - 1)? != 0 then
            return false
          end
          i = i - 1
        end
      end

      i = this_size.min(that_size)
      while i > 0 do
        if _digits(i - 1)? == that._digits(i - 1)? then
          i = i - 1
        else
          return _digits(i - 1)? < that._digits(i - 1)?
        end
      end
      return _digits(0)? < that._digits(0)?
    else
      Fail("Index (" + i.string() + ") not in range [0.." + this_size.string() +
           ") or " + "[0.." + that_size.string() + ")")
      false
    end


  fun abs_le(that: MPInt): Bool =>
    """
    Compare the absolute values of `this` and `that`. Return `true` when
    `this.abs() <= that.abs()`.
    """
    abs_lt(that) or abs_eq(that)


  fun abs_ge(that: MPInt): Bool =>
    """
    Compare the absolute values of `this` and `that`. Return `true` when
    `this.abs() >= that.abs()`.
    """
    not abs_lt(that)


  fun abs_gt(that: MPInt): Bool =>
    """
    Compare the absolute values of `this` and `that`. Return `true` when
    `this.abs() > that.abs()`.
    """
    not abs_le(that)


  fun abs_compare(that: MPInt): Compare =>
    """
    Compare the absolute values of `this` and `that`. Return `Less` when
    `this.abs() < that.abs()`; `Greater` when `this.abs() > that.abs()`;
    else `Equal` because they have the same value.
    """
    if abs_eq(that) then
      Equal
    elseif abs_lt(that) then
      Less
    else
      Greater
    end


  //- Miscellaneous functions -------------------------------------------------

  fun sign(): Compare =>
    """
    Return the sign of `this`. If `this < 0`, returns `Less`; if `this > 0`,
    returns `Greater`; and returns `Equal` when `this` is `NaN`.
    """
    if zero() then
      Equal
    elseif _negative then
      Less
    else
      Greater
    end


  fun get_base(): ULong =>
    """
    Get the value of the base used internally for encoding the `MPInt`. The
    current encoding use 65536 radix.
    """
    _base.ulong()


  //- Conversions -------------------------------------------------------------

  fun ilong(): ILong ? =>
    """
    Convert the `MPInt` to a `ILong`. Raises an error if it can't fit because
    of overflow.
    """
    var result: ILong = 0
    var i: USize = _digits.size()
    while i > 0 do
      result = (result *? _base.ilong()) +? _digits(i - 1)?.ilong()
      i = i - 1
    end
    if _negative then
      result = -result
    end
    result


  //- Miscellaneous functions -------------------------------------------------

  fun shl(n: USize): MPInt =>
    """
    Shift the current `MPint` by `n` base-digit positions to the left.
    Contrarily to small integers `shl`, the `n` argument is not a `MPInt` and
    the shift is not done on binary digits (bits). The result is a new
    allocated `MPInt`.

    Note that `shl` can be used to generate powers of `base`:
    `MPInt.from_ilong(1) << n` is `base^n< .
    """
    let size = _digits.size()
    let digits: Array[U16] iso = Array[U16].init(0, size + n)
    //TODO: How can we use _digits.copy_to(d, 0, n, size) instead?
    var i: USize = 0
    try
      while i < size do
        digits.update(i + n, _digits(i)?)?
	i = i + 1
      end
    else
      Fail("Index (" + i.string() + ") out of bounds [0.." +
	   (size + n).string() + ")")
    end
    _from_array(_negative, consume digits)


  fun shr(n: USize): MPInt =>
    """
    Shift the current `MPInt` by `n` base-digit positions to the right. The
    least important digit is lost in the operation.
    Contrarily to small integers `shr`, the `n` argument is not a `MPInt` and
    the shift is not done on bits (binary digits). The result is a new
    allocated `MPInt`. If no more base-digits are available in the `MPInt`,
    the result is 0.
    """
    let size = _digits.size()
    if n >= size then
      MPInt.from_ilong(0)
    else
      let digits: Array[U16] iso = Array[U16].init(0, size - n)
      // TODO: How can we use _digits.copy(digits, n, 0, size - n) instead?
      var i: USize = n
      try
        while i < size do
	  digits.update(i - n, _digits(i)?)?
	  i = i + 1
	end
      else
	Fail("Index (" + i.string() + ") out of bounds [" +
	     (size - n).string() + "..0)")
      end
      _from_array(_negative, consume digits)
    end


  //- Arithmetic --------------------------------------------------------------

  fun neg(): MPInt =>
    """
    Return a new `MPInt` with value `-this`.
    """
    let size = _digits.size()
    let digits: Array[U16] iso = Array[U16](size)
    var i: USize = 0

    try
      while i < size do
        digits.push(_digits(i)?)
	i = i + 1
      end
    else
      Fail("Index (" + i.string() + ") out of range [0.." + size.string() + ")")
    end
    _from_array(not _negative, consume digits)


  fun add(that: MPInt): MPInt =>
    """
    Return a new `MPInt` with the value `this + that`. The size of the result
    can increase of one `_base` digit if there is a carry.

    Depending on the signs of `this` and `that`, the function performs an
    addition or a substraction on the unsigned internal representation.

    Performance and memory size is O(N) where N is the size of the integers.
    """
    if _negative == that._negative then
      // this and that have the same sign: that's a true addition of unsigned
      // numbers.
      let result = clone()
      try
        result._add_in_place(that)?
      else
        Fail("Bug in _add_in_place algorithm")
      end
      result

    else
      // this and that have different signs: that's in fact a substraction.
      // When this is larger than that in absolute value, it gives its sign and
      // the number of digits of the result.

      // TODO: To be optimized to try to prevent allocation for temporaries.
      if this.abs_gt(that) then
        let result = abs()
        let that_abs = that.abs()
        try
          result._sub_in_place(consume that_abs)?
        else
          Fail("Bug in _sub_in_place algorithm - 1")
        end
        // this gives the sign
        result._negative = _negative
        result
      else
        let result = that.abs()
        let this_abs = abs()
        try
          result._sub_in_place(consume this_abs)?
        else
          Fail("Bug in _sub_in_place algorithm - 2")
        end
        // that gives the sign
        result._negative = that._negative
        result
      end
    end


  fun sub(that: MPInt): MPInt =>
    """
    Return a new `MPInt` whose value is equal to the substraction `this - that`.

    Depending on the signs of `this` and `that`, it can be in fact an addition
    of their absolute values.

    Performance and memory size is O(N) where N is the size of the integers 3/4
    of the time, and memory is O(2N) 1/4 of the time because of allocation of a
    temporary.
    """
    if _negative == that._negative then
      // this and that have the same sign: that's a true substraction of unsigned
      // numbers.
      if this.abs_gt(that) then
        // A - B = A - B
        // this gives the sign
        let result = clone()
        try
          result._sub_in_place(that)?
        else
          Fail("Bug in _sub_in_place algorithm - 1")
        end
        result
      
      else
        // A - B = -(B - A)
        // that gives the sign
        let result = that.clone()
        try
          // TODO: See if clone can be avoided
          result._sub_in_place(clone())?
        else
          Fail("Bug in _sub_in_place algorithm - 2")
        end
        result._negative = not result._negative
        result
      end

    else
      // this and that have different signs: it's in fact an addition
      // A - B = A + (-B) = -(-A + B)
      let result = clone()
      result._negative = not result._negative
      try
        result._add_in_place(that)?
      else
        Fail("Bug in _add_in_place algorithm")
      end
      result._negative = not result._negative
      result
    end


  fun ref _add_in_place(that: MPInt) ? =>
    """
    Add in-place the value of `that` to the current `MPInt`. Both integers must
    have the same sign else an error is raised. The size of `this` can increase
    by 1 digit depending if there is a carry.

    See Knuth 4.3.1.A
    """
    Assert(_negative == that._negative,
           "`this` and `that` must have the same sign for the _add_in_place operation", true)?

    let this_size = _digits.size()
    let that_size = that._digits.size()
    let size = this_size.min(that_size)
    var i: USize = 0
    _digits.reserve(this_size.max(that_size) + 1)

    try
      var carry: U32 = 0
        while i < size do
//Debug("Sum1 i=" + i.string())
          let s = _digits(i)?.u32() + that._digits(i)?.u32() + carry
          carry = _high(s).u32()
//Debug("Sum1     d(i)=" + _digits(i)?.string() + ", ta(i)=" + that._digits(i)?.string())
          _digits.update(i, _low(s))?
//Debug("Sum1 now d(i)=" + _digits(i)?.string() + ", carry=" + carry.string())
          i = i + 1
        end

        if this_size < that_size then
          while i < that_size do
//Debug("Sum2 i=" + i.string())
            let s = that._digits(i)?.u32() + carry
//Debug("Sum2     d(i)=0, ta(i)=" + that._digits(i)?.string())
            carry = _high(s).u32()
            _digits.push(_low(s))
//Debug("Sum2 now d(i)=" + _digits(i)?.string() + ", carry=" + carry.string())
            i = i + 1
          end
        else
          while i < this_size do
//Debug("Sum3 i=" + i.string())
            let s = _digits(i)?.u32() + carry
//Debug("Sum3     d(i)=" + _digits(i)?.string() + ", ta(i)=0")
            carry = _high(s).u32()
            _digits.update(i, _low(s))?
//Debug("Sum3 now d(i)=" + _digits(i)?.string() + ", carry=" + carry.string())
            i = i + 1
          end
	  if carry != 0 then
//Debug("Sum3 add carry=" + _low(carry).string())
	    _digits.push(_low(carry))
	  end
        end
    else
      Fail("Index (" + i.string() + ") out of range [0.." + this_size.string() +
           ") or [0.." + that_size.string() + ")")
    end


  fun ref _sub_in_place(that: MPInt) ? =>
    """
    Substract `that` from the current `MPInt`. Assumes that `this` is larger
    than `that` and that the numbers have the same sign else an error is
    raised.

    See Knuth TAOCP 4.3.1.S
    """
    Assert(_negative == that._negative,
           "`this` and `that` must have the same sign for the _sub_in_place operation", true)?
    Assert(this.abs_ge(that),
           "`this` must be larger or equal than `that` in absolute value for " +
           "_sub_in_place operation", true)?
//Debug("this.abs_ge(that)=" + (this.abs_ge(that)).string())
//Debug("this=" + this.string() + " | " + this.dump() + ", that=" + that.string() + " | " + that.dump())

    let this_size = _digits.size()
    let that_size = that._digits.size()

    // If the numbers are normalized (i.e. without leading 0s, then the size of
    // this is larger than the size of that.
    Assert(this_size >= that_size, "Integers are not normalized", true)?
    // So this should be useless
    _digits.reserve(this_size.max(that_size))

    let size = this_size.min(that_size)
    var i: USize = 0

    try
      var carry: U32 = 0
      var s: U32 = 0
//Debug("size=" + this_size.string() + ", tasize=" + that_size.string())
      while i < size do
//Debug("Sub1 i=" + i.string() + ", d(i)=" + _digits(i)?.string() + ", ta(i)=" + that._digits(i)?.string())
        if (that._digits(i)?.u32() + carry) > _digits(i)?.u32() then
          s = (_digits(i)?.u32() + _base) - (that._digits(i)?.u32() + carry)
          carry = 1
        else
          s = _digits(i)?.u32() - (that._digits(i)?.u32() + carry)
          carry = 0
        end
        _digits.update(i, _low(s))?
//Debug("Sub1 now d(i)-ta(i)=" + _digits(i)?.string() + ", carry=" + carry.string())
        i = i + 1
      end

      while i < this_size do
//Debug("Sub2 i=" + i.string() + ", d(i)=" + _digits(i)?.string() + ", carry=" + carry.string())
        if carry > _digits(i)?.u32() then
          s = (_digits(i)?.u32() + _base) - carry
          carry = 1
        else
          s = _digits(i)?.u32() - carry
          carry = 0
        end
        _digits.update(i, _low(s))?
//Debug("Sub2 now d(i)=" + _digits(i)?.string() + ", carry=" + carry.string())
        i = i + 1
      end
      ifdef debug then
        try
          Assert(carry == 0, "BUG in algorithm: `this` (" + this.string() +
          ") is not greater or equal to `that` (" + that.string() + ")", true)?
        end
      end
    else
      Fail("Index (" + i.string() + ") out of range [0.." + this_size.string() +
            ") or [0.." + that_size.string() + ")", true)
    end
    // We remove the leading 0s that could have appeared
    _normalize()


  fun mul(that: MPInt): MPInt =>
    """
    Multiply `this` by `that` and return a new `MPInt`. The number of digits
    of the result is `this.size() + that.size()` but when of the operands is
    zero.

    This multiplication algorithm is the traditional way learnt at school,
    with performance in O(N * M) and size in O(N + M).

    See Knuth TAOCP 4.3.1.M
    """
    // Multiplication by 0
    if zero() or that.zero() then
      return MPInt.from_ilong(0)
    end

    let this_size = _digits.size()
    let that_size = that._digits.size()
    let digits: Array[U16] iso = Array[U16].init(0, this_size + that_size)

    var i: USize = 0
    var j: USize = 0
    try
      i = 0
      while i < this_size do
        var carry: U32 = 0
        j = 0
        while j < that_size do
          let m = (_digits(i)?.u32() * that._digits(j)?.u32()) +
                  digits(i + j)?.u32() + carry
          carry = _high(m).u32()
          digits.update(i + j, _low(m))?
          j = j + 1
        end
        digits.update(i + that_size, _low(carry))?
        i = i + 1
      end
    else
      Fail("Index i (" + i.string() + ") must be in range [0.." +
           this_size.string() + ") and j (" + j.string() +
	   ") must be in range [0.." + that_size.string() + ")")
    end

    let negative = ((_negative and not that._negative) or
                    (not _negative and that._negative))
    _from_array(negative, consume digits)


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
    let this_size = _digits.size()
    let that_size = that._digits.size()

    // End recursion and use other multiplication algorithm when we have
    // small numbers.
    if (this_size <= 128) or (that_size <= 128) then
      mul(that)
    else
      // Now we do the divide and conquer Karatsuba
      let size = if this_size > that_size then this_size else that_size end
      let half = size / 2

      // Most of the time, these loops will raise an error for out of bound
      // index. It won't only when size if odd and this and that have the same
      // size.
      let this_l: Array[U16] iso = Array[U16](half)
      try
        for i in Range(0, half) do
	  this_l.push(_digits(i)?)
	end
      end
      let this_low = recover val _from_array(_negative, consume this_l) end

      let this_h: Array[U16] iso = Array[U16](size - half)
      try
        for i in Range(0, half + 1) do
	  this_h.push(_digits(i + half)?)
	end
      end
      let this_high = recover val _from_array(_negative, consume this_h) end

      let that_l: Array[U16] iso = Array[U16](half)
      try
        for i in Range(0, half) do
	  that_l.push(that._digits(i)?)
	end
      end
      let that_low = recover val _from_array(that._negative, consume that_l) end

      let that_h: Array[U16] iso = Array[U16](size - half)
      try
        for i in Range(0, half + 1) do
	  that_h.push(that._digits(i + half)?)
	end
      end
      let that_high = recover val _from_array(that._negative, consume that_h) end

      let z2 = this_high.karatsuba_mul(that_high)
      let z0 = this_low.karatsuba_mul(that_low)
      let z1 = (this_low + this_high).karatsuba_mul(that_low + that_high) - z2 - z0

        (z2 << (2 * half)) + (z1 << half) + z0
    end  


  fun fast_mul(that: MPInt): MPInt =>
    """
    Fast multiplication of `this` and `that` using FFT. The result has a size
    of `this.size() + that.size()`.

    See Knuth TAOCP 4.3.3.C
    """
    // Multiplication by 0
    if zero() or that.zero() then
      return MPInt.from_ilong(0)
    end

    let this_size = _digits.size()
    let that_size = that._digits.size()
    let max_size = this_size.max(that_size)
    let size = this_size + that_size
    // The minimal power of 2 for the Fourier transforms
    let pow2 = (max_size + 1).next_pow2()

    var i: USize = 0
    let digits: Array[U16] iso = Array[U16].init(0, size)

    try
      // Convert to arrays of F64
      var a: Array[Complex[F64]] = Array[Complex[F64]](pow2)
      i = 0
      while i < this_size do
        a.push(Complex[F64](_digits(i)?.f64()))
        i = i + 1
      end
      while i < pow2 do
        a.push(Complex[F64])
        i = i + 1
      end

a = cdump("astart=", consume a)

      var b: Array[Complex[F64]] = Array[Complex[F64]](pow2)
      i = 0
      while i < that_size do
        b.push(Complex[F64](that._digits(i)?.f64()))
        i = i + 1
      end
      while i < pow2 do
        b.push(Complex[F64])
        i = i + 1
      end

b = cdump("bstart=", consume b)

      // Convolution
      a = FFT.fourier(a)
      b = FFT.fourier(b)

a = cdump("fouriera=", consume a)
b = cdump("fourierb=", consume b)
      // Multiply in dual space with complex multiplication
      i = 0
      while i < max_size do
        b.update(i, a(i)? * b(i)?)?
        i = i + 1
      end

b = cdump("fouriermult=", consume b)

      // Inverse FFT
      b = FFT.fourier(b, true)

b = cdump("invfourmult=", consume b)

      // Convert back to integers
      let base: F64 = _base.f64()
      var carry: F64 = 0.0
      i = 0
      while i < size do
        // Add 0.5 to round conversion
        let temp = b(i)?.real() + carry + 0.5

        ifdef debug then
          Assert(temp > 0.0, "Bug in algorithm: FFT result b(" + i.string() +
                ")=" + b(i)?.real().string() + " must be positive", true)?
        end

        if temp >= base then
          carry = (temp / base).trunc()

          ifdef debug then
            Assert(temp >= (carry * base), "Value " + temp.string() +
                  " - " + (carry * base).string() + " must be positive", true)?
          end
          let m = temp - (carry * base)
          digits.update(i, m.u16())?
        else
          carry = 0.0
          digits.update(i, temp.u16())?
        end

Debug("b(" + i.string() + ")=" + b(i)?.string() + ", temp=" + temp.string() + ", carry=" + carry.string() + ", d(" + i.string() + ")=" + (temp - (carry * base)).string() + " / " + (temp - (carry * base)).u16().string())

        i = i + 1
      end

let d: Array[U16] = []
for i' in Range(0, digits.size()) do
  d.push(digits(i')?)
end
Debug("fastmult=" + idump(d))
//b = cdump("bresult=", consume b)

      ifdef debug then
        Assert(carry < base, "We can't have a carry (" + 
              carry.string() + ") higher than " + base.string(), true)?
      end
Debug("carry=" + carry.string() + " / " + carry.u16().string())

      let temp = carry.u16()
      if temp != 0 then
        try digits.update(size - 1, digits(size - 1)? + temp)? else Debug("Index error in update") end
      else
        // Non-significative 0
        try digits.pop()? else Debug("Error in pop") end
      end


    else
      Fail("Index (" + i.string() + ") out of bounds [0.." + size.string() + ")")
    end

    // Create result
    let negative = ((_negative and not that._negative) or
                    (not _negative and that._negative))
    _from_array(negative, consume digits)





/*
  fun fast_mul_backup(that: MPInt): MPInt =>
    """
    Fast multiplication of `this` and `that` using FFT. The result has a size
    of `this.size() + that.size()`.

    See Knuth TAOCP 4.3.3.C
    """
    // Multiplication by 0
    if zero() or that.zero() then
      return MPInt.from_ilong(0)
    end

    let this_size = _digits.size()
    let that_size = that._digits.size()
    let max_size = this_size.max(that_size)
    let size = this_size + that_size
    // The minimal power of 2 for the Fourier transforms
    let pow2 = (max_size + 1).next_pow2()

    var i: USize = 0
    let digits: Array[U16] iso = Array[U16].init(0, size)

    try
      // Convert to arrays of F64
      let a = Array[F64].init(0.0, pow2)
      i = 0
      while i < this_size do
        a.update(i, _digits(i)?.f64())?
        i = i + 1
      end

Debug("astart=" + fdump(a))

      let b = Array[F64].init(0.0, pow2)
      i = 0
      while i < that_size do
        b.update(i, that._digits(i)?.f64())?
        i = i + 1
      end

Debug("bstart=" + fdump(b))

      // Convolution
      FFT.fourier_real(a)
      FFT.fourier_real(b)

Debug("fouriera=" + fdump(a))
Debug("fourierb=" + fdump(b))
      // Multiply in dual space with complex multiplication
      // Special case of indexes 0 and 1
      b.update(0, b(0)? * a(0)?)?
      b.update(1, b(1)? * a(1)?)?
      i = 2
      while i < max_size do
        let temp = b(i)?
        b.update(i, (temp * a(i)?) - (b(i + 1)? * a(i + 1)?))?
        b.update(i + 1, (temp * a(i + 1)?) + (b(i + 1)? * a(i)?))?
        i = i + 2
      end

Debug("fouriermult=" + fdump(b))

      // Inverse FFT
      FFT.fourier_real(b, true)

Debug("invfourmult=" + fdump(b))

      // Convert back to integers
      let base: F64 = _base.f64()
      var carry: F64 = 0.0
      i = 0
      while i < (size - 1) do
        // Add 0.5 to round conversion
        let temp = b(i)? + carry + 0.5
        ifdef debug then
          Assert(temp > 0.0, "Bug in algorithm: FFT result b(" + i.string() +
                ")=" + b(i)?.string() + " must be positive", true)?
        end
        carry = (temp / base).trunc()

        ifdef debug then
          Assert(temp >= (carry * base), "Value " + temp.string() +
                " - " + (carry * base).string() + " must be positive", true)?
        end

Debug("b(" + i.string() + ")=" + b(i)?.string() + ", temp=" + temp.string() + ", carry=" + carry.string() + ", d(" + i.string() + ")=" + (temp - (carry * base)).string() + " / " + (temp - (carry * base)).u16().string())
        let m = temp - (carry * base)
        b.update(i, m)?
        i = i + 1
      end

Debug("bresult=" + fdump(b))

      ifdef debug then
        Assert(carry < base, "We can't have a carry (" + 
              carry.string() + ") higher than " + base.string(), true)?
      end
//Debug("d(" + (size - 1).string() + ")=" + carry.string() + " / " + carry.u16().string())

      // Prepare result
      i = 0
      while i < (size - 1) do
        digits.update(i, b(i)?.u16())?
        i = i + 1
      end
      let temp = carry.u16()
      if temp != 0 then
        digits.update(size - 1, temp)?
      else
        // Non-significative 0
        digits.pop()?
      end

let d: Array[U16] = []
for i' in Range(0, digits.size()) do
  d.push(digits(i')?)
end
Debug("fastmult=" + idump(d))

    else
      Fail("Index (" + i.string() + ") out of bounds [0.." + size.string() + ")")
    end

    // Create result
    let negative = ((_negative and not that._negative) or
                    (not _negative and that._negative))
    _from_array(negative, consume digits)
*/




fun cdump(t: String, a: Array[Complex[F64]]): Array[Complex[F64]] =>
  var s = ""
  for i in Range(0, a.size()) do
    s = s + " | " + try a(i)?.string() else "***" end
  end
  Debug(t + s)
  a

fun fdump(a: Array[F64]): String =>
  var s = ""
  for i in Range(0, a.size()) do
    s = s + " | " + try a(i)?.string() else "***" end
  end
  s

fun idump(a: Array[U16]): String =>
  var s = ""
  for i in Range(0, a.size()) do
    s = s + " | " + try a(i)?.string() else "***" end
  end
  s
