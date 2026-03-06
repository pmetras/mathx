// Multi-precision integers

use "debug"
use "format"

use "../assertx"


class MPInt
  """
  A multiple-precision integer that can be used for arbitrary precision
  calculation on integers.

  https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic

  This pure Pony implementation must not be used for cryptographic applications.

  A `MPInt` is represented by an array of `U16`, with the most important digit
  in index `0`. The digits are encoded in base `65536`. Because multiplication
  and division use a Fast Fourier Transform algorithm using `F64` values having
  a 53 bits mantissa, this limits the size of `MPInt` (and `MPFloat`) to around
  1e6 base-65536 digits. This should be enough for most uses of these classes.
  If you need more precision, use the GMP/MPFR implementations of `MPInt` and
  `MPFloat`.
  """

  let _base: U32 = 65536
    """
    The base used to encode the 'digits' of the number. The digits are coded as
    `U16` giving a `2^16` base.
    """


  let _digits: Array[U16] = []
    """
    The arbitrary-precision number encoded in base 65536. The least important
    digit is in index `0`.
    """


  var _negative: Bool = false
    """
    The sign of the integer, `true` when less than `0`.
    """


  new create(s: String = "", base: U32 = 10) ? =>
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
    value of the exponent is an `ILong`.

    The string format uses the following grammar:

    ```BNF
    mpint ::= sign? (digit | separator)+ exponent?
    sign ::= '+' | '-'
    digit ::= '0' | '1' | ... | '9' | 'a' | ... | 'z' | 'A' | ... | 'Z'
    separator ::= '_'
    exponent ::= [ 'e' | 'E' | '@' ] '+'? digit+
    ```

    The whole string `s` must represent an integer number. Leading spaces
    are accepted but not trailing spaces.
    """
    ifdef debug then
      Assert((2 <= base) and (base <= 36), "MPInt.create: The base (" +
            base.string() + ") must be in the range [2..36]", true)?
    end

    // Fast exit: empty string is 0
    if s == "" then
      _digits.push(0)
      return
    end
Debug("Fast exit")

    var i: USize = 0
    var c = s(i)?
    let size = s.size()

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
Debug("Leading space " + i.string() + "/" + c.string())

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
Debug("Sign " + i.string() + "/" + c.string())

    var m: U32 = 0
    while true do
      var n: U32 = 0

      // Separator
      if c == '_' then
        i = i + 1
        c = s(i)?
        continue
      // Exponent character
      elseif (c == '@') or (((c == 'E') or (c == 'e')) and (base <= 10)) then
        i = i + 1
        c = s(i)?
        break
      elseif ('0' <= c) and (c <= '9') then
        n = (c - '0').u32()
      elseif ('a' <= c) and (c <= 'z') then
        n = ((c - 'a') + 10).u32()
      elseif ('A' <= c) and (c <= 'Z') then
        n = ((c - 'A') + 10).u32()
      else
        Debug("MPInt.create: The character '" + c.string() + "' is not " +
              "accepted as a digit in base " + base.string())
        error
      end
Debug("Read digit " + i.string() + "/" + c.string())

      if ((m * base) + n) >= _base then
        (m, let r) = ((m * base) + n).divrem(_base)
Debug("Push r=" + r.string() + " / " + Format.int[U32](r, FormatHexBare))
        _digits.push(r.u16())
Debug("New m=" + m.string())
      else
        m = (m * base) + n
Debug("New m=" + m.string())
      end

      i = i + 1
      if i >= size then
        break
      else
        c = s(i)?
      end

Debug("Read next " + i.string() + "/" + c.string())
    end
Debug("End mantissa " + i.string() + "/" + c.string())

    if i < size then
      // Optional exponent with plus sign
      if c == '+' then
        i = i + 1
        if i >= size then
          Debug("MPInt.create: Missing exponent value")
          error
        else
          c = s(i)?
        end
      end
Debug("Sign exp " + i.string() + "/" + c.string())

      var exp: ILong = 0
      var seen: Bool = false
      while true do
        if c == '_' then
          None
        elseif ('0' <= c) and (c <= '9') then
Debug("Good exp digit " + c.string())
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
Debug("Exp=" + exp.string() + " " + i.string() + "/" + c.string())

      // Process exponent now
      while (exp > 0) and (m < _base) do
        exp = exp - 1
        m = 10 * m
      end
      while m >= _base do
        (exp, let m') = exp.divrem(_base.ilong())
        m = m'.u32()
        _digits.push(m.u16())
        while (exp > 0) and (m < _base) do
          m = 10 * m
          exp = exp - 1
        end
      end
    end
    _digits.push(m.u16())
Debug("Pushed exponent")
    _digits.reverse_in_place()


  new from_ilong(n: ILong) =>
    """
    Create a new arbitrary-precision integer with value `n`.
    """
    var q = n.abs()
    _negative = (q.ilong() != n)
    let base = _base.ulong()
    while q >= base do
      (q, let r) = q.divrem(base)
      _digits.push(r.u16())
    end
    _digits.push(q.u16())


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
      Assert(size > 0, "MPInt.string: [BUG] Non-initialized _digits array", true)?
    end
    if size == 0 then
      return "0".clone()
    end

    var result: String iso = String(5 * size)
    var i: USize = 0
    var m: U32 = 0
    try
      while i < size do
        m = (m * _base) + _digits(i)?.u32()
        repeat
          (m, let r) = m.divrem(10)
          let c: U8 = r.u8() + '0'
          result.push(c)
        until m < 10 end

        i = i + 1
      end
      let c: U8 = m.u8() + '0'
      result.push(c)

      if _negative then
        result.push('-')
      end
    else
      Debug("MPInt.string: Index " + i.string() + " out of bounds [0.." +
            size.string() + ")")
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
    var result: String = if _negative then "- " else "  " end
    try
      while i < size do
        if i != 0 then
          result = result + "_"
        end
        result = result + Format.int[U16](_digits(i)?, FormatHexBare)
        i = i + 1
      else
        result = result + "0"
      end
    else
      Debug("MPInt.dump: Index " + i.string() + " out of bounds [0.." +
            size.string() + ")")
    end
    result



