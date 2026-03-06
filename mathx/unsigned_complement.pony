// Unsigned arithmetic complementary functions


primitive UnsignedComp[U: UnsignedInteger[U] val = USize]
  """
  Unsigned arithmetic functions on unsigned integers.
  """

  fun pow(a: U, n: U): U =>
    """
    Calculate `a^n` (`a` to the power of `n`), that is `a * a * ... n times ... * a`.
    No overflow is done! Use with caution or considere using `pow_partial`.

    Uses [exponentiation by squaring](https://en.wikipedia.org/wiki/Exponentiation_by_squaring).

    See [Prime](prime.pony) for modular exponentiation `pow_mod`.
    """
    let zero = U.from[USize](0)
    let one = U.from[USize](1)
    let two = U.from[USize](2)
    var t = a
    var res = one
    var m = n
    while m > zero do
      if (m % two) == one then
        res = res * t
      end
      t = t * t
      m = m / two
    end
    res


  fun pow_partial(a: U, n: U): U ? =>
    """
    Calculate `a^n` (`a` to the power of `n`), that is `a * a * ... n times ... * a`.
    It raises an error if an overflow occurs.

    Uses [exponentiation by squaring](https://en.wikipedia.org/wiki/Exponentiation_by_squaring).
    """
    let zero = U.from[USize](0)
    let one = U.from[USize](1)
    let two = U.from[USize](2)
    var t = a
    var res = one
    var m = n
    while m > zero do
      if (m % two) == one then
        res = res *? t
      end
      m = m / two
      if m > zero then
        t = t *? t
      end
    end
    res


  fun pow2(n: U): U =>
    """
    Calculate `2^n` (2 at the power of `n`).
    """
    let zero = U.from[USize](0)
    let one = U.from[USize](1)
    if n == zero then
      return one
    end
    one << n


  fun log2(a: U): U =>
    """
    Calculate the logarithm of base 2 of `a`.

    The log2 of an integer `a` is the number of bits set in `a - 1`.

    This function is defined only for `a >= 1`. `log2(0)` returns the bitwidth of
    the type and can be considered an invalid value.
    """
    (a - U.from[USize](1)).popcount()


  fun hamming_distance(a: U, b: U): U =>
    """
    Calculate the Hamming distance between `a` and `b`, that is the minimal
    number of bit changes (`0 -> 1` or `1 -> 0`) required to convert `a` to
    `b`.

    [Hamming distance on Wikipedia](https://en.wikipedia.org/wiki/Hamming_distance)    
    """
    (a xor b).popcount()


