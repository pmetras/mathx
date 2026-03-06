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


  fun floor_log2(a: U): U =>
    """
    Returns `⌊log₂(a)⌋`, the position of the highest set bit (zero-indexed).
    This is the largest `k` such that `2^k <= a`.

    Returns the bitwidth of `U` as a sentinel when `a == 0` (mathematically
    undefined, since log(0) = -∞).

    Note: this is distinct from `log2`, which computes `popcount(a - 1)` and
    is only meaningful in Gray-code / binary-carry contexts.
    """
    let zero = U.from[USize](0)
    if a == zero then
      return zero.bitwidth()
    end
    (a.bitwidth() - U.from[USize](1)) - a.clz()


  fun is_power_of_two(a: U): Bool =>
    """
    Returns `true` if and only if `a` is an exact power of 2 (including 1 = 2^0).
    Returns `false` for 0.

    Uses the identity that powers of 2 have exactly one bit set:
    `a & (a - 1) == 0`.
    """
    let zero = U.from[USize](0)
    let one = U.from[USize](1)
    (a != zero) and ((a and (a - one)) == zero)


  fun next_power_of_two(a: U): U =>
    """
    Returns the smallest power of 2 that is greater than or equal to `a`.
    Returns 1 for `a <= 1`. Returns 0 on overflow (when no power of 2 >= `a`
    fits the type, e.g. `next_power_of_two(129)` on `U8`).

    Uses the bit-spreading technique: filling all bits below the highest set
    bit of `a - 1`, then adding 1.
    """
    let zero = U.from[USize](0)
    let one = U.from[USize](1)
    if a <= one then return one end
    var v = a - one
    var shift = one
    while shift < v.bitwidth() do
      v = v or (v >> shift)
      shift = shift + shift
    end
    v + one


  fun isqrt(a: U): U =>
    """
    Returns `⌊√a⌋`, the largest integer `k` such that `k * k <= a`.

    Uses Newton's method (integer variant): converges in O(log(a)) steps.
    Returns 0 for `a == 0` and 1 for `a == 1`.
    """
    let zero = U.from[USize](0)
    let one = U.from[USize](1)
    let two = U.from[USize](2)
    if a <= one then return a end
    // Initial overestimate: 2^ceil((bits+1)/2) where bits = floor_log2(a) + 1
    let bits = a.bitwidth() - a.clz()
    var x = one << ((bits + one) / two)
    var x1 = (x + (a / x)) / two
    while x1 < x do
      x = x1
      x1 = (x + (a / x)) / two
    end
    x


  fun parity(a: U): Bool =>
    """
    Returns `true` if `a` has an odd number of set bits (odd parity),
    `false` if it has an even number (even parity).

    Useful in error detection (e.g. parity bits) and Gray-code computations.
    """
    (a.popcount() % U.from[USize](2)) == U.from[USize](1)


  fun ilog(a: U, base: U): U =>
    """
    Returns `⌊log_base(a)⌋`, the largest `k` such that `base^k <= a`.

    Returns 0 for `a == 0` or `base < 2` (undefined cases).

    Computed by repeated division, so its complexity is O(log_base(a)).
    """
    let zero = U.from[USize](0)
    let one = U.from[USize](1)
    if (a == zero) or (base < U.from[USize](2)) then return zero end
    var result = zero
    var v = a
    while v >= base do
      v = v / base
      result = result + one
    end
    result


