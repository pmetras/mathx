"""
# Prime numbers

Prime numbers functionalities library, inspired by https://github.com/adri326/pony-primes

For information about prime numbers, visit the [PrimePages](https://primes.utm.edu/index.html)
"""

use "../assertx"
use "../bitsx"
use "random"


primitive Prime[A: UnsignedInteger[A] val = USize]
  """
  Primitive for regular prime operations: primality test, coprimality test, GCD
  and LCM, next prime and prime factorization
  """

  fun is_prime(num: A): Bool =>
    """
    Basic primality test, returns `true` if the given number is indeed prime. This
    test uses a simple method by [trial division](https://en.wikipedia.org/wiki/Primality_test#Simple_methods).

    **Usage:**
    ```pony
      if Prime.is_prime[U128](232862364312348451) then
        env.out.print("We got ourselves quite a big prime!")
      else
        env.out.print("Not a prime number :C")
      end
    ```
    """
    if (num == A.from[USize](2)) or (num == A.from[USize](3)) then
      return true
    end
    if (num <= A.from[USize](1)) or
       ((num % A.from[USize](2)) == A.from[USize](0)) or
       ((num % A.from[USize](3)) == A.from[USize](0)) then
      return false
    end
    
    var current: A = A.from[USize](5)
    while (current * current) <= num do
      if ((num % current) == A.from[USize](0)) or
         ((num % (current + A.from[USize](2))) == A.from[USize](0)) then
        return false
      end
      current = current + A.from[USize](6)
    end
    true


  fun _get_first_base(i: USize): A =>
    """
    Get the predefined base values. These predefined bases are sufficient to
    guarantee primality for numbers inferior to `U128.max_value()`.
    The match pattern should force LLVM to create a fix jump table (no memory
    allocation) instead of creating a Pony `Array` (allocating memory at each
    call).
    """
    match i
    | 0 => A.from[USize](2)
    | 1 => A.from[USize](3)
    | 2 => A.from[USize](5)
    | 3 => A.from[USize](7)
    | 4 => A.from[USize](11)
    | 5 => A.from[USize](13)
    | 6 => A.from[USize](17)
    | 7 => A.from[USize](19)
    | 8 => A.from[USize](23)
    | 9 => A.from[USize](29)
    | 10 => A.from[USize](31)
    | 11 => A.from[USize](37)
    | 12 => A.from[USize](41)
    else
      A.from[USize](0)
    end


  fun _get_base(i: USize, num: A): A =>
    """
    Get a base for the Rabin-Miller algorithm. According to the algorithm, it
    should be a random number, but some selected numbers are used first as they
    can guarantee primality. It returns the `i`-th base inferior to `num`.
    """
    if i <= 12 then // Maximum bases returned by _get_first_base
      let b = _get_first_base(i)
      if b >= num then
        // We have tried all interesting base values lower than `num`. It
        // must be a prime
        A.from[USize](0)
      else
        b
      end
    else
      let two = A.from[USize](2)
      let four = A.from[USize](4)
      let rand = Rand
      let b: A = two + (A.from[U64](rand.next()) % (num - four))
      b
    end


  fun is_probably_prime(num: A, k: USize = 13): Bool =>
    """
    Test if `num` seems to be a prime number, using probablistic primality tests.
    Return `false` when `num` is composed. Return `true` if `num` is probably
    prime.
    Parameter `k` determines the accuracy of the test. The greater the number of
    rounds, the more accurate the result. Default is probably more than sufficient
    for common `A` types.
    
    This version implements the
    [Miller-Rabin probabilistic primality test](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test).

    For type `A` having less than 128 bits width, this probabilistic test can
    find primes deterministically, but the algorithm implemented should work
    with larger types.
    """
    let zero = A.from[USize](0)
    let one = A.from[USize](1)
    let two = A.from[USize](2)
    let three = A.from[USize](3)
    let four = A.from[USize](4)

    // Handle small edge cases: 0 and 1 are not prime; 2 and 3 are prime
    if num < two then
      return false
    end
    if (num == two) or (num == three) then
      return true
    end
    if (num % two) == zero then
      return false
    end
    
    // Find the powers of 2 out of num - 1
    let s = (num - one).ctz()
    let d = (num - one) >> s
    
    var i: USize = 0
    while i < k do
      // Find the base a
      let a = _get_base(i, num)

      // We have tried all bases less than num; we have a prive
      if a == zero then
        return true
      end

      // Calculate x = a^d mod num
      var x = Modular[A].pow_mod(a, d, num)

      if (x == one) or (x == (num - one)) then
        // `a` in not a Miller-Rabin witness
        i = i + 1
        continue
      end

      var j = one
      while j < s do
        x = Modular[A].mul_mod(x, x, num)

        if x == (num - one) then
          // `a` is not a Miller-Rabin witness
          break
        end

        j = j + one
      end
      if j == s then
        // `num` is composite
        return false
      end
      
      // Try another base `a`
      i = i + 1
    end
    true


  fun is_coprime(a: A, b: A): Bool =>
    """
    Coprimality test: returns `true` if the only common divisor between `a` and `b` is `1`

    **Usage:**
    ```pony
      env.out.print(Prime.is_coprime(32, 33).string()) // should be 'true'
    ```
    """
    Modular[A].gcd2(a, b) == A.from[USize](1)


  fun next_prime(num: A): A =>
    """
    Returns the prime which follows `num`.
    """
    PrimeIterator[A].start_at(num).next()


  fun prime_factors(num: A): Array[A] ref =>
    """
    Returns an array of prime numbers representing all the prime factors of `num`.
    Prime factorization consists of splitting up a number into a series of prime
    numbers which all multiplied together give you that number.
    There only exists one prime factorization for every number (this is a
    property of the prime numbers).

    **Note:** `24 = 2^3 * 3` - in the case of a prime number occuring several
    times, they will simply be listed this amount of times in the result
    array: `[2; 2; 2; 3]`.
    
    `prime_factors(0)` returns `[]` instead of an error.
    `prime_factors(1)` returns `[]`
    """
    match num
    | A.from[USize](0) => return [] // By definition
    | A.from[USize](1) => return []
    | A.from[USize](2) => return [A.from[USize](2)]
    | A.from[USize](3) => return [A.from[USize](3)]
    end
    
    if Prime[A].is_prime(num) then
      return [num]
    end
    
    var new_num = num
    let iterator = PrimeIterator[A]
    var res = Array[A]
    for divisor in iterator do
      while (new_num % divisor) == A.from[USize](0) do
        new_num = new_num / divisor
        res.push(divisor)
      end
      if new_num == A.from[USize](1) then
        break
      end
      // If divisor > sqrt(new_num), the remaining new_num must be prime
      if (new_num / divisor) < divisor then
        res.push(new_num)
        break
      end
    end
    res


  fun prime_factors_unique(num: A): Array[A] ref =>
    """
    Returns an array of the distinct prime factors of `num`, each listed once.

    Unlike `prime_factors`, repeated factors are deduplicated:
    `prime_factors_unique(12)` returns `[2; 3]` (not `[2; 2; 3]`).

    `prime_factors_unique(0)` and `prime_factors_unique(1)` return `[]`.
    """
    match num
    | A.from[USize](0) => return []
    | A.from[USize](1) => return []
    | A.from[USize](2) => return [A.from[USize](2)]
    | A.from[USize](3) => return [A.from[USize](3)]
    end

    if Prime[A].is_prime(num) then
      return [num]
    end

    var new_num = num
    let iterator = PrimeIterator[A]
    var res = Array[A]
    for divisor in iterator do
      if (new_num % divisor) == A.from[USize](0) then
        res.push(divisor)
        while (new_num % divisor) == A.from[USize](0) do
          new_num = new_num / divisor
        end
      end
      if new_num == A.from[USize](1) then
        break
      end
      if (new_num / divisor) < divisor then
        res.push(new_num)
        break
      end
    end
    res


  fun euler_totient(n: A): A =>
    """
    Euler's totient function `φ(n)`: the count of integers in `[1, n]` that
    are coprime with `n`.

    Uses the formula `φ(n) = n * ∏ (p − 1) / p` over all distinct prime
    factors `p` of `n`. Each factor is applied as `result = result / p * (p − 1)`
    to stay in integer arithmetic.

    Returns 0 for `n == 0`. Returns 1 for `n == 1`.

    Fundamental in number theory: appears in RSA, discrete logarithms, and
    the multiplicative group of integers modulo n.
    """
    let zero = A.from[USize](0)
    let one = A.from[USize](1)
    if n == zero then
      return zero
    end
    var result = n
    for p in Prime[A].prime_factors_unique(n).values() do
      result = (result / p) * (p - one)
    end
    result


  fun radical(n: A): A =>
    """
    Returns the radical of `n`: the product of the distinct prime factors of `n`.

    Examples: `radical(12) = 2 * 3 = 6`, `radical(8) = 2`, `radical(30) = 30`.

    `radical(p^k) = p` for any prime `p` and `k >= 1`.

    Returns 1 for `n <= 1`.
    """
    let one = A.from[USize](1)
    if n <= one then
      return one
    end
    var result = one
    for p in Prime[A].prime_factors_unique(n).values() do
      result = result * p
    end
    result


  fun is_squarefree(n: A): Bool =>
    """
    Returns `true` if `n` is squarefree: no prime factor appears more than once.
    Equivalently, `n == radical(n)`.

    Returns `true` for `n <= 1` (vacuously squarefree).

    Examples: `is_squarefree(30) = true` (30 = 2*3*5), `is_squarefree(12) = false`
    (12 = 2²*3).
    """
    let one = A.from[USize](1)
    if n <= one then
      return true
    end
    Prime[A].radical(n) == n


  fun prev_prime(num: A): A =>
    """
    Returns the largest prime strictly less than `num`.

    Returns 0 if no such prime exists (i.e. `num <= 2`).
    """
    let zero = A.from[USize](0)
    let one = A.from[USize](1)
    let two = A.from[USize](2)
    if num <= two then
      return zero
    end
    if num == A.from[USize](3) then
      return two
    end
    // Start from the nearest odd number strictly below num
    var candidate = if (num % two) == zero then num - one else num - two end
    while not is_prime(candidate) do
      if candidate <= two then
        return zero
      end
      candidate = candidate - two
    end
    candidate


class PrimeIterator[A: UnsignedInteger[A] val = USize] is Iterator[A]
  """
  The prime iterator returns every primes up to a number (by default the maximum
  value of the type). Note that it will return the prime following the limit
  value.

  **Usage:**
  ```pony
    let iterator = PrimeIterator[U32](100)
    for prime in iterator do
      env.out.print(prime.string())
    end
    // will print: 2, 3, 5, 7, ..., 101
  ```
  """

  var _last: A
    """
    Last prime found.
    """
    
  let _limit: A
    """
    The upper limit of the iterator.
    """


  new create(limit: A = A.max_value()) =>
    """
    Create a new prime iterator that generates prime numbers up to the next prime
    after `limit`. The default limit is the maximum value of the type.
    """
    _last = A.from[USize](1)
    _limit = limit


  new start_at(starting_value: A = A.from[USize](1), limit: A = A.max_value()) =>
    """
    Create a new prime iterator starting from a particular value `starting_value`,
    allowing for calculation of primes up the next prime after `limit`.

    ```pony
      let iterator = PrimeIterator.start_at(10000)
      env.out.print(iterator.next().string()) // 10007
    ```
    """
    _last = if (starting_value % A.from[USize](2)) == A.from[USize](0) then
        starting_value - A.from[USize](1) // it isn't an issue if it already is 0
      else
        starting_value
      end
    _limit = limit


  fun has_next(): Bool =>
    """
    Do we need to find more prime number? Return `true` if the iterator hasn't
    reached the upper limit.
    
    As prime numbers are inifinite, there's always a next prime... Bu we don't
    know if this prime is before or after the upper limit of the iterator.
    `has_next` returns `true` when the last prime number found is over the
    upper limit.
    """
    _last < _limit


  fun ref next(): A =>
    """
    Find a new prime, up from the last one found.
    
    Finding the next prime requires testing for primality all the odd numbers
    upper from the last one found and this operation is calculation intensive
    for big numbers.
    """
    if _last <= A.from[USize](1) then
      _last = A.from[USize](2)
      return _last
    elseif _last == A.from[USize](2) then
      _last = A.from[USize](3)
      return _last
    end
    _last = _last + A.from[USize](2)
    while not Prime[A].is_prime(_last) do
      _last = _last + A.from[USize](2)
    end
    _last


class PrimeSieve
  """
  [Sieve of Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes):
  precomputes all primes up to a given `limit` in O(n log log n) time and O(n)
  bits of space, then answers primality queries in O(1).

  More efficient than `PrimeIterator` when many primes in a known range are
  needed, at the cost of upfront memory allocation. The sieve uses a `BitMap`
  data structure.

  **Usage:**
  ```pony
    let sieve = PrimeSieve(1_000_000)
    try
      if sieve.is_prime(999983)? then
        env.out.print("prime!")
      end
    end
    let ps = sieve.all_primes()  // all primes up to 1_000_000
    env.out.print(sieve.count().string() + " primes")
  ```

  If you need finding the primes on an interval [low, high], look at
  [`SegmentedSieve`](#SegmentedSieve).
  """

  let _sieve: BitMap
    """
    The bitmap sieve that contains the pre-calculated primes.
    """

  let _limit: USize
    """
    The upper limit of the sieve.
    """


  new create(n: USize) =>
    """
    Construct a sieve covering all integers in `[0, n]`.
    """
    _limit = n
    _sieve = BitMap(0, n + 1)

    if n >= 2 then
      // All numbers >= 2 start as prime candidates; 0 and 1 stay false
      _sieve.set_range_in_place(2, n + 1)
      // Sieve: for each prime p, mark all multiples starting from p*p
      let sqr = n.f64().sqrt().usize() + 1
      var i: USize = 2
      while i <= sqr do
        if _sieve(i) then
          var j: USize = i * i
          while j <= n do
            _sieve.unset(j)
            j = j + i
          end
        end
        i = i + 1
      end
    end


  fun is_prime(n: USize): Bool ? =>
    """
    Returns `true` iff `n` is prime. O(1) lookup.
    Returns `false` for `n < 2` (0 and 1 are definitively not prime).
    Raises an error if `n > limit`: the sieve has no information about numbers
    beyond its range, so silently returning `false` would be misleading.
    """
    if n > _limit then error end
    if n < 2 then return false end
    _sieve(n)


  fun primes(): Iterator[USize]^ =>
    """
    Returns an iterator over all primes up to (and including) the limit,
    in ascending order. O(1) space, lazy, single-pass.
    """
    _sieve.keys()


  fun all_primes(): Array[USize] ref =>
    """
    Returns all primes up to (and including) the limit, in ascending order,
    collected into an array.

    According to the [prime-counting function](https://en.wikipedia.org/wiki/Prime-counting_function),
    the size of the resulting array is approximately `π(n) ~ n / n.log()`. Add this object size
    to the size of the sieve.
    """
    let result = Array[USize](_sieve.cardinality())
    for p in _sieve.keys() do
      result.push(p)
    end
    result


  fun count(): USize =>
    """
    Returns π(limit): the number of primes ≤ limit.
    """
    _sieve.cardinality()


  fun limit(): USize =>
    """
    Returns the upper limit of the sieve.
    """
    _limit


class SegmentedSieve
  """
  [Segmented Sieve of Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes#Segmented_sieve):
  finds all primes in the closed interval `[low, high]` using memory
  proportional to `O(sqrt(high) + (high − low))` bits — much less than a
  full sieve up to `high`.

  A helper `PrimeSieve` covers `[0, sqrt(high)]` to provide the small
  sieving primes. The segment `BitMap` is created for the range `[low, high]`
  only, indexed by absolute position, so `is_prime`, `primes`, and
  `all_primes` all work with the original numbers without any offset
  arithmetic on the caller's side.

  **Usage:**
  ```pony
    let sieve = SegmentedSieve(1_000_000_000, 1_000_001_000)
    for p in sieve.primes() do
      env.out.print(p.string())
    end
    env.out.print(sieve.count().string() + " primes in range")
  ```

  When the `low` bound of the range is small (near to 0 compared with
  `high` value), consider using [`PrimeSieve`](#PrimeSieve) instead.
  """

  let _sieve: BitMap
    """
    Segment bitmap covering `[_low, _high]` (indexed by absolute position).
    """

  let _low: USize
    """
    Lower bound of the interval (inclusive).
    """

  let _high: USize
    """
    Upper bound of the interval (inclusive).
    """


  new create(low: USize, high: USize) =>
    """
    Construct a segmented sieve covering the closed interval `[low, high]`.

    Memory used: `O(sqrt(high))` bits for the helper sieve plus
    `O(high − low)` bits for the segment bitmap.
    """
    _low = low
    _high = high
    _sieve = BitMap(low, high + 1)

    if high >= 2 then
      // Mark all numbers in [max(2, low), high] as prime candidates.
      let start = if low < 2 then 2 else low end
      _sieve.set_range_in_place(start, high + 1)

      // Build a small sieve to find all primes up to sqrt(high).
      let sqr = high.f64().sqrt().usize()
      let small = PrimeSieve(sqr)

      // For each small prime p, strike its composites in [low, high].
      for p in small.primes() do
        // Smallest multiple of p that is >= low.
        let rem = low % p
        let first_mult = if rem == 0 then low else low + (p - rem) end
        // Advance past p itself: the first *composite* multiple is >= 2*p.
        var j = if first_mult <= p then p + p else first_mult end
        while j <= high do
          _sieve.unset(j)
          j = j + p
        end
      end
    end


  fun is_prime(n: USize): Bool ? =>
    """
    Returns `true` iff `n` is prime. O(1) lookup.
    Returns `false` for `n < 2`.
    Raises an error if `n` is outside `[low, high]`: the sieve has no
    information about numbers beyond its range, so silently returning `false`
    would be misleading.
    """
    if (n < _low) or (n > _high) then error end
    if n < 2 then return false end
    _sieve(n)


  fun primes(): Iterator[USize]^ =>
    """
    Returns an iterator over all primes in `[low, high]`, in ascending order.
    O(1) space, lazy, single-pass.
    """
    _sieve.keys()


  fun all_primes(): Array[USize] ref =>
    """
    Returns all primes in `[low, high]`, in ascending order, collected into
    an array. According to the
    [prime-counting function](https://en.wikipedia.org/wiki/Prime-counting_function),
    the approximate count is `π(high) − π(low) ~ (high − low) / ln(high)`.
    """
    let result = Array[USize](_sieve.cardinality())
    for p in _sieve.keys() do
      result.push(p)
    end
    result


  fun count(): USize =>
    """
    Returns the number of primes in `[low, high]`.
    """
    _sieve.cardinality()


  fun lower(): USize =>
    """
    Returns the lower bound of the sieve interval (inclusive).
    """
    _low


  fun limit(): USize =>
    """
    Returns the upper bound of the sieve interval (inclusive).
    """
    _high
