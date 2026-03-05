// Modular arithmetic

use "../assertx"

primitive Modular[A: UnsignedInteger[A] val = USize]
  """
  # Modular arithmetic

  Operations on modular arithmetic: unsigned integer parameters are considered
  modulo `m` (elements of Z/nZ).
  """

  fun add_mod(a: A, b: A, m: A): A =>
    """
    Calculate `(a + b) % m` i.e. `a` plus `b` modulo `m`, taking care of
    overflows.
    """
    let a' = a % m
    let b' = m - (b % m)
    if b' == m then
      a'
    elseif a' >= b' then
      a' - b'
    else
      (m - b') + a'
    end
    

  fun sub_mod(a: A, b: A, m: A): A =>
    """
    Calculate `(a - b) % m`, i.e. `a` minus `b` modulo `m`, taking care of
    overflows.
    """
    let a' = a % m
    let b' = b % m
    if a' >= b' then
      a' - b'
    else
      (m - b') + a'
    end


  fun neg_mod(a: A, m: A): A =>
    """
    Negation of `a` modulo `m`.
    """
    let zero = A.from[USize](0)
    let a' = a % m
    if a' == zero then
      zero
    else
      m - a'
    end


  fun mul_mod(a: A, b: A, m: A): A =>
    """
    Calculate `(a * b) % m`, i.e. `a` times `b` modulo `m`, taking care of
    overflows.

    Without overflows, `mul_mod(a, b, m) = ((a % m) * (b % m)) % m`.
    """
    let zero = A.from[USize](0)

    // If a and b are small enough, we won't overflow and we can use direct
    // multiplication
    if (a.clz() + b.clz()) > zero.bitwidth() then
      ((a % m) * (b % m)) % m
    else
      let one = A.from[USize](1)
      let two = A.from[USize](2)

      var res = zero
      var a' = a % m
      var b' = b
      while (b' > zero) do
        if ((b' % two) == one) then
          res = (res + a') % m
        end

        a' = (a' * two) % m
        b' = b' / two
      end

      res % m
    end


  fun inverse_mod(a: A, m: A): A =>
    """
    Calculate the inverse of `a` modulo `m`. This is possible only when
    `a` is coprime with `m`, that is `gcd(a, m) == 1`.

    When `a` has no inverse, returns the value `0`.

    See [explanations](https://en.wikipedia.org/wiki/Modular_multiplicative_inverse)
    and Knuth TAOCP 2.4.5.2 Algorithm X, Extended Euclid's algorithm.
    """
    // Initialize
    let zero = A.from[USize](0)
    let one = A.from[USize](1)
    var positive: Bool = true
    var u = one
    var u' = a
    var v = zero
    var v' = m
    
    // Loop while v' != 0
    while v' != zero do
      // Divide and "substract"
      let quotient = u' / v'
      // Swap
      (u, v, u', v') = (v, u + (quotient * v), v', u' % v')
      positive = not positive
    end
    
    // Make sure u' = gcd(a, m) == 1
    if u' != one then
      zero
    // Ensure a positive result
    elseif not positive then
      m - u
    else
      u
    end


  fun div_mod(a: A, b: A, m: A): A =>
    """
    Division of `a` by `b` modulo `m`. It is defined as the multiplication of
    `a` by `1/b`, the inverse of `b` modulo `m`. But not all values of `b`
    have an inverse value, only those where `b` is coprime with `m`, i.e.
    `gcd(b, m) == 1`.

    When `b` has no inverse, returns the value `0`.
    """
    let zero = A.from[USize](0)
    let inv = inverse_mod(b, m)
    if inv == zero then
      // Division is not possible
      zero
    else
      mul_mod(a, inv, m)
    end


  fun pow_mod(a: A, b: A, m: A): A =>
    """
    Calculate `a^b % m`, i.e. `a` to the power of `b` modulo `m`, taking care
    of overflows.

    This is done by using modular multiplication `mul_mod`.
    """
    let zero = A.from[USize](0)
    let one = A.from[USize](1)

    // a^0 = 1
    if b == zero then
      return one
    end

    var a' = a % m
    var res = zero

    if a' != zero then
      let two = A.from[USize](2)

      res = one
      var b' = b
      while b' > zero do
        if (b' % two) == one then
          res = mul_mod(res, a', m)
        end

        b' = b' / two
        a' = mul_mod(a', a', m)
      end
    end
    res


  fun gcd(a: A, b: A): A =>
    """
    Returns the greatest common divisor between `a` and `b`, that is, the product
    of all the shared prime factors of `a` and `b`. For example, `gcd(3, 4)`
    would be `1`, as they do not share any prime factor, while `gcd(12, 15)`
    would be `3`, as both are multiples of `3`.

    The GCD is calculated by the [Euclidean algorithm](https://en.wikipedia.org/wiki/Euclidean_algorithm).
    
    **Usage:**
    ```pony
      env.out.print(Modular.gcd(24, 78).string()) // should be 6
    ```
    """
    (var big, var small) = if a > b then (a, b) else (b, a) end
    let zero = A.from[USize](0)
    var res: A = small
    while small != zero do
      res = small
      small = big % small
      big = res
    end
    res


  fun gcd2(a: A, b: A): A =>
    """
    Calculate `a` and `b` GCD (Greatest Common Divisor) using the [Euclidean
    binary GCD algorithm](https://en.wikipedia.org/wiki/Binary_GCD_algorithm)
    that is a variant of Euclide's algorithm. Its complexity is
    O(a.log() + b.log())^2).
    """
    let zero = A.from[USize](0)
    if a == zero then
      return b
    elseif b == zero then
      return a
    end

    let i = a.ctz()
    var u = a.shr(i)
    let j = b.ctz()
    var v = b.shr(j)
    let k = if i < j then i else j end
    
    let one = A.from[USize](1)
    let two = A.from[USize](2)

    while true do
      try
        Assert(((u % two) == one) and ((v % two) == one), "u (" + u.string() +
               ") and v (" + v.string() + ") should be even...")?
      end
    
      if u > v then
        u = v = u // Swap u and v values
      end
      v = v - u
      
      if v == zero then
        break
      end
      
      v = v >> v.ctz()
    end
    u << k


  fun lcm(a: A, b: A): A =>
    """
    Returns the least common multiplier. This is the operation one would use
    when doing an addition of two fractions. This function returns the smallest
    number which is a multiplier of both `a` and `b`.

    For instance, `lcm(3, 4)` would be `12`, while `lcm(6, 8)` would be `24`,
    as `6*4 = 24` and `8*3 = 24`.

    **Usage:**
    ```pony
      env.out.print(Modular.lcm(24, 78).string()) // 312
    ```
    """
    (a * b) / gcd2(a, b)



