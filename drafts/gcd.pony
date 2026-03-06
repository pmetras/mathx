
// TODO: Write tests

primitive GCD[U: ((U8 | U16 | U32 | U64 | U128 | ULong) & UnsignedInteger[U])]
  """
  From [Wikipedia](https://en.wikipedia.org/wiki/Greatest_common_divisor), the
  Greatest Common Divisor (GCD) of two or more integers, which are not all 0, is
  the largest positive integer that divides each of the integers.
  """

  fun apply(a: U, b: U): U =>
    """
    Calculate `a` and `b` GCD (Greatest Common Divisor) using the [binary GCD
    algorithm](https://en.wikipedia.org/wiki/Binary_GCD_algorithm) that is a
    variant of Euclide's algorithm. Its complexity is O(a.log() + b.log())^2).

    This implementation is inspired from
    https://lemire.me/blog/2013/12/26/fastest-way-to-compute-the-greatest-common-divisor/
    """
    if a == U.from[USize](0) then
      return b
    end
    if b == U.from[USize](0) then
      return a
    end
    // Count the number of trailing zeroes
    let shift = (a or b).ctz()
    var u = a >> a.ctz()
    var v = b
    repeat
      v = v >> v.ctz()
      if u > v then
        let t = v
        v = u
        u = t
      end
      v = v - u
    until v == U.from[USize](0) end
    u << shift


  fun lcm(a: U, b: U): U =>
    """
    Calculate `a` and `b` LCM
    [(Least Common Multiple)](https://en.wikipedia.org/wiki/Least_common_multiple).
    This is done by calculating the GCD first.
    """
    if (a == U.from[USize](0)) or (b == U.from[USize](0)) then
      return U.from[USize](0)
    end
    (a / GCD[U](a, b)) * b
