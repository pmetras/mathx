
use "collections"

use "../assertx"


primitive NTT[N: UnsignedInteger[N] val = USize]

  fun naive_transform(invec: Array[N], root: N, mod: N, inverse: Bool = false): Array[N]^ =>
    """
    Returns the number-theoretic transform of the given vector `invec` with
    respect to the given primitive nth root of unity `root` under the given
    modulus `mod`.
    """
    let size = invec.size()
    let outvec = Array[N](size)
    
    try
      if not inverse then
        for i in Range[USize](0, size) do
          var sum: N = N.from[USize](0)
          for j in Range[USize](0, size) do
            let k = Modular[USize].mul_mod(i, j, size)
            let temp = Modular[N].mul_mod(invec(j)?, Modular[N].pow_mod(root, N.from[USize](k), mod), mod)
            sum = Modular[N].add_mod(sum, temp, mod)
          end
          outvec.push(sum)
        end
      else
        let scaler = Modular[N].inverse_mod(N.from[USize](size), mod)
        for i in Range[USize](0, size) do
          var sum: N = N.from[USize](0)
          for j in Range[USize](0, size) do
            let k = Modular[USize].mul_mod(i, j, size)
            let temp = Modular[N].mul_mod(invec(j)?, Modular[N].pow_mod(Modular[N].inverse_mod(root, mod), N.from[USize](k), mod), mod)
            sum = Modular[N].add_mod(sum, temp, mod)
          end
          outvec.push(Modular[N].mul_mod(sum, scaler, mod))
        end
      end
    else
      Fail("Index error when accessing array")
    end
    outvec
    
  
  fun transform(a: Array[N], root: N, mod: N, inverse: Bool = false): Array[N] =>
    """
    Computes the number-theoritic transform of vector `a` in-place, with respect
    to the given primitive nth root of unity `root` under modulus `mod`. It
    calculates the inverse transform when `inverse` is `true`.
    
    Size of `a` must be a power of 2.
    
    It uses in-place DFT algorithm that requires that `a` size must be a power
    of 2. It builds a table of `a.size()` of powers of `root` values that are
    used for the calculation. It runs in `O(log2(a.size()))` for execution.
    
    """
    let size = a.size()
    try
      Assert((size and (size - 1)) == 0, "The input array 'a' size (" +
             size.string() + ") must be a power of 2. Enlarge it to " +
             size.next_pow2().string())?
      Assert(size > 0, "Array 'a' can't be empty")?
    end

    // Prepare root power tables using recurrence
    let size2 = size / 2
    let root1 = Array[N](size2)
    var temp: N = N.from[USize](1)
    for i in Range[USize](0, size2) do
      root1.push(temp)
      temp = Modular[N].mul_mod(temp, root, mod)
    end

    // Reorganize the array content in bits order
    _rearrange[N](a)

    // Cooley-Tukey decimation-in-time radix-2 FFT
    try
      var step: USize = 2
      while step <= size do
        let half = step / 2
        let jump = size / step
        for group in Range[USize](0, size, step) do
          var k: USize = 0
          for i in Range[USize](group, group + half) do
            let j = i + half
            let left = a(i)?
            let right = Modular[N].mul_mod(a(j)?, root1(k)?, mod)
            a.update(i, Modular[N].add_mod(left, right, mod))?
            a.update(j, Modular[N].sub_mod(left, right, mod))?
            k = k + jump
          end
        end

        if step == size then
          break
        else
          step = 2 * step
        end
      end
    else
      Fail("Index out of bounds [0.." + size.string() + ")")
    end

    // If inverse, scale the result if normalize is required
    if inverse and normalize then
      _normalize(a)
    end

    // Return the changed array now containing FFT
    a







  fun _rearrange[T: Any val](a: Array[T]): Array[T] =>
    """
    Rearrange the indexes of array `a` in bits order.

    The indexes of `a`, when printed in binary, are reversed and exchanged.
    For instance, the element of `a` at index 4 is `100` in binary. When
    reversed, it becomes `001`, and the element is switched to index 1.
    
    This function also exists in `FFT` class.
    """
    let size = a.size()
    let width = size.bitwidth()
    let clz1 = 1 + size.clz()
    var j: USize = 0
    try
      for i in Range[USize](0, size) do
        j = i.bit_reverse() >> clz1
        if j > i then
          a.update(i, a.update(j, a(i)?)?)? // Swap a(i) and a(j) values
        end
      end
    else
      Fail("Index j (" + j.string() + ") out of bounds [0.." + size.string() + ")")
    end
    a


  fun _normalize(a: Array[Complex[F]]): Array[Complex[F]]^ =>
    """
    Normalize the value of items in array `a`, i.e. divide by size of `a`.
    """
    let size = a.size()
    let scale = Complex[F](F.from[USize](size))
    try
      for i in Range(0, size) do
        a.update(i, a(i)? / scale)?
      end
    else
      Fail("Index out of bound [0.." + size.string() + ") in normalization")
    end
    a


