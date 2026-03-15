//- NTT/INTT - Number Theoretic Transform and Inverse

use "collections"

use "../assertx"



primitive NTT[A: UnsignedInteger[A] val = USize]
  """
  The number-theoretic transform (NTT) is obtained by specializing the discrete
  Fourier transform to Z/p, the integers modulo a prime p.

  See https://en.wikipedia.org/wiki/Discrete_Fourier_transform_over_a_ring#Number-theoretic_transform
  """

  fun _p(): A =>
    """
    Modulo: NTT-friendly prime number of the form $k * 2^{n} + 1$.
    To be able to use an array of size $N$, we must have $2^n \ge N$.

    $998244353 - 1 = 119 \cdot 2^{23}$ (32 bits)

    $18446744069414584320 = 2^{64} - 2^{32} + 1 - 1 = 2^{32} \cdot 3 \cdot 5 \cdot 17 \cdot 257 \cdot 65537$ (64 bits)

    For 32 bits architecture, it means that this number is able to support NTT
    transforms with a size up to $2^{23}$ items (~8 millions elements in array).
    """
    ifdef ilp32 or lp64 then
      A.from[USize](998244353)
    elseif llp64 then
      A.from[USize](18446744069414584320)
    else // Be safe and consider as 32 bits
      A.from[USize](998244353)
    end


  fun _g(): A =>
    """
    Primitive root modulo `_p`: for any prime factor $q$ of `_p - 1`, we
    must have $$g^{(p - 1)/q} \not\equiv 1 \pmod p$$.
    """
    ifdef ilp32 or lp64 then
      A.from[USize](3)
    elseif llp64 then
      A.from[USize](7)
    else // Be safe and consider as 32 bits
      A.from[USize](7)
    end


  fun _precompute_roots(n: USize, inverse: Bool): Array[A] =>
    """
    Precompute the powers of omega $\omega^N \equiv 1 \pmod p$.
    """
    let one = A.from[USize](1)
    let n2 = A.from[USize](n / 2)
    let roots = Array[A](n / 2)
    let g_root = if inverse then Modular[A].inv_mod(_g(), _p()) else _g() end
    let w_n = Modular[A].pow_mod(g_root, (_p() - one) / A.from[USize](n), _p())
    
    var current_w = one
    for i in Range(0, n / 2) do
      roots.push(current_w)
      current_w = Modular[A].mul_mod(current_w, w_n, _p())
    end
    roots


  fun _rearrange(a: Array[A], n: USize) =>
    """
    Rearrange the elements in array `a` in bits order.

    The indexes of `a`, when printed in binary, are reversed and exchanged.
    For instance, the element of `a` at index 4 is `100` in binary. When
    reversed, it becomes `001`, and the element is switched to index 1.
    
    This function also exists in `FFT` class.
    """
    var j: USize = 0
    for i in Range(0, n) do
      if i < j then
        try
          // Swap values
          a(i)? = (a(j)? = a(i)?)
        end
      end
      var m = n >> 1
      while (m >= 1) and (j >= m) do
        j = j - m
        m = m >> 1
      end
      j = j + m
    end


  fun transform(a: Array[A], inverse: Bool) =>
    """
    Number-Theoretic Transform. Like a Fourier Transform, but on Z/p...
    Computes the number-theoretic transform of vector `a` in-place. It
    calculates the inverse transform when `inverse` is `true`.
    
    Size of `a` must be a power of 2.
    
    It uses in-place DFT algorithm that requires that `a` size must be a power
    of 2. It builds a table of `a.size()` of powers of `root` values that are
    used for the calculation. It runs in `O(log2(a.size()))` for execution.
    """
    let n = a.size()
    ifdef debug then
      try
        Assert((n and (n - 1)) == 0, "[NTT.transform] The input array 'a' size (" +
              n.string() + ") must be a power of 2. Enlarge it to " +
              n.next_pow2().string(), true)?
        Assert(n > 0, "[NTT.transform] Array 'a' can't be empty", true)?

        let zero = A.from[USize](0)
        for i in Range(0, n) do
          Assert(try a(i)? >= zero else false end, "[NTT.transform] Elements of array a(" +
                i.string() + ") = " + (try a(i)?.string() else "?" end) + " must be positive", true)?
        end
      end
    end


    _rearrange(a, n)
    
    // 1. Precompute roots of 1 on Z/n
    let roots = _precompute_roots(n, inverse)

    // 2. Butterfly loop
    var len: USize = 2
    while len <= n do
      let step = n / len
      var i: USize = 0
      while i < n do
        for j in Range(0, len / 2) do
          try
            // Use the pre-calculated Twiddle Factor
            let w = roots(j * step)? 
            let u = a(i + j)?
            let v = Modular[A].mul_mod(a(i + j + (len / 2))?, w, _p())
            
            a(i + j)? = Modular[A].add_mod(u, v, _p())
            a(i + j + (len / 2))? = Modular[A].sub_mod(u, v, _p())
          end
        end
        i = i + len
      end
      len = len << 1
    end

    if inverse then
      let n_inv = Modular[A].inv_mod(A.from[USize](n), _p())
      for k in Range(0, n) do
        try a(k)? = Modular[A].mul_mod(a(k)?, n_inv, _p()) end
      end
    end


  fun naive_transform(a: Array[A], inverse: Bool = false, root: A, mod: A): Array[A]^ =>
    """
    Returns the number-theoretic transform of the given vector `a`, or the
    inverse transform when `inverse` is set, with respect to the given
    primitive nth root of unity `root` under the given modulus `mod`.

    This naive implementation is not optimized but accepts an array whose size is
    not a power of 2.

    For good values, use `root = A.from[USize](3)` and `mod = A.from[USize](998244353)`.
    """
    let size = a.size()
    ifdef debug then
      try
        Assert(size > 0, "[NTT.naive_transform] Array 'a' can't be empty", true)?
        let zero = A.from[USize](0)
        for i in Range(0, size) do
          Assert(try a(i)? >= zero else false end, "[NTT.naive_transform] Elements of array a(" +
                i.string() + ") = " + (try a(i)?.string() else "?" end) + " must be positive", true)?
        end
      end
    end

    let result = Array[A](size)
    
    try
      if not inverse then
        for i in Range(0, size) do
          var sum: A = A.from[USize](0)
          for j in Range(0, size) do
            let k = Modular[USize].mul_mod(i, j, size)
            let temp = Modular[A].mul_mod(a(j)?, Modular[A].pow_mod(root, A.from[USize](k), mod), mod)
            sum = Modular[A].add_mod(sum, temp, mod)
          end
          result.push(sum)
        end
      else
        let scaler = Modular[A].inv_mod(A.from[USize](size), mod)
        for i in Range(0, size) do
          var sum: A = A.from[USize](0)
          for j in Range[USize](0, size) do
            let k = Modular[USize].mul_mod(i, j, size)
            let temp = Modular[A].mul_mod(a(j)?, Modular[A].pow_mod(Modular[A].inv_mod(root, mod), A.from[USize](k), mod), mod)
            sum = Modular[A].add_mod(sum, temp, mod)
          end
          result.push(Modular[A].mul_mod(sum, scaler, mod))
        end
      end
    else
      Fail("[NTT.naive_transform] Index error when accessing array")
    end
    result
