//- NTT/INTT - Number Theoretic Transform and Inverse

use "collections"

use "../assertx"
use "../formatx"


primitive NTT[A: (UnsignedInteger[A] & Any val) = USize]
  """
  The number-theoretic transform (NTT) is obtained by specializing the discrete
  Fourier transform to Z/p, the integers modulo a prime p.

  See https://en.wikipedia.org/wiki/Discrete_Fourier_transform_over_a_ring#Number-theoretic_transform
  """

  fun _p(): A =>
    """
    Modulo: NTT-friendly prime number of the form $k * 2^{n} + 1$.
    To be able to use an array of size $N$, we must have $2^n ≥ N$.

    998244353 - 1 = 119 ⋅ 2^{23} (32-bit prime, used for U32 / USize on ilp32)

    2^{64} - 2^{32} + 1 = 18446744069414584321, where
    p - 1 = 2^{32} ⋅ 3 ⋅ 5 ⋅ 17 ⋅ 257 ⋅ 65537 (64-bit prime, used for U64 /
    USize on lp64 / llp64, and all generic types such as MPInt)
    """
    let a_any: Any val = A.from[U8](0)
    match a_any
    | let _: U32 =>
      // p - 1 = 2^23 * 119
      A.from[U32](998244353)
    | let _: U64 =>
      // 2^64 - 2^32 + 1; p - 1 = 2^32 * 3 * 5 * 17 * 257 * 65537
      A.from[U64](18446744069414584321)
    | let _: USize =>
      ifdef ilp32 then
        // p - 1 = 2^23 * 119
        A.from[U32](998244353)
      else // lp64 or llp64 — both have 64-bit registers; use 64-bit prime
        // 2^64 - 2^32 + 1; p - 1 = 2^32 * 3 * 5 * 17 * 257 * 65537
        A.from[U64](18446744069414584321)
      end
    else
      // Generic type (e.g., MPInt): use the 64-bit NTT-friendly prime.
      // The caller's type must be capable of representing 18446744069414584321.
      // 2^64 - 2^32 + 1; p - 1 = 2^32 * 3 * 5 * 17 * 257 * 65537
      ifdef debug then
        (A.from[U64](18446744069414584321) == A.from[U128](18446744069414584321)) or
          Fail("[Modular[A]._p] Type `A` cannot represent 2^64 - 2^32 + 1. This type is not suitable for NTT.")
      end
      A.from[U64](18446744069414584321)
    end


  fun _g(): A =>
    """
    Primitive root modulo `_p`: for any prime factor $q$ of `_p - 1`, we
    must have $g^{(p - 1)/q} ≡ 1 mod p$.

    - 32-bit prime (998244353): primitive root = 3
    - 64-bit prime (2^64−2^32+1): primitive root = 7
    """
    let a_any: Any val = A.from[U8](0)
    match a_any
    | let _: U32 =>
      A.from[U32](3)
    | let _: U64 =>
      A.from[U32](7)
    | let _: USize =>
      ifdef ilp32 then
        A.from[U32](3)
      else // lp64 or llp64
        A.from[U32](7)
      end
    else
      // Generic type: primitive root for the 64-bit NTT prime (2^64−2^32+1) is 7.
      A.from[U32](7)
    end


  fun _mul_mod_p(a: A, b: A): A =>
    """
    Multiply `a * b mod _p()` using fast prime-specific reduction.

    For 32-bit types (U32 or USize on ilp32), the prime is p = 998244353 < 2^30.
    Both operands fit in 30 bits, so the U64 product fits in 60 bits and
    a single 64-bit division gives the exact remainder.

    For 64-bit types (U64 or USize on lp64/llp64), the prime is p = 2^64 − 2^32 + 1.
    Since 2^64 ≡ 2^32 − 1 (mod p), a 128-bit product is reduced in two
    steps, each collapsing the high 64 bits via that identity, followed by
    a single conditional subtraction. This replaces a 64-iteration
    binary-doubling loop in `Modular[A].mul_mod`.

    For generic types (e.g., MPInt), falls back to `Modular[A].mul_mod`.
    """
    let a_any: Any val = a
    match a_any
    | let _: U32 =>
      // p < 2^30; a*b < 2^60 fits in U64.
      A.from[U64]((a.u64() * b.u64()) % _p().u64())
    | let _: U64 =>
      _mul_mod_p_u64(a.u64(), b.u64())
    | let _: USize =>
      ifdef ilp32 then
        // p < 2^30; a*b < 2^60 fits in U64.
        A.from[U64]((a.u64() * b.u64()) % _p().u64())
      else
        _mul_mod_p_u64(a.u64(), b.u64())
      end
    else
      Modular[A].mul_mod(a, b, _p())
    end


  fun _mul_mod_p_u64(a: U64, b: U64): A =>
    """
    Reduction for p = 2^64 − 2^32 + 1 via two-step 128-bit collapse.
    """
    // Step 1: 128-bit product.
    let x: U128 = a.u128() * b.u128()
    let x_hi: U64 = (x >> 64).u64()  // < p < 2^64
    let x_lo: U64 = x.u64()

    // Step 2: first reduction. 2^64 ≡ 2^32 − 1 (mod p), so
    //   x ≡ x_hi*(2^32 − 1) + x_lo  (mod p)
    //     = (x_hi << 32) − x_hi + x_lo
    // x_hi < 2^64, so x_hi.u128() << 32 < 2^96 — fits in U128.
    let y: U128 = ((x_hi.u128() << 32) - x_hi.u128()) + x_lo.u128()

    // Step 3: second reduction. y < 2^96, so y >> 64 < 2^32.
    let y_hi: U64 = (y >> 64).u64()   // < 2^32
    let y_lo: U64 = y.u64()
    // y_hi*(2^32 − 1) = (y_hi << 32) − y_hi, fits in U64 since y_hi < 2^32.
    let z_part: U64 = (y_hi << 32) - y_hi
    (let z_sum, let z_carry) = z_part.addc(y_lo)

    // Step 4: z < 2p — one conditional subtraction.
    //   If z overflowed U64 (carry=1): z = 2^64 + z_sum ≥ 2^64 > p.
    //     result = z − p = z_sum + (2^64 − p) = z_sum + (2^32 − 1).
    //   Otherwise: result = z_sum − p  if z_sum ≥ p, else z_sum.
    let p: U64 = _p().u64()
    let result: U64 =
      if z_carry then
        z_sum + U32.max_value().u64()  // + (2^32 − 1) = 2^64 − p
      elseif z_sum >= p then
        z_sum - p
      else
        z_sum
      end
    A.from[U64](result)


  fun _add_mod_p(a: A, b: A): A =>
    """
    Compute `(a + b) mod _p()` for `a, b` already in `[0, p)`.

    For 32-bit types (U32 or USize on ilp32), widens to U64 to avoid overflow.
    For 64-bit types (U64 or USize on lp64/llp64), handles U64 carry with a
    single conditional add of `2^64 − p = 2^32 − 1` instead of a subtraction.
    For generic types (e.g., MPInt), falls back to `Modular[A].add_mod`.
    """
    let a_any: Any val = a
    match a_any
    | let _: U32 =>
      // p < 2^30; a + b < 2p < 2^31 — no overflow in U64.
      let sum: U64 = a.u64() + b.u64()
      let p: U64 = _p().u64()
      A.from[U64](if sum >= p then sum - p else sum end)
    | let _: U64 =>
      _add_mod_p_u64(a.u64(), b.u64())
    | let _: USize =>
      ifdef ilp32 then
        // p < 2^30; a + b < 2p < 2^31 — no overflow in U64.
        let sum: U64 = a.u64() + b.u64()
        let p: U64 = _p().u64()
        A.from[U64](if sum >= p then sum - p else sum end)
      else
        _add_mod_p_u64(a.u64(), b.u64())
      end
    else
      Modular[A].add_mod(a, b, _p())
    end


  fun _add_mod_p_u64(a: U64, b: U64): A =>
    """
    Add reduction for p = 2^64 − 2^32 + 1.
    carry=1 → 2^64 + sum ≥ 2^64 > p; result = sum + (2^64 − p) = sum + (2^32 − 1).
    """
    (let sum, let carry) = a.addc(b)
    let p: U64 = _p().u64()
    let result: U64 =
      if carry then
        sum + U32.max_value().u64()
      elseif sum >= p then
        sum - p
      else
        sum
      end
    A.from[U64](result)


  fun _sub_mod_p(a: A, b: A): A =>
    """
    Compute `(a - b) mod _p()` for `a, b` already in `[0, p)`.

    For 32-bit types (U32 or USize on ilp32), uses a simple comparison.
    For 64-bit types (U64 or USize on lp64/llp64), uses U64 `subc` with
    borrow correction.
    For generic types (e.g., MPInt), falls back to `Modular[A].sub_mod`.
    """
    let a_any: Any val = a
    match a_any
    | let _: U32 =>
      let au: U64 = a.u64(); let bu: U64 = b.u64()
      let p: U64 = _p().u64()
      A.from[U64](if au >= bu then au - bu else (au + p) - bu end)
    | let _: U64 =>
      _sub_mod_p_u64(a.u64(), b.u64())
    | let _: USize =>
      ifdef ilp32 then
        let au: U64 = a.u64(); let bu: U64 = b.u64()
        let p: U64 = _p().u64()
        A.from[U64](if au >= bu then au - bu else (au + p) - bu end)
      else
        _sub_mod_p_u64(a.u64(), b.u64())
      end
    else
      Modular[A].sub_mod(a, b, _p())
    end


  fun _sub_mod_p_u64(a: U64, b: U64): A =>
    """
    Sub reduction for p = 2^64 − 2^32 + 1.
    borrow=1 → a < b in U64; diff = a − b + 2^64.
    result = a − b + p = diff − (2^64 − p) = diff − (2^32 − 1).
    """
    (let diff, let borrow) = a.subc(b)
    let result: U64 = if borrow then
      diff - U32.max_value().u64()
    else
      diff
    end
    A.from[U64](result)


  fun _pow_mod_p(base: A, exp: A): A =>
    """
    Compute `base ^ exp mod _p()` using binary exponentiation with `_mul_mod_p`,
    avoiding the generic `Modular[A].mul_mod` slow path.
    """
    let zero = A.from[USize](0)
    let one  = A.from[USize](1)
    let two  = A.from[USize](2)
    if exp == zero then
      return one
    end

    var b = base % _p()
    var res = one
    var e = exp
    while e > zero do
      if (e % two) == one then
        res = _mul_mod_p(res, b)
      end
      e = e / two
      b = _mul_mod_p(b, b)
    end
    res


  fun _precompute_roots(n: USize, inverse: Bool): Array[A] =>
    """
    Precompute the powers of omega $ω^N ≡ 1 mod p$.
    """
    let one = A.from[USize](1)
    let roots = Array[A](n / 2)
    let g_root = if inverse then Modular[A].inv_mod(_g(), _p()) else _g() end
    let w_n = _pow_mod_p(g_root, (_p() - one) / A.from[USize](n))

    var current_w = one
    for i in Range(0, n / 2) do
      roots.push(current_w)
      current_w = _mul_mod_p(current_w, w_n)
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
    // Pre-conditions
    ifdef debug then
      ((n and (n - 1)) == 0) or
        Fail(Format("[NTT.transform] The input array 'a' size ({}) must be a" +
                    " power of 2. Enlarge it to {}", [n; n.next_pow2()]))
      (n > 0) or Fail("[NTT.transform] Array 'a' can't be empty")

      let zero = A.from[USize](0)
      try
        for i in Range(0, n) do
          (a(i)? >= zero) or
            Fail(Format("[NTT.transform] Array element must be non-negative." +
                        " a({}) = {}", [i; a(i)?]))
        end
      end
    end

    _rearrange(a, n)

    // 1. Precompute roots of unity in Z/p.
    let roots = _precompute_roots(n, inverse)

    // 2. Butterfly loop.
    var len: USize = 2
    while len <= n do
      let step = n / len
      var i: USize = 0
      while i < n do
        let half = len / 2
        for j in Range(0, half) do
          try
            let w = roots(j * step)?
            let u = a(i + j)?
            let v = _mul_mod_p(a(i + j + half)?, w)
            a(i + j)? = _add_mod_p(u, v)
            a(i + j + half)? = _sub_mod_p(u, v)
          end
        end
        i = i + len
      end
      len = len << 1
    end

    // 3. Inverse: divide all elements by n (mod p).
    if inverse then
      let n_inv = Modular[A].inv_mod(A.from[USize](n), _p())
      for k in Range(0, n) do
        try a(k)? = _mul_mod_p(a(k)?, n_inv) end
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
    // Pre-conditions
    let zero = A.from[USize](0)
    ifdef debug then
      (size > 0) or Fail("[NTT.naive_transform] Array 'a' can't be empty")
      try
        for i in Range(0, size) do
          (a(i)? >= zero) or
            Fail(Format("[NTT.naive_transform] Array element must be non-negative." +
                        " a({}) = {}", [i; a(i)?]))
        end
      end
    end

    let result = Array[A](size)

    try
      if not inverse then
        for i in Range(0, size) do
          var sum: A = zero
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
          var sum: A = zero
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
