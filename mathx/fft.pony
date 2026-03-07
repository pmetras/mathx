// Fast Fourier Transform

use "debug"
use "format"

use "collections"

use "../assertx"


primitive FFT[F: (Float & FloatingPoint[F]) = F64]
  """
  Discrete Fast Fourier Transform on complex numbers. A detailed explanation
  of the Cooley-Tukey algorithm used is available on the
  [Fast Fourier Transform article](https://www.librow.com/articles/article-10)
  from Sergey Chernenko.

  This primitive offers multiple implementations of the DFT (Discrete Fourier
  Transform) using the FFT (Fast Fourier Transform) algorithm:

  * `fourier`: In-place speed optimized DFT on `Complex[F]`.
  * `fourier_unsafe`: The same but faster with unsafe operations.
  * `fourier_complex`: In-place speed optimized DFT on `F` real and imaginary
    numbers.
  * `fourier_real`: In-place DFT on `F` real values.
  * `fourier2`: In-place DFT optimized for accuracy.
  * `bluestein`: In-place DFT using Bluestein's chir z-transform algorithm that
    does not require the input vector size to be a power-of-2.
  * `convolve`: A circular convolution implementation.

  The DFT definition that we use is:
  $$F_{k} = \sum_{t=0}^{n-1} f_{t} e^{-2 π i t k / n}$$

  The inverse DFT is:
  $$f_{k} = 1 / n \sum_{t=0}^{n-1} F_{t} e^{2 π i t k / n}$$

  Note that these implementations change the values in the original array, like
  most DFT implementations, and also return a result that is an alias of the
  original array so that they can be used with a functional syntax.

  Some of these algorithms were inspired from [Numerical Recipes in C++](http://numerical.recipes/).
  This source is mostly educational and suffers from
  [numerous critics](https://www.stat.uchicago.edu/~lekheng/courses/302/wnnr/nr.html).
  Check code before using it for critical calculations. Also note that Numerical
  Recipes uses inverse definitions of DTF.
  """

  fun _rearrange[T: Any val](a: Array[T]): Array[T] =>
    """
    Rearrange the indexes of array `a` in bits order.

    The indexes of `a`, when printed in binary, are reversed and exchanged.
    For instance, the element of `a` at index 4 is `100` in binary. When
    reversed, it becomes `001`, and the element is switched to index 1.
    """
    let size = a.size()
    var j: USize = 0
    try
      // Incremental bit-reversal: reuses previous j via cheap AND/XOR/OR/shift,
      // avoiding a full bit_reverse() call for each element.
      for i in Range[USize](1, size) do
        var bit = size >> 1
        while (j and bit) != 0 do
          j = j xor bit
          bit = bit >> 1
        end
        j = j or bit
        if i < j then
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


  fun fourier(a: Array[Complex[F]],
              inverse: Bool = false,
              normalize: Bool = true)
	     : Array[Complex[F]] =>
    """
    Calculate the discrete fast Fourier tranform of `a` when `inverse = false`
    and the inverse transform when `inverse = true`.

    When calculating inverse FFT, one can be interested only in relative
    values instead of the absolute values. If it is the case, set
    `normalize = false` to avoid the normalization step.

    This algorithm uses in-place transformation, and array `a` content is
    replaced. It runs in `O(log2(a.size()))` for execution. As it uses
    `Complex[F]` numbers that are not primitive, it allocates memory for
    calculation. Use `fourier_complex` for a memory optimized version.
    """
    let size = a.size()
    ifdef debug then
      try
        Assert((size and (size - 1)) == 0, "The input array 'a' size (" +
              size.string() + ") must be a power of 2. Enlarge it to " +
              size.next_pow2().string(), true)?
        Assert(size > 0, "Array 'a' can't be empty", true)?
        Assert(normalize or inverse, "`normalize` must be " +
              "set to `false` only with inverse FFT", true)?
      end
    end

    // Rearrange array indexes in reverse bits order.
    _rearrange[Complex[F]](a)

    // Some constants
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    let pi = if inverse then
        F.from[F64](F64.pi())
      else
        -F.from[F64](F64.pi())
      end

    var i: USize = 0
    var j: USize = 0
    try
      var step: USize = 1
      while step < size do
        let jump = 2 * step
        let theta = pi / F.from[USize](step)
        let temp = (theta / two).sin()
        let multiplier = Complex[F](-two * temp * temp, theta.sin())
        var factor = Complex[F](one)

        for group in Range(0, step) do
          i = group
          while i < size do
            j = i + step
            let ai = a(i)?
            let product = factor * a(j)?
            a.update(j, ai - product)?
            a.update(i, ai + product)?
            i = i + jump
          end

          // Trigonometric recurrence:
          // e^i*(phi + theta) = e^i*phi + e^i*phi * (-2 * sin^2(theta / 2) + i*sin(theta))
          // as cos(theta) = 1 - 2 * sin^2(theta / 2)
          // We use sinus instead of cosinus because theta is small and we have
          // a better accuracy when using floats (the error is smaller with
          // numbers like 0.374e-17 (i.e. sin) rather 0.99999999999999 (i.e. cos)
          factor = factor + (multiplier * factor)
        end

        step = 2 * step
      end
    else
      Fail("Index of of bounds [0.." + size.string() + ") i=" + i.string() +
           ", j=" + j.string())
    end

    // If inverse, scale the result if normalization is required
    if inverse and normalize then
      _normalize(a)
    end

    // Return the changed array now containing FFT
    a


  fun fourier_unsafe(a: Array[Complex[F]],
                     inverse: Bool = false,
                     normalize: Bool = true)
	            : Array[Complex[F]] =>
    """
    This version of the fourier transform calculation uses unsafe operations
    and silently ignores errors.

    See [FFT.fourier](#fourier) for documentation and explanations.
    """
    let size = a.size()
    ifdef debug then
      try
        Assert((size and (size - 1)) == 0, "The input array 'a' size (" +
              size.string() + ") must be a power of 2. Enlarge it to " +
              size.next_pow2().string(), true)?
        Assert(size > 0, "Array 'a' can't be empty", true)?
        Assert(normalize or inverse, "`normalize` must be " +
              "set to `false` only with inverse FFT", true)?
      end
    end

    // Rearrange array indexes in reverse bits order.
    _rearrange[Complex[F]](a)

    // Some constants
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    let fpi = if inverse then
        F.from[F64](F64.pi())
      else
        -~F.from[F64](F64.pi())
      end

    try
      var step: USize = 1
      while step < size do
        let jump = step *~ 2
        let theta = fpi /~ F.from[USize](step)
        let temp = (theta /~ two).sin()
        let multiplier = Complex[F](-~two *~ temp *~ temp, theta.sin())
        var factor = Complex[F](one)

        for group in Range(0, step) do
          for i in Range(group, size, jump) do
            let j = i +~ step
            let product = factor *~ a(j)?
            a.update(j, a(i)? -~ product)?
            a.update(i, a(i)? +~ product)?
          end

          // Trigonometric recurrence:
          factor = factor +~ (multiplier *~ factor)
        end

        step = jump
      end

      // If inverse, scale the result if normalize is required
      if inverse and normalize then
        let scale = Complex[F](F.from[USize](size))
        for k in Range(0, size) do
          a.update(k, a(k)? /~ scale)?
        end
      end
    end

    // Return the changed array now containing FFT
    a


  fun fourier_complex(a: Array[F val],
                      inverse: Bool = false,
                      normalize: Bool = true)
                     : Array[F val] =>
    """
    Calculate the discrete fast Fourier tranform of `a`. It calculates the
    inverse transform when `inverse = true`.

    The content of `a` is alternate values of real and imaginary parts of
    complex numbers: `a(2*k) = real` and `a(2*k+1) = imag`. The size of `a`
    divided by 2 must be a power of 2.

    When calculating the inverse transform, one can decide not to normalize
    the result by setting `normalize` parameter to `false`, and getting
    relative values.

    This algorithm uses in-place transformation, and array `a` content is
    replaced. It runs in `O(log2(a.size()))` for execution and `O(1)` for
    memory. Compared with [FFT.fourier](#fourier), this implementation does
    not allocate temporary memory for `Complex` objects and doesn't use
    `Complex` multiplication.
    """
    let size = a.size()
    let size2 = size / 2 // Number of complex data points
    ifdef debug then
      try
        Assert((size2 and (size2 - 1)) == 0,
              "The input array 'a' size (" + size.string() +
              ") must be a power of 4. Enlarge it to " +
              (size2.next_pow2() * 2).string(), true)?
        Assert(size > 0, "Array 'a' can't be empty", true)?
        Assert(normalize or inverse, "`normalize` must be " +
              "set to `false` only with inverse FFT", true)?
      end
    end

    // Reorganize the array content in bits order
    var j: USize = 1
    try
      for i in Range[USize](1, size, 2) do
        if j > i then
          // Swap values
          a.update(i, a.update(j, a(i)?)?)?
          a.update(i - 1, a.update(j - 1, a(i - 1)?)?)?
        end
        var m = size2
        while (m >= 2) and (j > m) do
          j = j - m
          m = m / 2
        end
        j = j + m
      end
    else
      Fail("Index out of bounds [0.." + size2.string() + ")")
    end

    // The constants
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    let pi2 = if inverse then
        two * F.from[F64](F64.pi())
      else
        -two * F.from[F64](F64.pi())
      end
    var step: USize = 2

    try
      while step < size do
        var jump = step * 2
        let theta = pi2 / F.from[USize](step)
        let temp = (theta / two).sin()
        let multiplier_r = -two * temp * temp
        let multiplier_i = theta.sin()
        var factor_r = one
        var factor_i = zero

        for group in Range(1, step, 2) do
          for i in Range(group, size, jump) do
            j = i + step
            let product_r = (factor_r * a(j - 1)?) - (factor_i * a(j)?)
            let product_i = (factor_r * a(j)?) + (factor_i * a(j - 1)?)
            let ai_re = a(i - 1)?
            let ai_im = a(i)?
            a.update(j - 1, ai_re - product_r)?
            a.update(j,     ai_im - product_i)?
            a.update(i - 1, ai_re + product_r)?
            a.update(i,     ai_im + product_i)?
          end

          // Trigonometric recurrence
          (factor_r, factor_i) = (factor_r + ((factor_r * multiplier_r) -
                                              (factor_i * multiplier_i)),
                                  factor_i + (factor_i * multiplier_r) +
                                             (factor_r * multiplier_i))
        end

        step = jump
      end
    else
      Fail("Index out of bounds [0.." + size.string() + ")")
    end

    // If inverse, scale the result if normalize is required
    if inverse and normalize then
      let scale = F.from[USize](size2)
      try
        for i in Range(0, size) do
          a.update(i, a(i)? / scale)?
        end
      else
        Fail("Index out of bound [0.." + size.string() + ") in normalization")
      end
    end

    // Return the changed array now containing FFT
    a


  fun fourier_real(a: Array[F],
                   inverse: Bool = false,
                   normalize: Bool = true)
                  : Array[F] =>
    """
    Calculate the Fourier transform on an array `a` of real numbers instead of
    complex number. Array `a` is changed to contain the transform result. It
    replaces these data by the positive frequency half of their complex Fourier
    transform. The real-valued ﬁrst and last components of the complex transform
    are returned as elements `a(0)` and `a(1)`, respectively. This is possible
    because the Fourier transform is
    [even symmetric](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#DFT_of_real_and_purely_imaginary_signals).

    Note that the Fourier transform of real data is complex data, and the result
    array contains complex numbers with real and imaginary parts alternate
    (i.e. a(2*k) = real; a(2*k+1) = imag).

    When `inverse = true`, it calculates the inverse transform of a complex
    array (where real and imaginary parts alternate) if it is the Fourier
    transform of real data. No check is done to validate that the input array
    is the complex DFT of a real array!

    The size of `a` must be a power of 2.

    This algorithm uses in-place transformation, and array `a` content is
    replaced. It runs in `O(1)` for memory and `O(log2(a.size()))` for execution.

    Contrarily to [Numerical Recipes 12.3.2](http://numerical.recipes/), and
    all functions in this module, the FFT is defined by
    $$
    F_n = \sum_{j=0}^{N-1} e^{-\frac{2 \pi i j n}{N}} f_j
    $$

    The FFT for $\{f_n\}_{n=0}^{N}$ real, we put
    $h_k = f_{2k} + i f_{2k + 1} \qquad k = 0, \dots, N/2 - 1$
    and taking the Fourier transform
    $$
    H_k &:= \sum_{j = 0}^{N/2 - 1} e^{-\frac{2 \pi i j k }{N/2}} h_k \\
        &= \sum_{j = 0}^{N/2 - 1} e^{-\frac{2 \pi i (2j) k}{N}} f_{2k} + i \sum_{j = 0}^{N/2 - 1} e^{-\frac{2 \pi i (2j) k }{N}} f_{2k+1} \\
        &:= F^e_k + i F^o_k.
    $$

    Because of the symetry of FFT, as $f$ is real,
    $F^{e/o}_n = \bar{F}^{e/o}_{N/2 - n}$. So
    $H_{N/2 - k} = F^e_{N/2 - k} + i F^o_{N/2 - k} = \bar{F}^e_k + i \bar{F}^o_k$
    and $\bar{H}_{N/2 - k} = F^e_k - i F^o_k$, with $H_k = F^e_k + i F^o_k$,
    one can solve for $F^e$ and $F^o$:
    $$
    F^e_k &= \frac{1}{2} (H_k + \bar{H}_{N/2 - k}) \\
    F^o_k &= - \frac{i}{2} (H_k - \bar{H}_{N/2 - k}).
    $$

    In the code, $\mathrm{h1r} = \Re F^e_k$, $\mathrm{h1i} = \Im F^e_k$,
    $\mathrm{h2r} = \Re F^o_k$ and $\mathrm{h2i} = \Im F^o_k$.

    From the definition
    $$
    F_n &:= \sum_{j = 0}^{N-1} e^{- \frac{2 \pi i j n}{N}} f_n \\
        &= \sum_{j = 0}^{N/2-1} e^{- \frac{2 \pi i (2j) n}{N}} f_{2j} + \sum_{j = 0}^{N/2 - 1} e^{-\frac{2 \pi i (2j+1) n}{N}} f_{2j+1} \\
        &= F^e_n + e^{-\frac{2 \pi i n}{N}} F^o_n,
    $$
    using $\theta_k = - \frac{2 \pi i k}{N}$ for lighter formulas
    $$
    \Re F_k &= \Re F^e_k + \cos( \theta_k) \Re F^o_k - \sin (\theta_k) \Im F^o_k \qquad \text{Line a.update(i1, h1r + tr)?}\\
    \Im F_k &= \Im F^e_k + \cos (\theta_k) \Im F^o_k + \sin (\theta_k) \Re F^o_k \qquad \text{Line a.update(i2, h1i + ti)?}.
    $$
    and because of
    $$
    F_{N/2 - k} &= F^e_{N/2 - k} + e^{-\frac{2 \pi i (N/2 - k)}{N}} F^o_{N/2 - k} \\
                &= \bar{F}^e_k - e^{+\frac{2 \pi i k}{N}} \bar{F}^o_k
    $$
    we get
    $$
    \Re F_{N/2 - k} &= \Re \bar{F}^e_k - \cos(-\theta_k) \Re \bar{F}^o_k + \sin(-\theta_k) \Im \bar{F}^o_k  \\
                    &= \Re F^e_k - \cos(\theta_k) \Re F^o_k + \sin(\theta_k) \Im F^o_k \qquad \text{Line a.update(i3, h1r - tr)?} \\
    \Im F_{N/2 - k} &= - \Im \bar{F}^e_k + \cos(-\theta_k) \Im \bar{F}^o_k - \sin(-\theta_k) \Re \bar{F}^o_k \\
                    &= - \Im F^e_k + \cos(\theta_k) \Im F^o_k + \sin(\theta_k) \Re F^o_k \qquad \text{Line a.update(i4, -h1i + ti)?}
    $$
    """
    let size = a.size()

    ifdef debug then
      try
        Assert((size and (size - 1)) == 0, "The input array 'a' size (" +
              size.string() + ") must be a power of 2. Enlarge it to " +
              size.next_pow2().string(), true)?
        Assert(size > 0, "Array 'a' can't be empty", true)?
        Assert(normalize or inverse, "`normalize` must be " +
              "set to `false` only with inverse FFT", true)?
      end
    end

    // Constants
    let pi = F.from[F64](F64.pi())
    let theta = if inverse then
        pi / F.from[USize](size / 2)
      else
        -pi / F.from[USize](size / 2)
      end
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    let c1 = one / two
    let c2 = if inverse then c1 else -c1 end

    if not inverse then
      fourier_complex(a)
    end

    let temp = (theta / two).sin()
    let multiplier_r = -two * temp * temp
    let multiplier_i = theta.sin()
    var factor_r = one + multiplier_r
    var factor_i = multiplier_i

    try
      for i in Range(1, size / 4) do
        let i1 = i + i
        let i2 = i1 + 1
        let i3 = size - i1
        let i4 = i3 + 1

        // The two transforms are separated
        let h1r = c1 * (a(i1)? + a(i3)?)
        let h1i = c1 * (a(i2)? - a(i4)?)
        let h2r = -c2 * (a(i2)? + a(i4)?)
        let h2i = c2 * (a(i1)? - a(i3)?)

        // And recombined to form the transform of the original data
        let tr = (factor_r * h2r) - (factor_i * h2i)
        let ti = (factor_r * h2i) + (factor_i * h2r)
        a.update(i1, h1r + tr)?
        a.update(i2, h1i + ti)?
        a.update(i3, h1r - tr)?
        a.update(i4, -h1i + ti)?

        // Trigonometric recurrence
        (factor_r, factor_i) = (factor_r + ((factor_r * multiplier_r) -
                                            (factor_i * multiplier_i)),
                                factor_i + ((factor_i * multiplier_r) +
                                            (factor_r * multiplier_i)))
      end

      let hr = a(0)?
      if inverse then
        a.update(0, c1 * (hr + a(1)?))?
        a.update(1, c1 * (hr - a(1)?))?
	      a.update((size / 2) + 1, -a((size / 2) + 1)?)?
        fourier_complex(a, true, false)
        if normalize then
          let scale = two / F.from[USize](size)
          for i in Range(0, size) do
            a.update(i, a(i)? * scale)?
          end
        end
      else
        a.update(0, hr + a(1)?)?
        a.update(1, hr - a(1)?)?
	      a.update((size / 2) + 1, -a((size / 2) + 1)?)?
      end
    else
      Fail("Index out of bounds [0.." + size.string() + ")")
    end
    // Return the array
    a


  fun fourier2(a: Array[Complex[F val]],
               inverse: Bool = false,
               normalize: Bool = true)
              : Array[Complex[F]] =>
    """
    Calculate the discrete fast Fourier tranform of `a`. It calculates the
    inverse transform when `inverse = true`. Calculation is done in-place,
    with the content of `a` being replaced by the result.

    The size of `a` must be a power of 2.

    When calculating the inverse transform, one can decide not to normalize
    the result by setting `normalize` parameter to `false`, and getting
    relative values.

    This version builds a table of `a.size()` trigonometric values that are
    used for the calculation. This consumes more memory but results should have
    less calculation errors because errors don't propagate in the trigonometric
    recurrence of `fourier`. It also allocates `Complex[F]` objects.
    """
    let size = a.size()
    ifdef debug then
      try
        Assert((size and (size - 1)) == 0, "The input array 'a' size (" +
              size.string() + ") must be a power of 2. Enlarge it to " +
              size.next_pow2().string(), true)?
        Assert(size > 0, "Array 'a' can't be empty", true)?
        Assert(normalize or inverse, "`normalize` must be " +
            "set to `false` only with inverse FFT", true)?
      end
    end

    // Prepare trigonometric tables: calculate the complex roots or 1
    let pi2 = if inverse then
        F.from[USize](2) * F.from[F64](F64.pi())
      else
        -F.from[USize](2) * F.from[F64](F64.pi())
      end
    let size2 = size / 2
    let root1 = Array[Complex[F]](size2)
    for i in Range[USize](0, size2) do
      let theta = (pi2 * F.from[USize](i)) / F.from[USize](size)
      root1.push(Complex[F](theta.cos(), theta.sin()))
    end

    // Reorganize the array content in bits order
    _rearrange[Complex[F val]](a)

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
            let ai = a(i)?
            let product = a(j)? * root1(k)?
            a.update(j, ai - product)?
            a.update(i, ai + product)?
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


  fun bluestein(a: Array[Complex[F]],
                inverse: Bool = false,
		normalize: Bool = true)
	       : Array[Complex[F]] =>
    """
    Calculate the discrete Fourier transform of `a` array that can have any size,
    not necessarily a power of 2 like with `fourier` methods. Calculation is
    done in-place and content of `a` is replaced by the result.
    
    Calculate the inverse DFT when `inverse` is `true`. The inverse is
    normalized when `normalize` is `true` (the default).

    This implementation uses the
    [Bluestein's Chirp Z-Transform algorithm](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm).
    It is not particularly memory efficient as it allocates multiple temporary
    vectors and uses the `Complex[F]]` class.
    """
    let size = a.size()
    // Find a power-of-2 convolution length m such that m >= 2 * size - 1
    let m = 2 * size.next_pow2()

    // Prepare trigonometric tables
    let pi = if inverse then
               F.from[F64](F64.pi())
	     else
	       -F.from[F64](F64.pi())
	     end
    let root1 = Array[Complex[F]](size)
    for i in Range[USize](0, size) do
      let theta = (pi * F.from[USize]((i * i) % (2 * size))) / F.from[USize](size)
      root1.push(Complex[F](theta.cos(), theta.sin()))
    end
   
    try
      // Temporary vectors
      let u = Array[Complex[F]].init(Complex[F], m)
      for i in Range[USize](0, size) do
        //let re = (a(i)?.real() * root1(i)?.real()) + (a(i)?.imag() * root1(i)?.imag())
      	//let im = (-a(i)?.real() * root1(i)?.imag()) + (a(i)?.imag() * root1(i)?.real())
        //u.update(i, Complex[F](re, im))?
	      u.update(i, a(i)? * root1(i)?)?
      end
      
      let v = Array[Complex[F]].init(Complex[F], m)
      v.update(0, root1(0)?)?
      for i in Range[USize](1, size) do
        let c = root1(i)?.conj()
        v.update(i, c)?
        v.update(m - i, c)?
      end
      
      // Convolution
      let w = convolve(u, v)

      // Postprocessing
      for i in Range[USize](0, size) do
        //let re = (w(i)?.real() * root1(i)?.real()) + (w(i)?.imag() * root1(i)?.imag())
	      //let im = (-w(i)?.real() * root1(i)?.imag()) + (w(i)?.imag() * root1(i)?.real())
        //a.update(i, Complex[F](re, im))?
	      a.update(i, w(i)? * root1(i)?)?
      end

      // If inverse, scale the result if normalize is required
      if inverse and normalize then
        _normalize(a)
      end
    else
      Fail("Index error when accessing arrays")
    end
    
    // Now a contains the DFT
    a


  fun convolve(u: Array[Complex[F]], v: Array[Complex[F]]): Array[Complex[F]]^ =>
    """
    Calculates the circular convolution of `u` and `v`. Both arrays must have
    the same length and be a power of 2.
    
    The circular convolution is
    $$result(j) = \sum_{k=0}^{n-1} u((j - k) mod n) v(k)$$
    
    More information on https://en.wikipedia.org/wiki/Circular_convolution
    """
    let size = u.size()
    ifdef debug then
      try
        Assert(v.size() == size, "Array sizes must be identical: u.size = " +
              size.string() + ", v.size = " + v.size().string(), true)?
      end
    end

    let u' = u.clone()
    let v' = v.clone()
    
    fourier(u')
    fourier(v')
    
    try
      for i in Range[USize](0, size) do
        u'.update(i, u'(i)? * v'(i)?)?
      end
    else
      Fail("Index access error")
    end
    
    // Inverse FFT
    fourier(u', true)
    
    // u' contains the convolution now
    u'
    