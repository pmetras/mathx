// Common interface for arbitrary-precision floating-point numbers.

use "../mathx"

interface val BigReal[T: BigReal[T] val]
  """
  Common interface for arbitrary-precision floating-point numbers.

  Both `mathx.MPFloat` (pure-Pony base-256 implementation) and
  `mathx/gmp.MPFloat` (MPFR binding) satisfy this interface, allowing
  generic algorithms to work with either backend:

  ```pony
  class val MPFloat is BigReal[MPFloat]
  ```

  The self-referential type parameter `T` preserves concrete return types so
  that arithmetic on `MPFloat` still returns `MPFloat`, not a boxed interface.
  Generic algorithms are written as:

  ```pony
  fun foo[T: BigReal[T] val](x: T, y: T): T =>
    x.add(y).sqrt()
  ```

  ## Rounding

  The interface declares methods without an explicit rounding-mode parameter.
  `gmp.MPFloat` methods carry `rnd: RoundingMode = RoundingNearest` as a
  default parameter; Pony structural typing allows a concrete method with extra
  default parameters to satisfy an interface method with fewer parameters.

  ## Special values

  All implementations represent IEEE-754-style special values: NaN, ±infinity,
  and ±zero. Use the `is_*` predicates to query these states.
  """

  //- Precision and metadata --------------------------------------------------

  fun get_precision(): USize
    """Precision of the mantissa in bits."""

  fun get_rounding_mode(): RoundingMode
    """The rounding mode in use."""

  fun get_base(): USize
    """Internal encoding base (256 for mathx.MPFloat, 2 for gmp.MPFloat)."""

  //- Special-value predicates ------------------------------------------------

  fun is_nan(): Bool
    """Is this value not-a-number?"""

  fun is_infinite(): Bool
    """Is this value ±infinity?"""

  fun is_finite(): Bool
    """Is this value finite (not NaN, not infinite)?"""

  fun is_zero(): Bool
    """Is this value zero (+0 or −0)?"""

  fun is_negative(): Bool
    """Is this value strictly negative?"""

  fun is_integer(): Bool
    """Is this value a finite integer?"""

  //- Conversions -------------------------------------------------------------

  fun string(): String iso^
    """Decimal string representation."""

  fun exact_string(base: U8 = 10, rnd: RoundingMode = RoundingNearest)
                  : (String iso^, I64, Bool)
    """
    Return `(mantissa, exponent, inexact)` where the value equals
    `0.mantissa × base^exponent`.
    """

  fun f64(): F64
    """Convert to F64 (may lose precision)."""

  fun f32(): F32
    """Convert to F32 (may lose precision)."""

  //- Sign and comparison -----------------------------------------------------

  fun sign(): Compare
    """Return Less / Equal / Greater for negative / zero / positive."""

  fun eq(that: T): Bool
    """Is `this == that`?"""

  fun ne(that: T): Bool
    """Is `this != that`?"""

  fun lt(that: T): Bool
    """Is `this < that`?"""

  fun le(that: T): Bool
    """Is `this <= that`?"""

  fun ge(that: T): Bool
    """Is `this >= that`?"""

  fun gt(that: T): Bool
    """Is `this > that`?"""

  fun compare(that: T): Compare
    """Three-way comparison: Less, Equal, or Greater."""

  //- Arithmetic --------------------------------------------------------------

  fun add(that: T): T
    """Addition `this + that`."""

  fun sub(that: T): T
    """Subtraction `this - that`."""

  fun mul(that: T): T
    """Multiplication `this × that`."""

  fun div(that: T): T
    """Division `this / that`."""

  fun divrem(that: T): (T, T)
    """Quotient and remainder of `this / that`."""

  fun neg(): T
    """Negation `-this`."""

  fun abs(): T
    """Absolute value `|this|`."""

  fun inv(): T
    """Multiplicative inverse `1/this`."""

  fun sqrt(): T
    """Square root."""

  fun cbrt(): T
    """Cube root."""

  fun rootn(n: USize): T
    """n-th root."""

  fun min(that: T): T
    """Minimum of `this` and `that`."""

  fun max(that: T): T
    """Maximum of `this` and `that`."""

  fun next_above(): T
    """Next representable value above `this`."""

  fun next_below(): T
    """Next representable value below `this`."""

  //- Rounding ----------------------------------------------------------------

  fun trunc(): T
    """Round toward zero."""

  fun floor(): T
    """Round toward −∞."""

  fun ceil(): T
    """Round toward +∞."""

  fun round(): T
    """Round to nearest integer, halfway away from zero."""

  //- Transcendentals ---------------------------------------------------------

  fun ln(): T
    """Natural logarithm ln(this)."""

  fun log(): T
    """Natural logarithm (alias for ln)."""

  fun log2(): T
    """Base-2 logarithm."""

  fun log10(): T
    """Base-10 logarithm."""

  fun exp(): T
    """e^this."""

  fun exp2(): T
    """2^this."""

  fun pow(that: T): T
    """this^that."""

  fun powi(n: ILong): T
    """this^n for integer exponent n."""

  fun sin(): T
    """Sine."""

  fun cos(): T
    """Cosine."""

  fun tan(): T
    """Tangent."""

  fun asin(): T
    """Arc-sine."""

  fun acos(): T
    """Arc-cosine."""

  fun atan(): T
    """Arc-tangent."""

  fun sinh(): T
    """Hyperbolic sine."""

  fun cosh(): T
    """Hyperbolic cosine."""

  fun tanh(): T
    """Hyperbolic tangent."""

  fun asinh(): T
    """Inverse hyperbolic sine."""

  fun acosh(): T
    """Inverse hyperbolic cosine."""

  fun atanh(): T
    """Inverse hyperbolic tangent."""
