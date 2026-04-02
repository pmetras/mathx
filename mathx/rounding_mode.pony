// Definition of the floating point rounding modes, based on the definitions
// of GMP and IEEE 754.


type RoundingMode is (RoundingNearest | RoundingNegInf | RoundingPosInf |
                      RoundingZero | RoundingAwayZ | RoundingFaithful )
  """
  The rounding mode that must be applied to operations.

  * MPFloat_RNDN: round to nearest, with the even rounding rule (roundTiesToEven in IEEE 754).
  * MPFloat_RNDD: round toward negative infinity (roundTowardNegative in IEEE 754).
  * MPFloat_RNDU: round toward positive infinity (roundTowardPositive in IEEE 754).
  * MPFloat_RNDZ: round toward zero (roundTowardZero in IEEE 754).
  * MPFloat_RNDA: round away from zero.
  * MPFloat_RNDF: faithful rounding. This feature is currently experimental.

  For full definition and background: https://www.mpfr.org/mpfr-current/mpfr.html#Rounding
  and https://en.wikipedia.org/wiki/Rounding
  """


primitive RoundingNearest
    """
    Round to the nearest, with the even rounding rule.
    A number `x` is rounded to the number `y` that is the closest to `x` such that
    in case the number to be rounded lies exactly in the middle between two consecutive
    representable numbers, it is rounded to the one with an even significand.
    """
  fun apply(): I32 =>
    0

  fun string(): String =>
    "RoundingNearest"


primitive RoundingNegInf
    """
    Round toward negative infinity.
    A number `x` is rounded to the number `y` that is the closest to `x` such that
    `y` is less than or equal to `x`.
    """
  fun apply(): I32 =>
    3

  fun string(): String =>
    "RoundingNegInf"


primitive RoundingPosInf
    """
    Round toward positive infinity.
    A number `x` is rounded to the number `y` that is the closest to `x` such that
    `y` is greater than or equal to `x`.
    """
  fun apply(): I32 =>
    2

  fun string(): String =>
    "RoundingPosInf"


primitive RoundingZero
    """
    Round toward zero.
    A number `x` is rounded to the number `y` that is the closest to `x` such that
    `abs(y)` is less than or equal to `abs(x)`.
    """
  fun apply(): I32 =>
    1

  fun string(): String =>
    "RoundingZero"


primitive RoundingAwayZ
    """
    Round away from zero.
    A number `x` is rounded to the number `y` that is the closest to `x` such that
    `abs(y)` is greater than or equal to `abs(x)`.
    """
  fun apply(): I32 =>
    4

  fun string(): String =>
    "RoundingAwayZ"


primitive RoundingFaithful
    """
    Faithfull rounding.
    A number `x` is rounded to the number `y` that is the closest to `x` such that
    the computed value is either that corresponding to `MPFloat_RNDD` or that
    corresponding to `MPFloat_RNDU`.
    """
  fun apply(): I32 =>
    5

  fun string(): String =>
    "RoundingFaithful"

