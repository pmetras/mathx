// Defines a `MPFloat` implementation

interface Precision
  """
  A `Precision` defines a category of `MPFloat` that can be used
  for calculations. It defines the precision of the floating point numbers.
  """

  fun precision(): ILong
    """
    Return the precision in bits used to represent the floating-point number.
    """


  fun interm_precision(name: String): ILong =>
    """
    Give the required intermediate precision required to get correct rounding when
    doing the calculations of the function `name`. The default implementation
    add 10-bits guard to the value of `precision`, for all functions (does not
    take `name` parameter into account).

    For example, here are the required intermediate precision values for correct
    minimal error with `F64`, for the following functions:

    * `exp`: 159
    * `cos`: 143
    * `sin`: 127
    * `asin`: 127
    * `tan`: 133
    * `atan`: 127
    * `sinh`: 127
    * `asinh`: 127
    * `log`: 119
    * `2^x`: 114
    * `log2`: 110
    * `acos`: 117
    * `cosh`: 112
    * `acosh`: 116
    """
    precision() + 10


class MPFloatB128 is Precision
  """
  Defines a category of `MPFloat` that is similar to a 128-bits floating point.
  """

  fun precision(): ILong =>
    """
    The significand of a 128-bits floating point number is 112 bits.
    """
    112


  fun interm_precision(name: String): ILong =>
    """
    Give the required precision for interpediate calculations of function
    `name` (the guard binary digits to prevent error propagation if we want to
    have correction `precision`).
    """
    let prec = precision()
    match name
    | "pi" => prec + 32
    | "pi_chudnovsky" => 2 * prec
    | "pi_pbb" => prec + 32
    else
      prec
    end