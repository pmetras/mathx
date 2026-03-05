"""
Display platform limits for F32 and F64 types
"""

use "debug"

use "../../mathx"


actor Main
  let _out: OutStream

  new create(env: Env) =>
    """
    Print the platform limits for floating-point numbers.
    """
    _out = env.out
    _out.print("Platform limits")
    _out.print("===============")
    let f32_limits = FLimits[F32]
    print_limits[F32](f32_limits)
    let f64_limits = FLimits[F64]
    print_limits[F64](f64_limits)

  fun print_limits[F: FloatingPoint[F] val](fl: FLimits[F]) =>
    """
    Prints the numerical limits of the given floating-point type.
    """
    iftype F <: F32 then
      _out.print("Calculated limits for type F32")
    end
    iftype F <: F64 then
      _out.print("Calculated limits for type F64")
    end
    _out.print("--------------------------------")
    _out.print("On binary hardware (radix == 2), unit = bit")

    _out.print("")
    _out.print("Radix            = " + fl.radix().string())
    _out.print("Digits           = " + fl.digit().string() + " units")
    _out.print("Guard digits     = " + fl.guard_digits().string())
    _out.print("\t0: if floating-point arithmetic rounds, or if it truncates" +
               " and only `Digit`-base `Radix` digits participate in the" +
               " post-normalization shift of the floating-point mantissa in" +
               " multiplication.")
    _out.print("\t1: if floating-point arithmetic truncates and more than" +
               " `Digit`-base `Radix` digits participate in the post-" +
               " normalization shift of the floating-point mantissa in" +
               " multiplication.")

    _out.print("Round style      = " + fl.round_style().string())
    _out.print("\t0: if floating-point addition chops")
    _out.print("\t1: if floating-point addition rounds, but not in the IEEE" +
               " style")
    _out.print("\t2: if floating-point addition rounds in the IEEE style")
    _out.print("\t3: if floating-point addition chops, and there is partial" +
               " underflow")
    _out.print("\t4: if floating-point addition rounds, but not in the IEEE" +
               " style, and there is partial underflow")
    _out.print("\t5: if floating-point addition rounds in the IEEE style, and" +
               " there is partial underflow")

    _out.print("Machep           = " + fl.machep().string())
    _out.print("\tThe exponent for the smallest power of `Radix` (but bounded" +
               " below by `Digits - 3`) whose sum with 1.0 is greater than 1.0.")
    _out.print("Negeps           = " + fl.negeps().string())
    _out.print("\tThe exponent for the smallest power of `Radix` (but bounded" +
               " below by `Digits - 3`) whose difference with 1.0 is less than 1.0.")
    _out.print("Exponent         = " + fl.exponent().string() + " units")
    _out.print("Minimal exponent = " + fl.min_exponent().string())
    _out.print("Maximal exponent = " + fl.max_exponent().string())
    _out.print("Epsilon          = " + fl.epsilon().string())
    _out.print("Negative epsilon = " + fl.negative_epsilon().string())
    _out.print("Minimal value    = " + fl.min_value().string())
    _out.print("Maximal value    = " + fl.max_value().string())

    _out.print("")
    iftype F <: F32 then
      _out.print("Compiler limits for type F32")
    end
    iftype F <: F64 then
      _out.print("Compiler limits for type F64")
    end
    _out.print("------------------------------")
    let one = F.from[ISize](1)

    _out.print("F.precision2     = " + one.precision2().string() + " bits")
    _out.print("F.precision10    = " + one.precision10().string())
    _out.print("F.min_exp2       = " + one.min_exp2().string())
    _out.print("F.min_exp10      = " + one.min_exp10().string())
    _out.print("F.max_exp2       = " + one.max_exp2().string())
    _out.print("F.max_exp10      = " + one.max_exp10().string())
    _out.print("F.epsilon        = " + F.epsilon().string())
    _out.print("F.min_value      = " + F.min_value().string())
    _out.print("F.max_value      = " + F.max_value().string())
    _out.print("F.min_normalised = " + F.min_normalised().string())
    _out.print("")

