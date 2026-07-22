// Collatz conjecture investigation with MPInt

use "../../bignum"


actor Main
  """
  Main entry point for the Collatz investigation.
  """

  new create(env: Env) =>
    """
    Initializes the Collatz investigation starting from N=1 and prints instructions.
    """
    env.out.print("Collatz (Syracuse) conjecture investigation using MPInt.")
    env.out.print("The program will run indefinitely. Use [Ctrl]+[C] to stop.")
    env.out.print("---------------------------------------------------------")

    let runner = CollatzRunner(env.out)

    runner.next(MPInt.from[ILong](1))


actor CollatzRunner
  """
  Analyzes the Collatz sequence for increasing values of N.
  """

  let _out: OutStream
    """
    The output stream to print recaps.
    """

  let _one: MPInt = MPInt.from[ILong](1)
    """
    The MPInt value 1.
    """

  let _two: MPInt = MPInt.from[ILong](2)
    """
    The MPInt value 2.
    """

  let _three: MPInt = MPInt.from[ILong](3)
    """
    The MPInt value 3.
    """

  let _ten: MPInt = MPInt.from[ILong](10)
    """
    The MPInt value 10.
    """

  var _max_n: MPInt = MPInt.from[ILong](1)
    """
    The value of N that produced the maximum number of steps so far.
    """

  var _max_steps: U64 = 0
    """
    The maximum number of steps found so far.
    """

  var _max_val_reached: MPInt = MPInt.from[ILong](1)
    """
    The maximum value ever reached in any sequence so far.
    """

  var _max_val_n: MPInt = MPInt.from[ILong](1)
    """
    The value of N that led to the maximum value reached so far.
    """

  var _next_power_of_10: MPInt = MPInt.from[ILong](10)
    """
    The next power of 10 at which to print a recap.
    """

  var _power: U64 = 1
    """
    The current exponent of the power of 10.
    """


  new create(out: OutStream) =>
    """
    Creates a new CollatzRunner.
    """
    _out = out


  be next(n: MPInt) =>
    """
    Analyzes the Collatz sequence for N and schedules the next one.
    """
    if n == _next_power_of_10 then
      _print_recap()
      _next_power_of_10 = _next_power_of_10 * _ten
      _power = _power + 1
    end

    (let steps, let local_max) = _analyze_sequence(n)

    if steps > _max_steps then
      _max_n = n
      _max_steps = steps
    end

    if local_max > _max_val_reached then
      _max_val_reached = local_max
      _max_val_n = n
    end

    next(n + _one)


  fun ref _analyze_sequence(n: MPInt): (U64, MPInt) =>
    """
    Calculates the number of steps and the maximum value reached for the
    Collatz sequence starting at N.
    """
    var current = n
    var steps: U64 = 0
    var max_val = n

    while not current.is_one() do
      if (current.u64() % 2) == 0 then
        current = current / _two
      else
        current = (current * _three) + _one
      end

      if current > max_val then
        max_val = current
      end

      steps = steps + 1
    end

    (steps, max_val)


  fun _print_recap() =>
    """
    Prints the summary of the maximum steps and maximum value reached for N 
    less than the current power of 10.
    """
    let msg = recover
      let s = String
      s.append("less than ")

      if _power < 4 then
        s.append(_next_power_of_10.string())
      else
        s.append("10^")
        s.append(_power.string())
      end

      s.append(" is ")
      s.append(_max_n.string())
      s.append(", which has ")
      s.append(_max_steps.string())
      s.append(" steps; max value reached was ")
      s.append(_max_val_reached.string())
      s.append(" for N=")
      s.append(_max_val_n.string())
      s.append(",")
      s
    end

    _out.print(consume msg)
