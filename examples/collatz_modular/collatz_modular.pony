// Collatz conjecture investigation with Modular[ULong]

use "../../mathx"


actor Main
  """
  Main entry point for the Collatz investigation using modular arithmetic.
  """

  new create(env: Env) =>
    """
    Initializes the Collatz investigation starting from N=1 and prints instructions.
    """
    env.out.print("Collatz (Syracuse) conjecture investigation using Modular[ULong].")
    env.out.print("The program will run indefinitely. Use [Ctrl]+[C] to stop.")
    env.out.print("---------------------------------------------------------")

    let runner = CollatzRunner(env.out)

    runner.next(1)


actor CollatzRunner
  """
  Analyzes the Collatz sequence for increasing values of N using ULong.
  """

  let _out: OutStream
    """
    The output stream to print recaps.
    """

  let _mod: Modular[ULong] = Modular[ULong]
    """
    The modular arithmetic utility for ULong.
    """

  let _limit: ULong = ULong.max_value()
    """
    The modulus limit for calculations, representing the range of ULong.
    """

  var _max_n: ULong = 1
    """
    The value of N that produced the maximum number of steps so far.
    """

  var _max_steps: U64 = 0
    """
    The maximum number of steps found so far.
    """

  var _max_val_reached: ULong = 1
    """
    The maximum value ever reached in any sequence so far.
    """

  var _max_val_n: ULong = 1
    """
    The value of N that led to the maximum value reached so far.
    """

  var _next_power_of_10: ULong = 10
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


  be next(n: ULong) =>
    """
    Analyzes the Collatz sequence for N and schedules the next one.
    """
    if n == _next_power_of_10 then
      _print_recap()
      _next_power_of_10 = _next_power_of_10 * 10
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

    next(n + 1)


  fun ref _analyze_sequence(n: ULong): (U64, ULong) =>
    """
    Calculates the number of steps and the maximum value reached for the
    Collatz sequence starting at N using Modular[ULong].
    """
    var current = n
    var steps: U64 = 0
    var max_val = n

    while current != 1 do
      // Use modular arithmetic to check parity: if (current % 2) == 0
      if _mod.sub_mod(current, 0, 2) == 0 then
        current = current / 2
      else
        // current = 3n + 1
        // Perform the step modulo ULong.max_value()
        current = _mod.add_mod(_mod.mul_mod(current, 3, _limit), 1, _limit)
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
