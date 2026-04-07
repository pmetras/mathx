// Parallel Collatz conjecture investigation with MPInt

use "../../mathx"
use "collections"


actor Main
  """
  Main entry point for the parallel Collatz investigation.
  Analyzes the Collatz sequence for increasing values of N using multiple worker actors.
  """

  new create(env: Env) =>
    """
    Prints initial instructions and starts the Master coordinator.
    """
    env.out.print("Parallel Collatz investigation using MPInt and Pony actors.")
    env.out.print("Using standard integer operations for stability.")
    env.out.print("The program will run indefinitely. Use [Ctrl]+[C] to stop.")
    env.out.print("---------------------------------------------------------")

    let master = Master(env.out)
    master.start()


actor Master
  """
  Coordinates worker actors to analyze Collatz sequences in parallel.
  Dispatches work in batches and aggregates results.
  """

  let _out: OutStream
    """
    The output stream to print recaps.
    """

  let _ten: MPInt = MPInt.from[ILong](10)
    """
    The MPInt value 10.
    """

  let _batch_size: USize = 500
    """
    Number of sequences to analyze in each worker batch.
    """

  let _max_pending: USize = 32
    """
    Maximum number of concurrent worker actors to limit memory usage.
    """

  var _current_n: MPInt = MPInt.from[ILong](1)
    """
    The next N to be analyzed.
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
    Creates a new Master with the given output stream.
    """
    _out = out


  be start() =>
    """
    Dispatches the initial set of worker actors.
    """
    for _ in Range(0, _max_pending) do
      _dispatch_batch()
    end


  be report_results(batch_max_n: MPInt, 
                    batch_max_steps: U64, 
                    batch_max_val: MPInt, 
                    batch_max_val_n: MPInt,
                    last_n: MPInt) =>
    """
    Aggregates results from a worker and dispatches a new worker for the next batch.
    """
    // Update global maximums
    if batch_max_steps > _max_steps then
      _max_n = batch_max_n
      _max_steps = batch_max_steps
    end

    if batch_max_val > _max_val_reached then
      _max_val_reached = batch_max_val
      _max_val_n = batch_max_val_n
    end

    // Check for periodic recap
    if last_n >= _next_power_of_10 then
      _print_recap()
      _next_power_of_10 = _next_power_of_10 * _ten
      _power = _power + 1
    end

    _dispatch_batch()


  fun ref _dispatch_batch() =>
    """
    Dispatches the next batch of numbers to a new worker actor.
    """
    let start_n = _current_n
    let end_n = _current_n + MPInt.from[ILong](_batch_size.ilong())
    _current_n = end_n

    let worker = Worker(this)
    worker.analyze_range(start_n, end_n)


  fun _print_recap() =>
    """
    Prints the summary of the maximum steps and maximum value reached.
    """
    let msg = recover
      let s = String
      s.append("progress recap: analyzed around ")

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


actor Worker
  """
  Calculates Collatz sequences for a range of N. 
  A new worker is created for each batch to maximize scheduling flexibility.
  """

  let _master: Master
    """
    The master actor to report results back to.
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


  new create(master: Master) =>
    """
    Creates a new Worker attached to the given master.
    """
    _master = master


  be analyze_range(start_n: MPInt, end_n: MPInt) =>
    """
    Analyzes all N in the range [start_n, end_n).
    Reports results back to the master when finished.
    """
    var max_n = start_n
    var max_steps: U64 = 0
    var max_val = start_n
    var max_val_n = start_n

    var n = start_n
    while n < end_n do
      (let steps, let local_max) = _analyze_sequence(n)

      if steps > max_steps then
        max_n = n
        max_steps = steps
      end

      if local_max > max_val then
        max_val = local_max
        max_val_n = n
      end

      n = n + _one
    end

    _master.report_results(max_n, max_steps, max_val, max_val_n, end_n)


  fun _analyze_sequence(n: MPInt): (U64, MPInt) =>
    """
    Calculates the number of steps and the maximum value reached for the
    Collatz sequence starting at N using standard integer arithmetic.
    """
    var current = n
    var steps: U64 = 0
    var peak = n

    while not current.is_one() do
      // Use u64() for parity check as it is reliable and fast for the low bit
      if (current.u64() % 2) == 0 then
        // even: current = n / 2
        current = current / _two
      else
        // odd: current = 3n + 1
        current = (current * _three) + _one
      end

      if current > peak then
        peak = current
      end

      steps = steps + 1
    end

    (steps, peak)
