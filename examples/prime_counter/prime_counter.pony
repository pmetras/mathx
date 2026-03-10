// Demonstration of the prime number theorem

use "../../mathx"
use "math"

use @fprintf[I32](stream: Pointer[U8] tag, fmt: Pointer[U8] tag, ...)
use @pony_os_stdout[Pointer[U8]]()


actor Main
  """
  Inifinite calculation of the number of prime numbers.
  """

  let _out: OutStream
    """
    The output stream to print prime counting status.
    """


  new create(env: Env) =>
    """
    Initializes the prime counter and selects the calculation strategy.
    """
    _out = env.out

    _out.print("Prime counting function. Experimental proof of the prime number theorem...")
    _out.print("--------------------------------------------------------------------------")
    _out.print("Arguments:")
    _out.print("\tNone:\tProbabilistic test with actor messages")
    _out.print("\t1:\tProbabilistic prime test in loop")
    _out.print("\t2:\tPony stdlib primality test in loop")
    _out.print("")
    _out.print("The number of prime numbers in the interval [1 .. n], named `Pi(n)`, tends to")
    _out.print("`n/log(n)` when `n` tends to the infinite.")
    _out.print("")
    _out.print("Stop execution with [Ctrl]+[C] when tired of the results.")

    // If we give an argument to the command line, then we run the fast version
    // that does not send messages (but leaks memory!).
    if env.args.size() > 1 then
      try
        match env.args(1)?
        | "1" => fast_find(Prime[U128]~is_probably_prime(where k = 13))
        | "2" => fast_find(IsPrime[U128]~apply())
        else
          _out.print("Argument must be 1 for probablistic prime test")
          _out.print("Argument must be 2 for Pony primality test using 6k ± 1 method.")
        end
      end
    else
      find_primes(3, 1)
    end


  be find_primes(n: U128, pi_n: U128) =>
    """
    Find a batch of prime numbers in the range `[n .. n + 100_000)`,
    print results and call itself to continue.

    Every 1_000_000 numbers, it prints the status of the calculation. As `print`
    send a message to another actor, we could saturate the queue of that
    actor if we run an infinite loop, and in the end consume all memory by
    filling up that actor's queue with print messages. We break the infinite
    loop into a tail-recursive behaviour: we do the calculation ony on a batch
    of numbers, then send a message to itself to continue the calculation. We
    call the print actor only 1% of the time we call ourself, so Pony backpressure
    can cope with the different actors.
    """
    // Let's us be ambitious!
    var n' = n
    var pi_n' = pi_n
    let limit = n + 100_000

    while n' < limit do
      if Prime[U128].is_probably_prime(n') then
        pi_n' = pi_n' + 1
      end
      if (n' % 1_000_000) == 1 then
        let n_log_n: F64 = F64.from[U128](n') / F64.from[U128](n').log()
        let pi_n_log_n: F64 = F64.from[U128](pi_n') / n_log_n
        _out.print("n = " + n'.string() + "\tPi(n) = " + pi_n'.string() +
                   "\tn/log(n) = " + n_log_n.string() +
                   "\t1 - Pi(n)/(n/log(n)) = " + (1 - pi_n_log_n).string())
      end

      // No need to test for even numbers
      n' = n' + 2
    end

    // Continue work from limit
    find_primes(n', pi_n')


  fun fast_find(primer: {(U128): Bool}) =>
    """
    Without sending messages to other actors (the OutStream), we can find primes
    much faster.
    """
    @fprintf(@pony_os_stdout(), "FAST non-actors version\n".cstring())

    var n: U128 = 3
    var pi_n: U128 = 1

    while true do
      if primer(n) then
        pi_n = pi_n + 1
      end
      if (n % 1_000_000) == 1 then
        let n_log_n: F64 = F64.from[U128](n) / F64.from[U128](n).log()
        let pi_n_log_n: F64 = F64.from[U128](pi_n) / n_log_n
        let msg = "n = " + n.string() + "\tPi(n) = " + pi_n.string() +
                  "\tn/log(n) = " + n_log_n.string() +
                  "\t1 - Pi(n)/(n/log(n)) = " + (1 - pi_n_log_n).string()
        @fprintf(@pony_os_stdout(), "%s\n".cstring(), msg.cstring())
      end
      n = n + 2

      // The while true loop in fast_find is a synchronous fun that never
      // returns to the Pony scheduler. All String allocations (from .string()
      // calls and + concatenation) accumulate without GC, growing at a rate
      // of ~several allocations per million numbers tested. The breakpoint
      // at n > 1_000_000_000 is the workaround.
      if n > 1_000_000_000 then
        @fprintf(@pony_os_stdout(), "Completed at n = %s\n".cstring(), n.string().cstring())
        break
      end
    end

