"""
A Fibonacci generator
"""
actor Main
  new create(env: Env) =>
    let n = try
      env.args(1)?.usize()?
    else
      10
    end
    env.out.print("Fibonacci(" + n.string() + ") = " + fibonacci(n).string())

  fun fibonacci(n: USize): USize =>
    """
    Calculate the Fibonnaci sequence:
      F(n) = F(n -1) + F(n - 2)
    with F(0) = 0 and F(1) = 1
    """
    if n == 0 then
      0
    elseif n == 1 then
      1
    else
      fibonacci(n - 1) + fibonacci(n - 2)
    end

"""
A distributed Fibonacci generator
"""
actor Main

  new create(env: Env) =>
    let n = try
        env.args(1)?.usize()?
      else
        3
      end

    // Distributed Fibonacci calculation on multiple actors, splitting the load
    // on all CPU cores.
    var fibo = Fibonacci
    fibo.calculate(n, {(x: USize) => env.out.print("Fibonacci(" + n.string() + ") = " + x.string())})

actor Fibonacci
  var _result: USize = 0                               // Result of calculation
  var _nb_add: USize = 0                               // Count partial results
  var _call_me: {(USize)} val = {(x: USize) => None}   // Dummy callback for init

  be accumulate(n: USize) =>
    """
    Called by an actor who as completed its part of the calculation. We add the
    `n` value to the result and we check if we have gathered all sub-calculations.
    If that the case, we can report the result to the initial caller.
    """
    _result = _result + n
    _nb_add = _nb_add + 1
    if _nb_add == 2 then
      // The two children have reported their part of the calculation
      _call_me(_result)
    end

  be calculate(n: USize, callback: {(USize)} val) =>
    if n == 0 then
      // Tell caller that Fibonacci(0) = 0
      callback(0)
    elseif n == 1 then
      // Tell caller that Fibonacci(1) = 1
      callback(1)
    else
      // We don't know the result: delegate the sub-calculations to 2 actors
      _call_me = callback
      let child1 = Fibonacci
      let child2 = Fibonacci
      child1.calculate(n - 1, {(x: USize)(self: Fibonacci tag = this) => self.accumulate(x)})
      child2.calculate(n - 2, {(x: USize)(self: Fibonacci tag = this) => self.accumulate(x)})
    end
  