# Counting the number of prime numbers

The objective of this example is to calculate the **prime counting function** `Pi(n)` that counts the number of prime numbers. It counts the number of prime numbers starting from 1, using the probabilistic test, and evaluates the ratio `Pi(n)/(x / log(n))` to see if it converges to 1, as expected by the prime number theorem, as demonstrated by Jacques Adamard and Charles de la Vallée Poussin. It prints the partial results every 1,000,000 numbers, until you stop the program with Ctrl^C.

* [Wikipedia page](https://en.wikipedia.org/wiki/Prime-counting_function)
* [Counting primes](http://numbers.computation.free.fr/Constants/Primes/countingPrimes.html)

```
Prime counting function. Experimental proof of the prime number theorem...
--------------------------------------------------------------------------
The number of prime numbers in the interval [1 .. n], named `Pi(n)`, tends to
`n/log(n)` when `n` tends to the infinite.

n = 1000001	Pi(n) = 78498	n/log(n) = 72382.5	1 - Pi(n)/(n/log(n)) = -0.0844889
n = 2000001	Pi(n) = 148933	n/log(n) = 137849	1 - Pi(n)/(n/log(n)) = -0.0804085
n = 3000001	Pi(n) = 216816	n/log(n) = 201152	1 - Pi(n)/(n/log(n)) = -0.0778732
n = 4000001	Pi(n) = 283146	n/log(n) = 263127	1 - Pi(n)/(n/log(n)) = -0.0760823
n = 5000001	Pi(n) = 348513	n/log(n) = 324150	1 - Pi(n)/(n/log(n)) = -0.0751588
...
```

## Command line

Prime test is done by the statistical test that is run 13 times by number.

Without arguments, it uses a tail-recursive loop as explained below, to let Pony runtime manage actors and memory.

Using the argument `1`, it does the same tests but using an infinite `while true...` loop. The prime calculations have been optimized not to use memory so that they can run without Pony garbage collector. That is that they don't use any `class` or `actor` instances, but only `primitive`. As I can't guarantee that they don't consume memory (though checked with Balgrind Massif), a maximum execution limit with `n < 1_000_000_000` has been set.

Using the argument `2`, the prime test uses `math/IsPrime` object instead and it still use an infinite `while` loop. For small values of `n`, this test is more performant than the probabilistic test. To assert that a number is a prime, the statistical test must be run a certain number of times (13 times in this example, giving `(1/4)^n` chances that the number is compound), while the exact test runs `O(sqrt(n))` at the maximum, and much less in practice. But for large values of `n`, the probabilistic test should take the advantage.

For both infinite loop versions, we don't use the `OutStream` actor to print status but we bypass using another actor by direct calls to FFP `fprintf`.

## Long running process

It is dangerous to write an infinite loop in Pony as it prevents Pony mechanism managing the runtime (backpressure, garbage collection, etc.). If messages are sent to other actors (through behaviour calls), then it can prevent these actors processing the messages and they stack up in the actors queues, consumming the memory and could kill the OS. This can happen when you call `env.out.print()` for instance as printing messages is done by sending messages to the `OutStream` actor that can't keep up with the stream of messages.

An infinite loop must be rewritten with tail-recursion.

```pony
  var i: USize = 1
  while true do
    ... do something depending on i ...
    i = i + 1
  end
```

becomes

```pony
  do_work(1)

  be do_work(i: USize)
    ... do something depending on i ...
    do_work(i + 1)
```

The backpressure mechanism must kickoff to share CPU time between actors and let the garbage collector reclaim memory.

