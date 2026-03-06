// Tests for primes

use "collections"
use "random"
use "time"

use "../mathx"
use "../pony_testx"


class _TestPrimeTest is UnitTest
  // everyone should know them
  let primes_expected: Array[USize] val = [2; 3; 5; 7; 11; 13; 17; 19; 23; 29
    31; 37; 41; 43; 47]

  fun name(): String =>
    "prime/test"
    
  fun apply(h: TestHelper) =>
    var primes: Array[USize] trn = []
    for n in Range[USize](0, 50) do
      if Prime.is_prime(n) then primes.push(n) end
    end
    h.assert_array_eq[USize](primes_expected, consume val primes, "Prime.is_prime error!")


class _TestPrimeTestExtended is UnitTest
  // you don't have to learn these :)
  let primes_expected: Array[USize] val = [10007; 10009; 10037; 10039; 10061
    10067; 10069; 10079; 10091; 10093; 10099; 10103; 10111; 10133; 10139; 10141
    10151; 10159; 10163; 10169; 10177; 10181; 10193; 10211; 10223; 10243; 10247
    10253; 10259; 10267; 10271; 10273; 10289; 10301; 10303; 10313; 10321; 10331
    10333; 10337; 10343; 10357; 10369; 10391; 10399; 10427; 10429; 10433; 10453
    10457; 10459; 10463; 10477; 10487; 10499; 10501; 10513; 10529; 10531; 10559
    10567; 10589; 10597; 10601; 10607; 10613; 10627; 10631; 10639; 10651; 10657
    10663; 10667; 10687; 10691; 10709; 10711; 10723; 10729; 10733; 10739; 10753
    10771; 10781; 10789; 10799; 10831; 10837; 10847; 10853; 10859; 10861; 10867
    10883; 10889; 10891; 10903; 10909; 10937; 10939; 10949; 10957; 10973; 10979
    10987; 10993]

  fun name(): String =>
    "prime/test/advanced"
    
  fun apply(h: TestHelper) =>
    var primes: Array[USize] trn = []
    for n in Range[USize](10000, 11000) do
      if Prime.is_prime(n) then primes.push(n) end
    end
    h.assert_array_eq[USize](primes_expected, consume val primes, "Prime.is_prime error!")


class _TestNextPrime is UnitTest
  fun name(): String =>
    "prime/nextprime"
    
  fun apply(h: TestHelper) =>
    for n in Range[USize](1, 100) do
      let from = n * 100
      let next_prime = Prime.next_prime(from)
      for n' in Range[USize](from + 1, next_prime) do
        h.assert_false(Prime.is_prime(n'), "Prime.next_prime error! " + n'.string() + " - " + next_prime.string())
      end
    end


class _TestPrimeFactorsStatic is UnitTest
  fun name(): String =>
    "prime/factors/static"
    
  fun apply(h: TestHelper) =>
    h.assert_array_eq[USize]([], Prime.prime_factors(1))
    h.assert_array_eq[USize]([as USize: 2], Prime.prime_factors(2))
    h.assert_array_eq[USize]([as USize: 3], Prime.prime_factors(3))
    h.assert_array_eq[USize]([as USize: 2; 2], Prime.prime_factors(4))
    h.assert_array_eq[USize]([as USize: 2; 2; 7; 11; 17; 31],
                             Prime.prime_factors(162316))
    h.assert_array_eq[USize]([as USize: 2; 3; 5; 7; 11; 13; 17; 19; 23; 29],
                             Prime.prime_factors(6469693230))
    h.assert_array_eq[USize]([as USize: 2; 3; 5; 7; 7; 17; 31; 31; 43; 47; 10007],
                             Prime.prime_factors(485690777622330))
    h.assert_array_eq[U128]([as U128: 2; 3; 5; 7; 11; 13; 17; 19; 23; 29; 31; 37; 41; 43; 47],
                            Prime[U128].prime_factors(614889782588491410))
    h.assert_array_eq[U128]([as U128: 2; 3; 5; 7; 11; 13; 17; 19; 23; 29
                            31; 37; 41; 43; 47; 53; 59; 61; 67; 71; 73; 79
                            83; 89; 97; 101],
                            Prime[U128].prime_factors(232862364358497360900063316880507363070))


class _TestPrimeFactorsRandom is UnitTest
  fun name(): String =>
    "prime/factors/random"

  fun apply(h: TestHelper) =>
    let now = Time.now()
    let dice = Dice(XorOshiro128Plus(now._1.u64(), now._2.u64()))
    
    for x in Range[USize](0, 8) do
      let to: USize = dice(1, 4).usize()
      let iterator = PrimeIterator
      var primes: Array[USize] trn = []
      var composite: USize = 1
      for y in Range[USize](0, to) do
        var value = iterator.next()
        for z in Range[USize](0, dice(1, 6).usize()) do
          value = iterator.next()
        end
        for z in Range[USize](0, dice(1, 2).usize()) do
          primes.push(value)
          composite = composite * value
        end
      end
      h.assert_array_eq[USize](Prime.prime_factors(composite), consume primes)
    end


class _TestPrimeIterator is UnitTest
  let primes_expected: Array[USize] box = [10007; 10009; 10037; 10039; 10061
    10067; 10069; 10079; 10091; 10093; 10099; 10103; 10111; 10133; 10139; 10141
    10151; 10159; 10163; 10169; 10177; 10181; 10193; 10211; 10223; 10243; 10247
    10253; 10259; 10267; 10271; 10273; 10289; 10301; 10303; 10313; 10321; 10331
    10333; 10337; 10343; 10357; 10369; 10391; 10399; 10427; 10429; 10433; 10453
    10457; 10459; 10463; 10477; 10487; 10499; 10501; 10513; 10529; 10531; 10559
    10567; 10589; 10597; 10601; 10607; 10613; 10627; 10631; 10639; 10651; 10657
    10663; 10667; 10687; 10691; 10709; 10711; 10723; 10729; 10733; 10739; 10753
    10771; 10781; 10789; 10799; 10831; 10837; 10847; 10853; 10859; 10861; 10867
    10883; 10889; 10891; 10903; 10909; 10937; 10939; 10949; 10957; 10973; 10979
    10987; 10993]

  fun name(): String => "prime/iterator"

  fun apply(h: TestHelper) =>
    let iterator = PrimeIterator
    h.assert_eq[USize](2, iterator.next(), "PrimeIterator error!")
    h.assert_eq[USize](3, iterator.next(), "PrimeIterator error!")
    h.assert_eq[USize](5, iterator.next(), "PrimeIterator error!")

    let iterator' = PrimeIterator.start_at(10000)
    let primes: Array[USize] = []
    var value = iterator'.next()
    while value < 11000 do
      primes.push(value = iterator'.next())
    end
    h.assert_array_eq[USize](primes_expected, primes, "PrimeIterator error!")


class _TestCoprime is UnitTest
  fun name(): String =>
    "prime/coprime"

  fun apply(h: TestHelper) =>
    h.assert_true(Prime.is_coprime(1, 1))
    h.assert_true(Prime.is_coprime(1, 2))
    h.assert_true(Prime.is_coprime(2, 3))
    h.assert_true(Prime.is_coprime(3, 5))
    h.assert_true(Prime.is_coprime(7, 12))
    h.assert_false(Prime.is_coprime(2, 6))
    h.assert_false(Prime.is_coprime(3, 93))
    h.assert_false(Prime.is_coprime(12, 93))


class _TestProbablyPrime is UnitTest
  fun name(): String =>
    "prime/probably prime"
    
  fun apply(h: TestHelper) =>
    let one = U128(1)

    // Fermat F6 is composite
    h.assert_false(Prime[U128].is_probably_prime(18446744073709551617))
    // Mersenne M11
    h.log("M11 = " + (UnsignedComp[U128].pow2(107) - one).string())
    h.assert_true(Prime[U128].is_probably_prime(UnsignedComp[U128].pow2(107) - one))
    // Mersenne M12
    h.log("M12 = " + (UnsignedComp[U128].pow2(127) - one).string())
    h.assert_true(Prime[U128].is_probably_prime(UnsignedComp[U128].pow2(127) - one))
    // Perfect number P9
    h.log("P9 = " + (UnsignedComp[U128].pow2(60) * (UnsignedComp[U128].pow2(61) - one)).string())
    h.assert_false(Prime[U128].is_probably_prime(UnsignedComp[U128].pow2(60) * (UnsignedComp[U128].pow2(61) - one)))
    // Prime with all digits
    h.assert_true(Prime[U128].is_probably_prime(100123456789))
    // No prime in this range
    for i in Range[ULong](1671781 + 1, 1671907) do
      h.assert_false(Prime[ULong].is_probably_prime(i, 5))
    end



class _TestProbablyPrimeEdgeCases is UnitTest
  fun name(): String =>
    "prime/probably prime/edge cases"

  fun apply(h: TestHelper) =>
    // Bug fix: 0 and 1 are not prime
    h.assert_false(Prime.is_probably_prime(0))
    h.assert_false(Prime.is_probably_prime(1))
    // Bug fix: 2 and 3 are prime
    h.assert_true(Prime.is_probably_prime(2))
    h.assert_true(Prime.is_probably_prime(3))
    // Even numbers > 2 are not prime
    h.assert_false(Prime.is_probably_prime(4))
    h.assert_false(Prime.is_probably_prime(100))


class _TestPrimeFactorsPerf is UnitTest
  fun name(): String =>
    "prime/factors/large prime"

  fun apply(h: TestHelper) =>
    // 2 * large_prime: previously would iterate PrimeIterator all the way to
    // large_prime; now terminates after sqrt(large_prime) steps
    h.assert_array_eq[USize](
      [as USize: 2; 10007],
      Prime.prime_factors(2 * 10007))
    h.assert_array_eq[USize](
      [as USize: 2; 3; 10007],
      Prime.prime_factors(6 * 10007))
    // Product of two large primes
    h.assert_array_eq[USize](
      [as USize: 10007; 10009],
      Prime.prime_factors(10007 * 10009))


class _TestPrimeFactorsUnique is UnitTest
  fun name(): String =>
    "prime/factors/unique"

  fun apply(h: TestHelper) =>
    h.assert_array_eq[USize]([], Prime.prime_factors_unique(0))
    h.assert_array_eq[USize]([], Prime.prime_factors_unique(1))
    h.assert_array_eq[USize]([as USize: 2], Prime.prime_factors_unique(2))
    h.assert_array_eq[USize]([as USize: 2], Prime.prime_factors_unique(4))
    h.assert_array_eq[USize]([as USize: 2], Prime.prime_factors_unique(8))
    h.assert_array_eq[USize]([as USize: 2; 3], Prime.prime_factors_unique(12))
    h.assert_array_eq[USize]([as USize: 2; 3], Prime.prime_factors_unique(18))
    h.assert_array_eq[USize]([as USize: 2; 3; 5], Prime.prime_factors_unique(30))
    h.assert_array_eq[USize]([as USize: 7], Prime.prime_factors_unique(49))
    // Same as prime_factors for squarefree numbers
    h.assert_array_eq[USize]([as USize: 2; 3; 5; 7], Prime.prime_factors_unique(210))


class _TestEulerTotient is UnitTest
  fun name(): String =>
    "prime/euler_totient"

  fun apply(h: TestHelper) =>
    // φ(1) = 1
    h.assert_eq[USize](1, Prime.euler_totient(1))
    // φ(p) = p - 1 for prime p
    h.assert_eq[USize](1, Prime.euler_totient(2))
    h.assert_eq[USize](2, Prime.euler_totient(3))
    h.assert_eq[USize](4, Prime.euler_totient(5))
    h.assert_eq[USize](6, Prime.euler_totient(7))
    // φ(p^k) = p^(k-1) * (p-1)
    h.assert_eq[USize](2, Prime.euler_totient(4))   // φ(2²) = 2
    h.assert_eq[USize](4, Prime.euler_totient(8))   // φ(2³) = 4
    h.assert_eq[USize](6, Prime.euler_totient(9))   // φ(3²) = 6
    // φ(mn) = φ(m)*φ(n) when gcd(m,n)=1
    h.assert_eq[USize](2, Prime.euler_totient(6))   // φ(6) = φ(2)*φ(3) = 1*2 = 2
    h.assert_eq[USize](4, Prime.euler_totient(10))  // φ(10) = φ(2)*φ(5) = 1*4 = 4
    h.assert_eq[USize](4, Prime.euler_totient(12))  // φ(12) = φ(4)*φ(3) = 2*2 = 4
    h.assert_eq[USize](8, Prime.euler_totient(24))  // φ(24) = φ(8)*φ(3) = 4*2 = 8


class _TestRadical is UnitTest
  fun name(): String =>
    "prime/radical"

  fun apply(h: TestHelper) =>
    h.assert_eq[USize](1, Prime.radical(0))  // by convention (no prime factors)
    h.assert_eq[USize](1, Prime.radical(1))
    h.assert_eq[USize](2, Prime.radical(2))
    h.assert_eq[USize](2, Prime.radical(4))   // 2²
    h.assert_eq[USize](2, Prime.radical(8))   // 2³
    h.assert_eq[USize](6, Prime.radical(12))  // 2²*3 → radical = 2*3
    h.assert_eq[USize](30, Prime.radical(30)) // squarefree: radical = itself
    h.assert_eq[USize](6, Prime.radical(36))  // 2²*3² → radical = 2*3
    h.assert_eq[USize](7, Prime.radical(49))  // 7²


class _TestIsSquarefree is UnitTest
  fun name(): String =>
    "prime/is_squarefree"

  fun apply(h: TestHelper) =>
    h.assert_true(Prime.is_squarefree(1))
    h.assert_true(Prime.is_squarefree(2))
    h.assert_true(Prime.is_squarefree(3))
    h.assert_true(Prime.is_squarefree(6))   // 2*3
    h.assert_true(Prime.is_squarefree(30))  // 2*3*5
    h.assert_false(Prime.is_squarefree(4))  // 2²
    h.assert_false(Prime.is_squarefree(8))  // 2³
    h.assert_false(Prime.is_squarefree(12)) // 2²*3
    h.assert_false(Prime.is_squarefree(18)) // 2*3²
    h.assert_false(Prime.is_squarefree(36)) // 2²*3²


class _TestPrevPrime is UnitTest
  fun name(): String =>
    "prime/prev_prime"

  fun apply(h: TestHelper) =>
    // No prime before 2
    h.assert_eq[USize](0, Prime.prev_prime(0))
    h.assert_eq[USize](0, Prime.prev_prime(1))
    h.assert_eq[USize](0, Prime.prev_prime(2))
    // prev_prime(3) = 2
    h.assert_eq[USize](2, Prime.prev_prime(3))
    h.assert_eq[USize](3, Prime.prev_prime(4))
    h.assert_eq[USize](3, Prime.prev_prime(5))
    h.assert_eq[USize](5, Prime.prev_prime(6))
    h.assert_eq[USize](5, Prime.prev_prime(7))
    h.assert_eq[USize](7, Prime.prev_prime(8))
    h.assert_eq[USize](7, Prime.prev_prime(9))
    h.assert_eq[USize](7, Prime.prev_prime(10))
    h.assert_eq[USize](11, Prime.prev_prime(12))
    // Verify prev_prime and next_prime are inverses around a prime.
    // next_prime(2) = 2 (PrimeIterator.start_at(2) quirk), so we skip n=3.
    for n in Range[USize](5, 100) do
      if Prime.is_prime(n) then
        let prev = Prime.prev_prime(n)
        let next = Prime.next_prime(prev)
        h.assert_eq[USize](n, next, "prev/next inverse failed at " + n.string())
      end
    end


class _TestPrimeSieve is UnitTest
  let primes_to_50: Array[USize] val = [2; 3; 5; 7; 11; 13; 17; 19; 23; 29
    31; 37; 41; 43; 47]

  fun name(): String =>
    "prime/sieve"

  fun apply(h: TestHelper) ? =>
    // Basic sieve up to 50
    let sieve = PrimeSieve(50)
    h.assert_eq[USize](50, sieve.limit())
    h.assert_array_eq[USize](primes_to_50, sieve.all_primes())
    h.assert_eq[USize](15, sieve.count())

    // is_prime queries within range
    h.assert_false(sieve.is_prime(0)?)
    h.assert_false(sieve.is_prime(1)?)
    h.assert_true(sieve.is_prime(2)?)
    h.assert_true(sieve.is_prime(47)?)
    h.assert_false(sieve.is_prime(48)?)
    h.assert_false(sieve.is_prime(49)?)
    h.assert_false(sieve.is_prime(50)?)

    // Out of range raises an error
    h.assert_true(
      try sieve.is_prime(51)?; false else true end,
      "is_prime(51) should error for sieve(50)")

    // Agrees with is_prime (trial division) for [0, 50]
    for n in Range[USize](0, 51) do
      h.assert_eq[Bool](
        Prime.is_prime(n),
        sieve.is_prime(n)?,
        "sieve disagrees with is_prime at " + n.string())
    end

    // Edge cases: empty sieves
    let sieve0 = PrimeSieve(0)
    h.assert_eq[USize](0, sieve0.count())
    h.assert_true(
      try sieve0.is_prime(2)?; false else true end,
      "is_prime(2) on sieve(0) should error")

    let sieve1 = PrimeSieve(1)
    h.assert_eq[USize](0, sieve1.count())

    let sieve2 = PrimeSieve(2)
    h.assert_eq[USize](1, sieve2.count())
    h.assert_true(sieve2.is_prime(2)?)

    // Larger sieve: count of primes <= 1000 is 168
    let sieve1000 = PrimeSieve(1000)
    h.assert_eq[USize](168, sieve1000.count())
    // Verify all primes found match is_prime
    for p in sieve1000.primes() do
      h.assert_true(Prime.is_prime(p),
        "sieve returned non-prime " + p.string())
    end
