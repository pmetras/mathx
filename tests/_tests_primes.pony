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

