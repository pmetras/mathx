// Tests for modular arithmetic

use "collections"
use "random"
use "time"

use "../mathx"
use "../pony_testx"


class _TestGCDLCM is UnitTest
  fun name(): String =>
    "modular/gcd-lcm"

  fun apply(h: TestHelper) =>
    // --- Modular.gcd (Euclidean) ---
    h.assert_eq[USize](21, Modular.gcd(1071, 462))
    h.assert_eq[USize](28, Modular.gcd(2156, 3220))
    h.assert_eq[USize](5, Modular.gcd(153215, 2100235))
    h.assert_eq[USize](16, Modular.gcd(16, 32))
    h.assert_eq[USize](1, Modular.gcd(2, 3))
    // zero arguments
    h.assert_eq[USize](5, Modular.gcd(0, 5))
    h.assert_eq[USize](5, Modular.gcd(5, 0))
    h.assert_eq[USize](0, Modular.gcd(0, 0))
    // equal arguments
    h.assert_eq[USize](7, Modular.gcd(7, 7))
    // commutativity
    h.assert_eq[USize](Modular.gcd(1071, 462), Modular.gcd(462, 1071))

    // --- Modular.gcd2 (binary GCD) ---
    h.assert_eq[USize](21, Modular.gcd2(1071, 462))
    h.assert_eq[USize](28, Modular.gcd2(2156, 3220))
    h.assert_eq[USize](5, Modular.gcd2(153215, 2100235))
    h.assert_eq[USize](16, Modular.gcd2(16, 32))
    h.assert_eq[USize](1, Modular.gcd2(2, 3))
    // zero arguments
    h.assert_eq[USize](5, Modular.gcd2(0, 5))
    h.assert_eq[USize](5, Modular.gcd2(5, 0))
    h.assert_eq[USize](0, Modular.gcd2(0, 0))
    // equal arguments
    h.assert_eq[USize](7, Modular.gcd2(7, 7))
    // commutativity
    h.assert_eq[USize](Modular.gcd2(1071, 462), Modular.gcd2(462, 1071))

    // --- Modular.lcm ---
    h.assert_eq[USize](23562, Modular.lcm(1071, 462))
    h.assert_eq[USize](247940, Modular.lcm(2156, 3220))
    h.assert_eq[USize](64357501105, Modular.lcm(153215, 2100235))
    h.assert_eq[USize](32, Modular.lcm(16, 32))
    h.assert_eq[USize](6, Modular.lcm(2, 3))
    // zero arguments
    h.assert_eq[USize](0, Modular.lcm(0, 0))
    h.assert_eq[USize](0, Modular.lcm(0, 5))
    h.assert_eq[USize](0, Modular.lcm(5, 0))
    // equal arguments
    h.assert_eq[USize](7, Modular.lcm(7, 7))


class _TestModularInverse is UnitTest
  fun name(): String =>
    "modular/inverse"

  fun apply(h: TestHelper) =>
    // Inverse modulo prime number
    for i in Range[U32](1, 997) do
      h.assert_ne[U32](0, Modular[U32].inverse_mod(i, 997), "Modulo prime")
    end


class _TestModularMPIntBasic is UnitTest
  """
  Verifies every Modular operation using MPInt, with the same values as the
  USize tests so results can be cross-checked by inspection.
  """
  fun name(): String =>
    "modular/mpint/basic"

  fun mp(n: USize): MPInt => MPInt.from[USize](n)

  fun apply(h: TestHelper) =>
    let m0 = mp(0)
    let m1 = mp(1)
    let m7 = mp(7)

    // add_mod
    h.assert_eq[MPInt](mp(1), Modular[MPInt].add_mod(mp(3), mp(5), m7))   // (3+5)%7=1
    h.assert_eq[MPInt](mp(5), Modular[MPInt].add_mod(mp(6), mp(6), m7))   // (6+6)%7=5
    h.assert_eq[MPInt](m0,    Modular[MPInt].add_mod(m0, m0, m7))

    // sub_mod
    h.assert_eq[MPInt](mp(2), Modular[MPInt].sub_mod(mp(5), mp(3), m7))   // 5-3=2
    h.assert_eq[MPInt](mp(5), Modular[MPInt].sub_mod(mp(3), mp(5), m7))   // 3-5 ≡ 5 (mod 7)

    // neg_mod
    h.assert_eq[MPInt](mp(4), Modular[MPInt].neg_mod(mp(3), m7))          // -3 ≡ 4 (mod 7)
    h.assert_eq[MPInt](m0,    Modular[MPInt].neg_mod(m0, m7))

    // mul_mod — MPInt always takes the binary-doubling path (clz always 0)
    h.assert_eq[MPInt](m1,    Modular[MPInt].mul_mod(mp(3), mp(5), m7))   // 15%7=1
    h.assert_eq[MPInt](m1,    Modular[MPInt].mul_mod(mp(6), mp(6), m7))   // 36%7=1

    // gcd (Euclidean)
    h.assert_eq[MPInt](mp(21), Modular[MPInt].gcd(mp(1071), mp(462)))
    h.assert_eq[MPInt](m1,     Modular[MPInt].gcd(mp(2), mp(3)))
    h.assert_eq[MPInt](m0,     Modular[MPInt].gcd(m0, m0))
    // commutativity
    h.assert_eq[MPInt](Modular[MPInt].gcd(mp(1071), mp(462)),
                       Modular[MPInt].gcd(mp(462), mp(1071)))

    // gcd2 (binary GCD)
    h.assert_eq[MPInt](mp(21), Modular[MPInt].gcd2(mp(1071), mp(462)))
    h.assert_eq[MPInt](m1,     Modular[MPInt].gcd2(mp(2), mp(3)))
    // zero edge cases
    h.assert_eq[MPInt](mp(5),  Modular[MPInt].gcd2(m0, mp(5)))
    h.assert_eq[MPInt](mp(5),  Modular[MPInt].gcd2(mp(5), m0))

    // lcm
    h.assert_eq[MPInt](mp(23562), Modular[MPInt].lcm(mp(1071), mp(462)))
    h.assert_eq[MPInt](m0,        Modular[MPInt].lcm(m0, m0))
    h.assert_eq[MPInt](m0,        Modular[MPInt].lcm(m0, mp(5)))

    // inverse_mod: 3 * 5 ≡ 1 (mod 7)
    h.assert_eq[MPInt](mp(5), Modular[MPInt].inverse_mod(mp(3), m7))
    // no inverse when not coprime
    h.assert_eq[MPInt](m0, Modular[MPInt].inverse_mod(mp(2), mp(4)))

    // div_mod: 6 / 3 ≡ 2 (mod 7)
    h.assert_eq[MPInt](mp(2), Modular[MPInt].div_mod(mp(6), mp(3), m7))

    // pow_mod
    h.assert_eq[MPInt](mp(24), Modular[MPInt].pow_mod(mp(2), mp(10), mp(100))) // 2^10%100=24
    h.assert_eq[MPInt](m1,     Modular[MPInt].pow_mod(mp(2), m0, m7))          // a^0 = 1


class _TestModularMPIntBig is UnitTest
  """
  Demonstrates Modular[MPInt] handling values beyond any fixed-width integer.

  Ground truth is provided by two number-theory identities that can be verified
  without computing the intermediate products:

  1. Fermat's little theorem: a^(p-1) ≡ 1 (mod p) for prime p, gcd(a,p)=1.

  2. Exponent reduction: since p ≡ 1 (mod p-1), we have p^3 ≡ 1 (mod p-1),
     so 2^(p^3) = 2^(k*(p-1)+1) ≡ 2 (mod p).
     p^3 ≈ 10^27 cannot be stored in any fixed-width integer type.
  """
  fun name(): String =>
    "modular/mpint/big"

  fun apply(h: TestHelper) =>
    let p   = MPInt.from[U64](1_000_000_007)   // well-known prime
    let one = MPInt.from[USize](1)
    let two = MPInt.from[USize](2)

    // Fermat's little theorem: 2^(p-1) ≡ 1 (mod p)
    h.assert_eq[MPInt](one,
      Modular[MPInt].pow_mod(two, p - one, p))

    // Exponent p^3 ≈ 10^27 exceeds U64 capacity; result is 2 by Fermat + exponent reduction
    h.assert_eq[MPInt](two,
      Modular[MPInt].pow_mod(two, p * p * p, p))

    // Modular inverse round-trip: a * a^(-1) ≡ 1 (mod p)
    let a   = MPInt.from[U64](123_456_789)
    let inv = Modular[MPInt].inverse_mod(a, p)
    h.assert_ne[MPInt](MPInt.from[USize](0), inv, "123456789 must be invertible mod p")
    h.assert_eq[MPInt](one, Modular[MPInt].mul_mod(a, inv, p))


