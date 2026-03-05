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
    h.assert_eq[USize](21, Modular.gcd(1071, 462))
    h.assert_eq[USize](28, Modular.gcd(2156, 3220))
    h.assert_eq[USize](5, Modular.gcd(153215, 2100235))
    h.assert_eq[USize](16, Modular.gcd(16, 32))
    h.assert_eq[USize](1, Modular.gcd(2, 3))

    h.assert_eq[USize](21, Modular.gcd2(1071, 462))
    h.assert_eq[USize](28, Modular.gcd2(2156, 3220))
    h.assert_eq[USize](5, Modular.gcd2(153215, 2100235))
    h.assert_eq[USize](16, Modular.gcd2(16, 32))
    h.assert_eq[USize](1, Modular.gcd2(2, 3))

    h.assert_eq[USize](23562, Modular.lcm(1071, 462))
    h.assert_eq[USize](247940, Modular.lcm(2156, 3220))
    h.assert_eq[USize](64357501105, Modular.lcm(153215, 2100235))
    h.assert_eq[USize](32, Modular.lcm(16, 32))
    h.assert_eq[USize](6, Modular.lcm(2, 3))


class _TestModularInverse is UnitTest
  fun name(): String =>
    "modular/inverse"

  fun apply(h: TestHelper) =>
    // Inverse modulo prime number
    for i in Range[U32](1, 997) do
      h.assert_ne[U32](0, Modular[U32].inverse_mod(i, 997), "Modulo prime")
    end


