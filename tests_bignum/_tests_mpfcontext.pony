use "../bignum"
use "../pony_testx"


class iso _TestMPFContextCreate is UnitTest
  """
  Verify that `MPFContext` stores precision and rounding mode correctly and that
  `p_bytes()` converts bits to bytes (rounding up).
  """

  fun name(): String => "MPFContext/create"


  fun apply(h: TestHelper) =>
    // Default: 128 bits = 16 bytes.
    let ctx = MPFContext
    h.assert_eq[USize](ctx.precision, 128)
    h.assert_eq[USize](ctx.p_bytes(), 16)

    // 113 bits rounds up to 15 bytes.
    let ctx2 = MPFContext(113)
    h.assert_eq[USize](ctx2.p_bytes(), 15)

    // 8 bits = 1 byte.
    let ctx3 = MPFContext(8)
    h.assert_eq[USize](ctx3.p_bytes(), 1)

    // 256 bits = 32 bytes.
    let ctx4 = MPFContext(256)
    h.assert_eq[USize](ctx4.p_bytes(), 32)


class iso _TestMPFContextWorkingBytes is UnitTest
  """
  Verify that `working_bytes` adds the correct guard bytes for known operation
  names and falls back to the default guard for unknown names.
  """

  fun name(): String => "MPFContext/working_bytes"


  fun apply(h: TestHelper) =>
    // MPFContext(112): p_bytes = 14, log_guard = floor(log₂(14))+1 = 4.
    let ctx = MPFContext(112)
    let p: USize = ctx.p_bytes()  // 14
    // At p=14, log_guard=4 is >= every floor except ln(6) and exp/trig/pi(8).
    h.assert_eq[USize](ctx.working_bytes("add"),     p + 4)  // max(2, 4)
    h.assert_eq[USize](ctx.working_bytes("sub"),     p + 4)  // max(2, 4)
    h.assert_eq[USize](ctx.working_bytes("mul"),     p + 4)  // max(4, 4)
    h.assert_eq[USize](ctx.working_bytes("inv"),     p + 4)  // max(4, 4)
    h.assert_eq[USize](ctx.working_bytes("sqrt"),    p + 4)  // max(4, 4)
    h.assert_eq[USize](ctx.working_bytes("ln"),      p + 6)  // max(6, 4)
    h.assert_eq[USize](ctx.working_bytes("exp"),     p + 8)  // max(8, 4)
    h.assert_eq[USize](ctx.working_bytes("trig"),    p + 8)  // max(8, 4)
    h.assert_eq[USize](ctx.working_bytes("pi"),      p + 8)  // max(8, 4)
    h.assert_eq[USize](ctx.working_bytes("unknown"), p + 4)  // max(4, 4)

    // MPFContext(10000): p_bytes = 1250, log_guard = 11.
    // At p=1250, log_guard=11 dominates all floors.
    let ctx2 = MPFContext(10000)
    let p2: USize = ctx2.p_bytes()  // 1250
    h.assert_eq[USize](ctx2.working_bytes("add"),     p2 + 11)
    h.assert_eq[USize](ctx2.working_bytes("mul"),     p2 + 11)
    h.assert_eq[USize](ctx2.working_bytes("ln"),      p2 + 11)
    h.assert_eq[USize](ctx2.working_bytes("pi"),      p2 + 11)
