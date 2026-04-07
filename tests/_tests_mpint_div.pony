
use "../mathx"
use "../pony_testx"
use "collections"
use "random"

class iso _TestMPIntDivision is UnitTest
  """
  Comprehensive tests for MPInt division (Algorithm D).
  """

  fun name(): String =>
    "MPInt/division"

  fun apply(h: TestHelper) =>
    let rand = Rand()
    
    // 1. Basic properties: (a / b) * b + (a % b) == a
    for _ in Range(0, 100) do
      let a: MPInt = _random_mpint(rand, 10 + (rand.u16() % 50).usize())
      let b: MPInt = _random_mpint(rand, 1 + (rand.u16() % 10).usize())
      
      if b.is_zero() then continue end
      
      (let q, let r) = a.divrem(b)
      
      // Property check
      let check = (q * b) + r
      if check != a then
        h.log("FAIL Property: a=" + a.dump() + " (" + a.string() + ")")
        h.log("              b=" + b.dump() + " (" + b.string() + ")")
        h.log("              q=" + q.dump() + " (" + q.string() + ")")
        h.log("              r=" + r.dump() + " (" + r.string() + ")")
        let qb = q * b
        h.log("              q*b=" + qb.dump())
        h.log("              (q*b)+r=" + check.dump())
        h.assert_true(check == a, "Division property: (a/b)*b + r == a")
      end
      
      // Remainder magnitude check
      if not r.abs_lt(b.abs()) then
        h.log("FAIL Magnitude: |r| >= |b|")
        h.log("               r=" + r.dump() + " (" + r.string() + ")")
        h.log("               b=" + b.dump() + " (" + b.string() + ")")
        h.assert_true(r.abs_lt(b.abs()), "Remainder magnitude: |r| < |b|")
      end
      
      // Remainder sign check (Truncated division: remainder has same sign as dividend)
      if not r.is_zero() then
        let z = MPInt.from[ILong](0)
        h.assert_true((r < z) == (a < z), "Remainder sign: sign(r) == sign(a)")
      end
    end

    // 2. Division by 1 and -1
    let large: MPInt = _random_mpint(rand, 100)
    let one = MPInt.from[ILong](1)
    let m_one = MPInt.from[ILong](-1)
    
    (let q1, let r1) = large.divrem(one)
    h.assert_true(q1 == large, "a / 1 == a")
    h.assert_true(r1.is_zero(), "a % 1 == 0")
    
    (let q2, let r2) = large.divrem(m_one)
    h.assert_true(q2 == large.neg(), "a / -1 == -a")
    h.assert_true(r2.is_zero(), "a % -1 == 0")

    // 3. Small numbers (compared with ILong arithmetic)
    for _ in Range(0, 100) do
      let ai: I64 = rand.i64() / 2 
      let bi: I64 = (rand.i64() % 1000).abs().i64() + 1
      
      let a: MPInt = MPInt.from[ILong](ai.ilong())
      let b: MPInt = MPInt.from[ILong](bi.ilong())
      
      (let q, let r) = a.divrem(b)
      h.assert_true(q == MPInt.from[ILong]((ai / bi).ilong()), "Small div matches ILong")
      h.assert_true(r == MPInt.from[ILong]((ai % bi).ilong()), "Small rem matches ILong")
    end

    // 4. Dividend < Divisor
    try
      let a_small: MPInt = MPInt.from[ILong](500)
      let b_large: MPInt = MPInt.from_string("1000000000000000000000000")?
      (let q3, let r3) = a_small.divrem(b_large)
      h.assert_true(q3.is_zero(), "Small / Large == 0")
      h.assert_true(r3 == a_small, "Small % Large == Small")
    else
      h.fail("Error in large number string conversion")
    end

    // 5. Division by zero (should not crash, follows divrem implementation)
    let zero = MPInt.from[ILong](0)
    (let q0, let r0) = large.divrem(zero)
    h.assert_true(q0.is_zero() and r0.is_zero(), "Division by zero returns (0,0)")

  fun _random_mpint(rand: Rand, size: USize): MPInt =>
    var res = MPInt.from[ILong](rand.ilong())
    for i in Range(1, size) do
      res = res + (MPInt.from[ILong](rand.ilong()).abs().digit_shl(i))
    end
    if (rand.next() % 2) == 0 then
      res.neg()
    else
      res
    end
