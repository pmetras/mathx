// Multi-precision float context: precision policy and rounding.

use "../assertx"


class val MPFContext
  """
  Execution context for arbitrary-precision floating-point arithmetic.

  An `MPFContext` binds together two orthogonal pieces of policy:

  - **`precision`** — the number of significant bits in any result produced
    under this context.  Stored in bits (user-facing unit); converted to bytes
    internally via `p_bytes()`.
  - **`rounding`** — the IEEE 754 / MPFR rounding mode applied when a
    computed result is shortened to `precision` bits.

  ## Usage pattern

  Every arithmetic operation follows this two-phase discipline:

  1. Compute at *working precision* `w = working_bytes(op)` bytes — enough
     guard bytes to bound the rounding error of the operation to ≤ 1 ULP.
  2. Round back to *output precision* `p = p_bytes()` bytes via `_round_to`.

  `MPFloat` holds an `MPFContext` and delegates this discipline automatically.
  Power users can also use `MPFContext` directly with `MPFRep` operands once
  the `MPFContext` arithmetic API is complete (Step 4 of the architectural
  split).

  ## Sentinel `tag` values

  `p_bytes()` and `working_bytes` require a live instance.  Any code path that
  truly needs a context without an instance (e.g. `tag` factory methods) should
  construct a default `MPFContext` with `MPFContext.create()`.

  ## Default precision

  The default precision is 112 bits (≈ IEEE 754 binary128 / quad-precision),
  giving about 34 significant decimal digits.
  """

  let precision: ULong
  """
  The output precision in bits.  Converted to bytes via `p_bytes()`.
  The minimum meaningful value is 8 (one base-256 digit ≈ 2.4 decimal digits).
  """

  let rounding: RoundingMode
  """
  The IEEE 754 / MPFR rounding mode applied at every operation boundary.
  Defaults to `RoundingNearest` (round to nearest, ties to even).
  """


  new val create(prec: ULong = 112, rnd: RoundingMode = RoundingNearest) =>
    """
    Create an `MPFContext` with the given precision `prec` (bits) and rounding
    mode `rnd`.  The default (no arguments) produces a 112-bit context with
    round-to-nearest, which matches IEEE 754 binary128.
    """
    precision = prec
    rounding  = rnd


  fun val p_bytes(): USize =>
    """
    Output precision in base-256 bytes, computed as `⌈precision / 8⌉`.

    This is the length of the `_digits` array in any `MPFRep` produced by an
    operation under this context.
    """
    (precision.usize() + 7) / 8


  fun val working_bytes(op: String): USize =>
    """
    Working precision in bytes for the named operation `op`.

    Equals `p_bytes() + guard`, where `guard` is a per-operation constant that
    gives enough extra digits to keep the rounding error of the full computation
    ≤ 1 ULP at output precision `p_bytes()`.

    | `op`      | guard | Reason                                                  |
    |-----------|-------|---------------------------------------------------------|
    | `"add"`   |   2   | Catastrophic cancellation loses at most 1 ULP           |
    | `"sub"`   |   2   | Same as `"add"`                                         |
    | `"mul"`   |   2   | Product has 2× digits; FFT rounding errors              |
    | `"inv"`   |   2   | Newton error ≤ 1 ULP at convergence                     |
    | `"sqrt"`  |   2   | Newton error ≤ 1 ULP at convergence                     |
    | `"ln"`    |   6   | Argument reduction introduces cancellation              |
    | `"exp"`   |   8   | Argument reduction + Taylor series                      |
    | `"trig"`  |   8   | Near-cancellation for large arguments                   |
    | `"pi"`    |   8   | Borwein / BBP / Chudnovsky convergence                  |
    | default   |   4   | Conservative fallback for unlisted operations           |
    """
    let p = p_bytes()
    let guard: USize =
      if   (op == "add") or (op == "sub") then 2
      elseif (op == "mul") then 2
      elseif (op == "inv") then 2
      elseif (op == "sqrt") then 2
      elseif (op == "ln") then 6
      elseif (op == "exp") then 8
      elseif (op == "trig") then 8
      elseif (op == "pi") then 8
      else 4
      end
    p + guard


  fun val _round_to(r: MPFRep, n: USize, mode: RoundingMode): MPFRep =>
    """
    Round `r` to `n` most-significant base-256 digits using `mode`.

    This is the keystone operation that gives semantic meaning to
    `rounding`: it replaces every `._trunc(n)` call at operation output
    boundaries.  Inside Newton loops `._trunc` (= `RoundingZero`) remains
    correct and intentional.

    ## Algorithm

    Let `p = r._size()`.  If `p ≤ n` the value already fits; return a copy
    with `_digits` zero-padded to `n` (no rounding needed).

    Otherwise, let `d[n]` be the first discarded byte (index `n` in
    `_digits`), and let *sticky* be `true` when any byte at index `> n` is
    nonzero.

    ### Mode behaviour when discarded bytes are nonzero

    | Mode              | Action                                                     |
    |-------------------|------------------------------------------------------------|
    | `RoundingZero`    | Truncate: keep top `n` bytes unchanged                     |
    | `RoundingNearest` | Round-half-to-even: increment if `d[n] > 128`, or if      |
    |                   | `d[n] == 128` and (sticky or `d[n-1]` is odd)             |
    | `RoundingNegInf`  | If negative, increment magnitude; if positive, truncate    |
    | `RoundingPosInf`  | If positive, increment magnitude; if negative, truncate    |
    | `RoundingAwayZ`   | Always increment magnitude if any discarded byte nonzero  |
    | `RoundingFaithful`| Same as `RoundingNearest` (implementation-defined)        |

    "Increment magnitude" means adding 1 to the last kept digit and
    propagating the carry, which may shift the exponent by 1 when all `n`
    digits overflow (e.g. `[255,255,…]` → `[1,0,0,…]` with exponent + 1).
    """
    let p = r._size()

    // Already fits: zero-pad if needed, no rounding.
    if p <= n then
      if p == n then
        return r
      end
      let padded: Array[U8] val = recover
        let d = Array[U8].create(n)
        var i: USize = 0
        while i < p do
          try d.push(r.raw_digits()(i)?) end
          i = i + 1
        end
        while d.size() < n do
          d.push(0)
        end
        d
      end
      return MPFRep._create(r.sign_bit(), r.is_nan(), r.is_infinite(),
        r.exponent(), padded)
    end

    // Gather rounding information.
    let digits: Array[U8] val = r.raw_digits()
    let first_discarded: U8 = try digits(n)? else 0 end
    var sticky: Bool = false
    var j: USize = n + 1
    while j < p do
      if (try digits(j)? else 0 end) != 0 then
        sticky = true
        break
      end
      j = j + 1
    end

    let last_kept: U8 = try digits(n - 1)? else 0 end

    // Decide whether to increment the magnitude.
    let increment: Bool =
      match mode
      | RoundingZero    => false
      | RoundingNearest =>
        if first_discarded > 128 then
          true
        elseif first_discarded == 128 then
          sticky or ((last_kept and 1) == 1)
        else
          false
        end
      | RoundingFaithful =>
        if first_discarded > 128 then
          true
        elseif first_discarded == 128 then
          sticky or ((last_kept and 1) == 1)
        else
          false
        end
      | RoundingNegInf  => r.sign_bit() and ((first_discarded != 0) or sticky)
      | RoundingPosInf  => (not r.sign_bit()) and ((first_discarded != 0) or sticky)
      | RoundingAwayZ   => (first_discarded != 0) or sticky
      end

    // Truncate to n digits.
    var truncated: MPFRep = r._trunc(n)

    if not increment then
      return truncated
    end

    // Increment magnitude: add 1 to the last kept digit with carry propagation.
    // We work on a mutable copy of the truncated digits.
    let new_exp: I64 = truncated.exponent()
    let new_digits: Array[U8] val = truncated.raw_digits()

    // Check if we overflow (all n digits are 255).
    var all_max: Bool = true
    var k: USize = 0
    while k < n do
      if (try new_digits(k)? else 0 end) != 255 then
        all_max = false
        break
      end
      k = k + 1
    end

    if all_max then
      // Overflow: result is 1.000…0 × 256^(new_exp+1).
      let overflow_digits: Array[U8] val = recover
        let d = Array[U8].init(0, n)
        try d(0)? = 1 end
        d
      end
      MPFRep._create(truncated.sign_bit(), false, false,
        new_exp + 1, overflow_digits)
    else
      // Normal carry propagation from the last digit toward the first.
      let incremented: Array[U8] val = recover
        let d = Array[U8].create(n)
        var i: USize = 0
        while i < n do
          try d.push(new_digits(i)?) end
          i = i + 1
        end
        var pos: USize = n - 1
        var carry: Bool = true
        while carry do
          let cur: U8 = try d(pos)? else 0 end
          if cur == 255 then
            try d(pos)? = 0 end
          else
            try d(pos)? = cur + 1 end
            carry = false
          end
          if (pos == 0) and carry then
            // Should not reach here (all_max handled above), but be safe.
            carry = false
          elseif carry then
            pos = pos - 1
          end
        end
        d
      end
      MPFRep._create(truncated.sign_bit(), false, false,
        new_exp, incremented)
    end
