// Multi-precision float context: precision policy and rounding only.

use "../mathx"
use "collections"


class val MPFContext
  """
  Precision and rounding policy for arbitrary-precision floating-point arithmetic.

  An `MPFContext` is a pure policy record: it holds two orthogonal pieces of
  policy and provides helpers to derive working precision and apply output
  rounding. It has **no arithmetic methods** — all computation lives in
  `_MPFAlgo`.

  - **`precision`** — the number of significant bits in any result produced
    under this context. Stored in bits (user-facing unit); converted to bytes
    internally via `p_bytes()`.
  - **`rounding`** — the IEEE 754 / MPFR rounding mode applied when a
    computed result is shortened to `precision` bits.

  ## Usage pattern

  `MPFloat` holds an `MPFContext` and calls `_MPFAlgo.op(_ctx, ...)` directly.
  Every operation follows this two-phase discipline:

  1. Compute at *working precision* `w = working_bytes(op)` bytes — enough
     guard bytes to bound the rounding error of the operation to ≤ 1 ULP.
  2. Round back to *output precision* `p = p_bytes()` bytes via `_round_to`.

  ## Default precision

  The default precision is 128 bits, giving about 38 significant decimal
  digits. This covers the full range of `U128`/`I128` integer values without
  loss, which is why 128 bits is preferred over the 112-bit IEEE 754 binary128
  mantissa width.
  """

  let precision: USize
  """
  The output precision in bits. Converted to bytes via `p_bytes()`.
  The minimum meaningful value is 8 (one base-256 digit ≈ 2.4 decimal digits).
  """

  let rounding: RoundingMode
  """
  The IEEE 754 / MPFR rounding mode applied at every operation boundary.
  Defaults to `RoundingNearest` (round to nearest, ties to even).
  """


  new val create(prec: USize = 128, rnd: RoundingMode = RoundingNearest) =>
    """
    Create an `MPFContext` with the given precision `prec` (bits) and rounding
    mode `rnd`. The default (no arguments) produces a 128-bit context with
    round-to-nearest, giving ~38 significant decimal digits.
    """
    precision = prec
    rounding  = rnd


  fun val p_bytes(): USize =>
    """
    Output precision in base-256 bytes, computed as `⌈precision / 8⌉`.

    This is the length of the `_digits` array in any `MPFRep` produced by an
    operation under this context.
    """
    (precision + 7) / 8


  fun val working_bytes(op: String): USize =>
    """
    Working precision in bytes for the named operation `op`.

    Delegates to `_MPFAlgo._guard_bytes`, which owns the authoritative
    guard-byte table next to the algorithms that need it. See that method
    for the per-operation guard values and their rationale.
    """
    _MPFAlgo._guard_bytes(this, op)


  fun val _round_to(rep: MPFRep box, prec: USize, mode: RoundingMode): MPFRep val =>
    """
    Round `rep` to `prec` most-significant base-256 digits using `mode`.

    ## Algorithm

    Let `p = rep._size()`. If `p ≤ prec` the value already fits; return a copy
    with `_digits` zero-padded to `prec` (no rounding needed).

    Otherwise, let `d[prec]` be the first discarded byte (index `prec` in
    `_digits`), and let *sticky* be `true` when any byte at index `> prec` is
    nonzero (See https://en.wikipedia.org/wiki/Floating-point_arithmetic#Addition_and_subtraction
    for information on rounding and sticky bit)

    ### Mode behaviour when discarded bytes are nonzero

    | Mode              | Action                                                   |
    |-------------------|----------------------------------------------------------|
    | `RoundingZero`    | Truncate: keep top `prec` bytes unchanged                |
    | `RoundingNearest` | Round-half-to-even: increment if `d[prec] > 128`, or if  |
    |                   | `d[prec] == 128` and (sticky or `d[prec - 1]` is odd)    |
    | `RoundingNegInf`  | If negative, increment magnitude; if positive, truncate  |
    | `RoundingPosInf`  | If positive, increment magnitude; if negative, truncate  |
    | `RoundingAwayZ`   | Always increment magnitude if any discarded byte nonzero |
    | `RoundingFaithful`| Same as `RoundingNearest` (implementation-defined)       |

    "Increment magnitude" means adding 1 to the last kept digit and propagating
    the carry, which may shift the exponent by 1 when all `prec` digits overflow
    (e.g. `[255, 255, ...]` → `[1, 0, 0, ...]` with exponent + 1).
    """
    let p = rep._size()
    let midb: U8 = 128
    let base_max: U8 = 255

    // Already fits: zero-pad if needed, no rounding.
    if p <= prec then
      let src = rep.raw_digits()
      let padded: Array[U8] val = recover
        let d = Array[U8].create(prec)
        var i: USize = 0
        while i < p do
          try d.push(src(i)?) end
          i = i + 1
        end
        while d.size() < prec do
          d.push(0)
        end
        d
      end
      return MPFRep._create(rep.sign_bit(), rep.is_nan(), rep.is_infinite(),
        rep.exponent(), padded)
    end

    // Gather rounding information.
    let digits: Array[U8] val = rep.raw_digits()
    let first_discarded: U8 = try digits(prec)? else 0 end
    var sticky: Bool = false
    for j in Range(prec + 1, p) do
      if (try digits(j)? else 0 end) != 0 then
        sticky = true
        break
      end
    end

    let last_kept: U8 = try digits(prec - 1)? else 0 end

    // Decide whether to increment the magnitude.
    let increment: Bool = match \exhausitve\ mode
      | RoundingZero => false
      | RoundingNearest => if first_discarded > midb then
          true
        elseif first_discarded == midb then
          sticky or ((last_kept and 1) == 1)
        else
          false
        end
      | RoundingFaithful => if first_discarded > midb then
          true
        elseif first_discarded == midb then
          sticky or ((last_kept and 1) == 1)
        else
          false
        end
      | RoundingNegInf  => rep.sign_bit() and ((first_discarded != 0) or sticky)
      | RoundingPosInf  => (not rep.sign_bit()) and ((first_discarded != 0) or sticky)
      | RoundingAwayZ   => (first_discarded != 0) or sticky
      end

    // Truncate to prec digits.
    let truncated: MPFRep val = rep._trunc(prec)

    if not increment then
      return truncated
    end

    // Increment magnitude: add 1 to the last kept digit with carry propagation.
    // We work on a mutable copy of the truncated digits.
    let new_exp: I64 = truncated.exponent()
    let new_digits: Array[U8] val = truncated.raw_digits()

    // Check if we overflow (all prec digits are 255).
    var all_max: Bool = true
    for k in Range(0, prec) do
      if (try new_digits(k)? else 0 end) != base_max then
        all_max = false
        break
      end
    end

    if all_max then
      // Overflow: result is 1.000...0 × 256^(new_exp+1).
      let overflow_digits: Array[U8] val = recover
        let d = Array[U8].init(0, prec)
        try d(0)? = 1 end
        d
      end
      MPFRep._create(truncated.sign_bit(), false, false, new_exp + 1, overflow_digits)
    else
      // Normal carry propagation from the last digit toward the first.
      let incremented: Array[U8] val = recover
        let d = Array[U8].create(prec)
        for i in Range(0, prec) do
          try d.push(new_digits(i)?) end
        end
        var pos: USize = prec - 1
        var carry: Bool = true
        while carry do
          let cur: U8 = try d(pos)? else 0 end
          if cur == base_max then
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
      MPFRep._create(truncated.sign_bit(), false, false, new_exp, incremented)
    end

