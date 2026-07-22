// A library of various mathematical functions usable by MPFRep


primitive MPFMathLibRep
    """
    The `MPFMathLibRep` mathematical library contains various mathematical
    algorithms that deal with `MPFRep` instances. If you want to use more
    friendly `MPFloat` objects, use the `MPFMathLib`  primitive instead.
    """

  fun pi_machin(ctx: MPFContext): MPFRep val =>
    """
    Compute π at output precision `ctx.p_bytes()` using
    [Machin formula](https://en.wikipedia.org/wiki/Machin-like_formula):
    `π = 16·arctan(1/5) − 4·arctan(1/239)`.

    Uses Taylor `arctan(x) = x − x³/3 + x⁵/5 − ...`

    **Rounding**: faithfully rounded (≤ 1 ULP error).
    """
    let w = ctx.working_bytes("pi")
    let p: USize = w + 4
    let one: MPFRep val = MPFRep.from[F64](1.0, p)
    let five: MPFRep val = MPFRep.from[F64](5.0, p)
    let two39: MPFRep val = MPFRep.from[F64](239.0, p)
    let x5: MPFRep val = _MPFAlgo._inv(five, p)
    let x239: MPFRep val = _MPFAlgo._inv(two39, p)

    // arctan(1/5) via Taylor series
    let neg_x52: MPFRep val = _MPFAlgo._neg_rep(_MPFAlgo._mul(x5, x5, p))
    var pow5: MPFRep val = x5
    var atan5: MPFRep val = x5
    var k5: USize = 1
    var iters5: USize = 0
    let max_iters: USize = p * 30
    while iters5 < max_iters do
      pow5 = _MPFAlgo._mul(pow5, neg_x52, p)
      let d5: USize = (2 * k5) + 1
      let term5: MPFRep val = if d5 <= 255 then
          (let q5: MPFRep val, _) = pow5._short_div(d5.u8())
          q5
        else
          _MPFAlgo._div(pow5, MPFRep.from[F64](d5.f64(), p), p)
        end
      let new_a5: MPFRep val = _MPFAlgo._add(atan5, term5, p)

      if _MPFAlgo._converged(new_a5, atan5) then
        break
      end

      atan5 = new_a5
      k5 = k5 + 1
      iters5 = iters5 + 1
    end

    // arctan(1/239) via Taylor series
    let neg_x2392: MPFRep val = _MPFAlgo._neg_rep(_MPFAlgo._mul(x239, x239, p))
    var pow239: MPFRep val = x239
    var atan239: MPFRep val = x239
    var k239: USize = 1
    var iters239: USize = 0
    while iters239 < max_iters do
      pow239 = _MPFAlgo._mul(pow239, neg_x2392, p)
      let d239: USize = (2 * k239) + 1
      let term239: MPFRep val = if d239 <= 255 then
          (let q239: MPFRep val, _) = pow239._short_div(d239.u8())
          q239
        else
          _MPFAlgo._div(pow239, MPFRep.from[F64](d239.f64(), p), p)
        end
      let new_a239: MPFRep val = _MPFAlgo._add(atan239, term239, p)

      if _MPFAlgo._converged(new_a239, atan239) then
        break
      end

      atan239 = new_a239
      k239 = k239 + 1
      iters239 = iters239 + 1
    end

    // π = 16·arctan(1/5) − 4·arctan(1/239)
    let sixteen: MPFRep val = MPFRep.from[F64](16.0, p)
    let four: MPFRep val = MPFRep.from[F64](4.0, p)
    let result = _MPFAlgo._sub(_MPFAlgo._mul(sixteen, atan5, p), _MPFAlgo._mul(four, atan239, p), p)._trunc(w)
    ctx._round_to(result, ctx.p_bytes(), ctx.rounding)


  fun pi_bbp(ctx: MPFContext): MPFRep val =>
    """
    Compute π using the Bailey–Borwein–Plouffe formula at output precision
    `ctx.p_bytes()`.
    [BBP formula](https://en.wikipedia.org/wiki/Bailey%E2%80%93Borwein%E2%80%93Plouffe_formula):
    `π = Σ_{k=0}^∞ (1/16^k) × [4/(8k+1) − 2/(8k+4) − 1/(8k+5) − 1/(8k+6)]`.
 
    **Rounding**: faithfully rounded (≤ 1 ULP error).

    The `1/16^k` factor is maintained as a running product to avoid `powi`.
     """
    let w = ctx.working_bytes("pi")
    let p: USize = w + 4
    let k_1: MPFRep val = MPFRep.from[F64](1.0, p)
    let k_2: MPFRep val = MPFRep.from[F64](2.0, p)
    let k_4: MPFRep val = MPFRep.from[F64](4.0, p)
    let k_5: MPFRep val = MPFRep.from[F64](5.0, p)
    let k_6: MPFRep val = MPFRep.from[F64](6.0, p)
    let k_8: MPFRep val = MPFRep.from[F64](8.0, p)
    let inv16: MPFRep val = _MPFAlgo._inv(MPFRep.from[F64](16.0, p), p)
    var k: USize = 0
    var result: MPFRep val = MPFRep._create(false, false, false, 0, Array[U8].init(0, p))
    var prev_res: MPFRep val = result
    var pow16k: MPFRep val = k_1
    repeat
      let t0: MPFRep val = _MPFAlgo._mul(k_8, MPFRep.from[ULong](k.ulong(), p), p)
      let t1: MPFRep val = _MPFAlgo._div(k_4, _MPFAlgo._add(t0, k_1, p), p)
      let t2: MPFRep val = _MPFAlgo._div(k_2, _MPFAlgo._add(t0, k_4, p), p)
      let t3: MPFRep val = _MPFAlgo._div(k_1, _MPFAlgo._add(t0, k_5, p), p)
      let t4: MPFRep val = _MPFAlgo._div(k_1, _MPFAlgo._add(t0, k_6, p), p)
      let inner: MPFRep val = _MPFAlgo._sub(_MPFAlgo._sub(_MPFAlgo._sub(t1, t2, p), t3, p), t4, p)
      let term: MPFRep val = _MPFAlgo._mul(pow16k, inner, p)
      prev_res = result
      result = _MPFAlgo._add(result, term, p)
      pow16k = _MPFAlgo._mul(pow16k, inv16, p)
      k = k + 1
    // BBP delivers log2(16)=4 bits per term, so p bytes=8p bits requires ~2p
    // terms. Safety limit is 2*p+8.
    until _MPFAlgo._converged(result, prev_res) or (k > ((2 * p) + 8)) end
    
    ctx._round_to(result._trunc(w), ctx.p_bytes(), ctx.rounding)

