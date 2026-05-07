// A library of mathematical functions usable by MPFloat

primitive MPFMathLib
    """
    A library of mathematical functions dealing with `MPFloat` objects.
    """
  fun pi_machin(prec: USize = 128, rnd: RoundingMode = RoundingNearest): MPFloat =>
    """
    The pi constant, calculated with the specified `prec` accuracy (number of
    bits, default 128 giving ~38 decimal digits).

    Uses Machin's formula:
      π = 16·arctan(1/5) − 4·arctan(1/239)
    where each arctan is computed by the Taylor series
      arctan(x) = x − x³/3 + x⁵/5 − x⁷/7 + ⋯
    Both arguments are small (< 0.25), giving fast geometric convergence.
    """
    let ctx = MPFContext(prec, rnd)
    let rep = MPFMathLibRep.pi_machin(ctx)
    MPFloat._from(rep, ctx)


  fun pi_bbp(prec: USize = 128, rnd: RoundingMode = RoundingNearest): MPFloat =>
    """
    Compute π using the Bailey–Borwein–Plouffe (BBP) formula:
    `π = Σ_{k=0}^∞ (1/16^k) × [4/(8k+1) − 2/(8k+4) − 1/(8k+5) − 1/(8k+6)]`.
    """
    let ctx = MPFContext(prec, rnd)
    let rep = MPFMathLibRep.pi_bbp(ctx)
    MPFloat._from(rep, ctx)


