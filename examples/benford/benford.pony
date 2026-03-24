/*
  Benford's Law verification example.
  Checks the frequency of leading digits in various mathematical sequences.
*/

use "../../mathx"
use "collections"
use "format"
use "random"
use "time"

actor Main
  """
  Main entry point for Benford's Law verification.
  Calculates frequencies for Fibonacci, Factorials, and Linear Recurrence series.
  """
  // Number of terms to calculate for each series
  let _terms: USize = 1_000

  new create(env: Env) =>
    """
    Run the Benford's Law verification on all requested series.
    """
    env.out.print("Benford's Law Verification (" + _terms.string() + " terms per series)")
    env.out.print("================================================================")

    // 1. Fibonacci
    run_fibonacci(env)

    // 2. Factorials
    run_factorials(env)

    // 3. Linear Recurrence x_n+1 = 2*x_n + 1
    // 3a. With 5 random numbers
    run_random_recurrence(env)

    // 3b. With specific x_0 = 9.94962308959395941219332124109226
    //run_specific_recurrence(env)


  fun run_fibonacci(env: Env) =>
    """
    Calculate Fibonacci numbers and analyze their digits.
    """
    env.out.print("\nSerie: Fibonacci Numbers (F_n = F_{n-1} + F_{n-2})")
    let collector = StatCollector(5)
    
    var f0 = MPInt.from_ulong(ULong.from[U64](0))
    var f1 = MPInt.from_ulong(ULong.from[U64](1))
    
    // F_1 is 1
    collector.add_number(f1)
    
    for i in Range[USize](1, _terms) do
      let f2 = f0 + f1
      f0 = f1
      f1 = f2
      collector.add_number(f1)
    end
    
    collector.print_table(env.out)


  fun run_factorials(env: Env) =>
    """
    Calculate factorials and analyze their digits.
    """
    env.out.print("\nSerie: Factorials (n!)")
    let collector = StatCollector(5)
    
    var fact = MPInt.from_ulong(ULong.from[U64](1))
    collector.add_number(fact)
    
    for i in Range[U64](2, _terms.u64() + 1) do
      fact = fact * MPInt.from_ulong(ULong.from[U64](i))
      collector.add_number(fact)
    end
    
    collector.print_table(env.out)


  fun run_random_recurrence(env: Env) =>
    """
    Run the x_n+1 = 2*x_n + 1 series with 5 random starting values.
    """
    let mt = MT(Time.nanos())
    
    for i in Range[USize](0, 5) do
      let x0_u64 = mt.next()
      env.out.print("\nSerie: x_{n+1} = 2*x_n + 1, x_0 = " + x0_u64.string())
      let collector = StatCollector(5)
      
      var x = MPInt.from_ulong(ULong.from[U64](x0_u64))
      collector.add_number(x)
      
      let two = MPInt.from_ulong(ULong.from[U64](2))
      let one = MPInt.from_ulong(ULong.from[U64](1))
      
      for _ in Range[USize](1, _terms) do
        x = (x * two) + one
        collector.add_number(x)
      end
      
      collector.print_table(env.out)
    end


  fun run_specific_recurrence(env: Env) =>
    """
    Run the x_n+1 = 2*x_n + 1 series with x_0 = 9.94962308959395941219332124109226.
    """
    let x0_str = "9.94962308959395941219332124109226"
    env.out.print("\nSerie: x_{n+1} = 2*x_n + 1, x_0 = " + x0_str)
    
    try
      // Use high precision for the series (128 bytes ≈ 308 decimal digits)
      let prec: USize = 128
      var x = MPFloat.from_string(x0_str, prec)?
      let collector = StatCollector(5)
      
      collector.add_number(x)
      
      let two = MPFloat.from_f64(2.0, prec)
      let one = MPFloat.from_f64(1.0, prec)
      
      for _ in Range[USize](1, _terms) do
        x = x.mul(two).add(one)
        collector.add_number(x)
      end
      
      collector.print_table(env.out)
    else
      env.out.print("Error parsing x_0: " + x0_str)
    end


class StatCollector
  """
  Collects and prints digit frequencies for Benford's Law analysis.
  """
  let _max_pos: USize
  let _counts: Array[Array[U64]]
  var _total: U64 = 0

  new create(max_pos: USize) =>
    """
    Initialize a collector for a given number of digit positions.
    """
    _max_pos = max_pos
    _counts = Array[Array[U64]](_max_pos)
    for i in Range[USize](0, _max_pos) do
      _counts.push(Array[U64].init(0, 10))
    end


  fun ref add_number(n: (MPInt | MPFloat)) =>
    """
    Extract the first few digits of the number and update frequencies.
    """
    _total = _total + 1
    
    let s: String = match n
    | let mi: MPInt =>
      // Convert to MPFloat for easy extraction of leading digits
      MPFloat.from_mpint(mi, 10).string()
    | let mf: MPFloat =>
      mf.string()
    end

    // Extract digits, ignoring sign, decimal point, leading zeros, and 'e'
    var digits = String(_max_pos)
    var found_first = false
    
    for c in s.values() do
      if (c >= '1') and (c <= '9') then
        found_first = true
        digits.push(c)
      elseif (c == '0') and found_first then
        digits.push(c)
      elseif (c == 'e') or (c == 'E') then
        break
      end
      if digits.size() == _max_pos then break end
    end
    
    for i in Range[USize](0, digits.size()) do
      try
        let d_char = digits(i)?
        let d = (d_char - '0').usize()
        let row = _counts(i)?
        row(d)? = row(d)? + 1
      end
    end


  fun print_table(out: OutStream) =>
    """
    Print a table comparing observed frequencies with Benford's Law.
    """
    let benford = BenfordLaw
    let col_width: USize = 17 // " 100.00 / 100.00 "
    let digit_width: USize = 7 // " Digit "

    // Top border: ┌───────┬─────────────────┬...┐
    let top = recover val
      let t = String
      t.append("┌")
      for i in Range(0, digit_width) do t.append("─") end
      for pos in Range(0, _max_pos) do
        t.append("┬")
        for i in Range(0, col_width) do t.append("─") end
      end
      t.append("┐")
      t
    end
    out.print(top)

    // Header: │ Digit │      Pos 1      │...│
    let header = recover val
      let h = String
      h.append("│ Digit │")
      for pos in Range[USize](1, _max_pos + 1) do
        let p_str = "Pos " + pos.string()
        let space = col_width - p_str.size()
        let left = space / 2
        let right = space - left
        for i in Range(0, left) do h.append(" ") end
        h.append(consume p_str)
        for i in Range(0, right) do h.append(" ") end
        h.append("│")
      end
      h
    end
    out.print(header)

    // Separator: ├───────┼─────────────────┼...┤
    let mid = recover val
      let m = String
      m.append("├")
      for i in Range(0, digit_width) do m.append("─") end
      for pos in Range(0, _max_pos) do
        m.append("┼")
        for i in Range(0, col_width) do m.append("─") end
      end
      m.append("┤")
      m
    end
    out.print(mid)

    // Body rows
    for digit in Range[USize](0, 10) do
      var line: String = "│   " + digit.string() + "   │"
      for pos in Range[USize](0, _max_pos) do
        let count = try _counts(pos)?(digit)? else 0 end
        let observed: F64 = if _total > 0 then (count.f64() / _total.f64()) * 100.0 else 0.0 end
        let expected: F64 = benford.expected(pos + 1, digit) * 100.0
        
        let obs_str = Format.float[F64](observed where fmt = FormatFix, prec = 2, width = 6)
        let exp_str = Format.float[F64](expected where fmt = FormatFix, prec = 2, width = 6)
        
        let cell = (consume obs_str) + " / " + (consume exp_str)
        let space = col_width - cell.size()
        let left = space / 2
        let right = space - left
        
        var cell_padded = String(col_width)
        for i in Range(0, left) do cell_padded.push(' ') end
        cell_padded.append(consume cell)
        for i in Range(0, right) do cell_padded.push(' ') end
        
        line = (consume line) + (consume cell_padded) + "│"
      end
      out.print(consume line)
    end
    
    // Bottom border: └───────┴─────────────────┴...┘
    let bottom = recover val
      let b = String
      b.append("└")
      for i in Range(0, digit_width) do b.append("─") end
      for pos in Range(0, _max_pos) do
        b.append("┴")
        for i in Range(0, col_width) do b.append("─") end
      end
      b.append("┘")
      b
    end
    out.print(bottom)
    out.print("Format: Observed % / Expected %")


primitive BenfordLaw
  """
  Provides theoretical Benford's Law probabilities.
  """
  fun expected(pos: USize, digit: USize): F64 =>
    """
    Returns the expected probability of 'digit' at 'pos' (1-indexed).
    """
    if digit > 9 then return 0.0 end
    if (pos == 1) and (digit == 0) then return 0.0 end
    
    if pos == 1 then
      // P(d) = log10(1 + 1/d)
      (1.0 + (1.0 / digit.f64())).log10()
    elseif pos == 2 then
      // P(d) = sum_{k=1..9} log10(1 + 1/(10k + d))
      var sum: F64 = 0.0
      for k in Range[USize](1, 10) do
        sum = sum + (1.0 + (1.0 / ((10.0 * k.f64()) + digit.f64()))).log10()
      end
      sum
    elseif pos == 3 then
      // P(d) = sum_{k=10..99} log10(1 + 1/(10k + d))
      var sum: F64 = 0.0
      for k in Range[USize](10, 100) do
        sum = sum + (1.0 + (1.0 / ((10.0 * k.f64()) + digit.f64()))).log10()
      end
      sum
    elseif pos == 4 then
      var sum: F64 = 0.0
      for k in Range[USize](100, 1000) do
        sum = sum + (1.0 + (1.0 / ((10.0 * k.f64()) + digit.f64()))).log10()
      end
      sum
    elseif pos == 5 then
      var sum: F64 = 0.0
      for k in Range[USize](1000, 10000) do
        sum = sum + (1.0 + (1.0 / ((10.0 * k.f64()) + digit.f64()))).log10()
      end
      sum
    else
      // For pos > 1, the distribution quickly approaches 10% for each digit.
      0.1
    end
