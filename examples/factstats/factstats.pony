// Digits repartition of factorials up to 1000

use "../../mathx"
use "collections"
use "format"

actor Main
  """
  Analyzes the distribution of digits in decimal representation of factorials.
  Uses `MPInt` for arbitrary-precision factorial calculation.
  """

  new create(env: Env) =>
    """
    Analyzes factorial distributions from 0 to 1000 and prints results in a table.
    """
    let out = env.out
    
    // Header for the table
    print_header(out)
    
    var fact = MPInt.from[ILong](ILong.from[U64](1))
    // total_counts(0..9) for digits, total_counts(10) for total number of digits
    let total_counts = Array[U64].init(0, 11)
    
    // N = 0
    process_fact(0, fact, total_counts, out)
    
    // Loop through N = 1 to 1000
    for n in Range[U64](1, 1001) do
      fact = fact * MPInt.from[ILong](ILong.from[U64](n))
      process_fact(n, fact, total_counts, out)
      
      // Periodic summaries every 50 lines
      if (n > 0) and ((n % 50) == 0) then
        print_recap(total_counts, out)
        if n < 1000 then print_header(out) end
      end
    end


  fun print_header(out: OutStream) =>
    """
    Prints the table header with column names.
    """
    // Fixed widths: N (4), space (1), fact(N) (124), space (1), Digits 0-9 (10 * 8 = 80), space (1), Total Digits (8)
    var h = pad("N", 4, false) + " " + 
            pad("fact(N) (up to 120 msb)", 124) + " "
    
    for i in Range[USize](0, 10) do
      h = (consume h) + pad(i.string(), 8, false)
    end
    
    h = (consume h) + " " + pad("Digits", 8, false)
    out.print(consume h)


  fun print_recap(total_counts: Array[U64], out: OutStream) =>
    """
    Prints cumulative totals (SUM) and relative occurrences (%) of each digit.
    """
    // Cumulative Sums
    // Offset before digit columns: 4 (N) + 1 (space) + 124 (fact) + 1 (space) = 130
    var sum_line = pad("SUM", 130)
    for i in Range(0, 10) do
      let s_val = try total_counts(i)?.string() else "0" end
      sum_line = (consume sum_line) + pad(s_val, 8, false)
    end
    // Add sum of total digits
    let total_s = try total_counts(10)?.string() else "0" end
    sum_line = (consume sum_line) + " " + pad(total_s, 8, false)
    out.print(consume sum_line)
    
    // Total digits encountered so far (sum of counts for digits 0-9)
    var total_digits_sum: U64 = 0
    for i in Range(0, 10) do
      total_digits_sum = total_digits_sum + (try total_counts(i)? else 0 end)
    end
    
    // Relative Percentages (2 decimals)
    var rel_line = pad("%", 130)
    for i in Range(0, 10) do
      let r_val = try
        let rel = (total_counts(i)?.f64() / total_digits_sum.f64()) * 100.0
        Format.float[F64](rel where fmt = FormatFix, prec = 2)
      else
        "0.00"
      end
      rel_line = (consume rel_line) + pad(r_val, 8, false)
    end
    // Relative percentage for the total column is 100.00%
    rel_line = (consume rel_line) + " " + pad("100.00", 8, false)
    out.print((consume rel_line) + "\n")


  fun process_fact(n: U64, fact: MPInt, total_counts: Array[U64], out: OutStream) =>
    """
    Calculates digit counts for a single factorial, updates globals and prints the data line.
    """
    let s: String = fact.string()
    let digit_count = s.size()
    
    // Update total digit count sum
    try
      total_counts(10)? = total_counts(10)? + digit_count.u64()
    end
    
    // Digits analysis (0-9)
    let counts = Array[U64].init(0, 10)
    for c in s.values() do
      if (c >= '0') and (c <= '9') then
        let digit = (c - '0').usize()
        try
          counts(digit)? = counts(digit)? + 1
          total_counts(digit)? = total_counts(digit)? + 1
        end
      end
    end
    
    // Display string: keep only the 120 most significant digits and append '...' if needed
    let display_s: String = if s.size() > 120 then
      s.substring(0, 120) + "..."
    else
      s
    end
    
    // Assemble and print the line
    // Use 124 for the fact(N) column to accommodate 120 digits + "..."
    var line = pad(n.string(), 4, false) + " " + pad(display_s, 124) + " "
    for i in Range(0, 10) do
      let c_val = try counts(i)?.string() else "0" end
      line = (consume line) + pad(c_val, 8, false)
    end
    line = (consume line) + " " + pad(digit_count.string(), 8, false)
    out.print(consume line)


  fun pad(s: String, width: USize, left: Bool = true): String =>
    """
    Right or left-pads a string with spaces to reach the specified width.
    """
    if s.size() >= width then
      s
    else
      let p_size = width - s.size()
      let p: String = recover
        let s' = String(p_size)
        var k: USize = 0
        while k < p_size do
          s'.push(' ')
          k = k + 1
        end
        s'
      end
      if left then
        s + p
      else
        p + s
      end
    end
