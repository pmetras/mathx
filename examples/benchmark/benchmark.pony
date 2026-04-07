// MPInt benchmark suite

use "../../mathx"
use "random"
use "time"


actor Main
  """
  Benchmark suite for `MPInt` operations.

  Each run emits one CSV row per `(benchmark, size)` pair to stdout.
  Redirect or append the output to a persistent log file, then analyse it
  with a spreadsheet or with Python / pandas:

  ```sh
  # Initialise a new log (write the header once)
  ./build/release/benchmark --header-only > bench.csv

  # Append a release run
  ./build/release/benchmark --no-header --iter=10 >> bench.csv

  # Append a debug baseline at a subset of sizes
  ./build/debug/benchmark --no-header --bench=mul,mul_fft --size=100,500 >> bench.csv

  # Analyse in Python
  # import pandas as pd
  # df = pd.read_csv('bench.csv')
  # df['datetime'] = pd.to_datetime(df['ts_s'], unit='s')
  # df.pivot_table(values='ns_per_iter', index='n_digits',
  #                columns=['build','bench']).plot(logy=True)
  ```

  ## Output columns

  | Column       | Description                                              |
  |--------------|----------------------------------------------------------|
  | `ts_s`       | Unix timestamp (seconds) of the run start               |
  | `build`      | `debug` or `release` (compiled in at build time)        |
  | `bench`      | Benchmark name (see below)                              |
  | `n_digits`   | Operand size in base-65536 digits                       |
  | `n_iter`     | Number of timed iterations                              |
  | `elapsed_ns` | Total wall-clock nanoseconds for all `n_iter` iterations|
  | `ns_per_iter`| `elapsed_ns / n_iter` — per-operation cost              |
  | `checksum`   | XOR of result bytes; detects silent correctness changes |
  | `comment`    | Free-text note from `--comment=TEXT`; empty if not given |

  ## Benchmarks

  | Name           | Operation           | Complexity  | Notes                       |
  |----------------|---------------------|-------------|-----------------------------|
  | `mul`          | `a * b`             | O(n²) now   | Tracks dispatch improvement |
  | `mul_karatsuba`| `a.mul_karatsuba(b)`| O(n^1.585)  | Karatsuba algorithm         |
  | `mul_fft`      | `a.mul_fft(b)`      | O(n log n)  | FFT; limited to ~1M digits  |
  | `string`       | `a.string()`        | O(n²) now   | Decimal conversion          |
  | `from_string`  | `MPInt.from_string` | O(n²) now   | Decimal parsing             |
  | `pow2`         | `a * a`             | O(n²) now   | Squaring; proxy for `pow`   |
  | `divrem`       | `a.divrem(b)`       | O(n²) now   | Tracks division improvement |
  | `isqrt`        | `(a*a+1).isqrt()`   | O(n² log n) | Newton; uses divrem loop    |

  ## CLI flags

  | Flag                 | Default              | Meaning                             |
  |----------------------|----------------------|-------------------------------------|
  | `--bench=n1,n2,...`  | all benchmarks       | Comma-separated bench names to run  |
  | `--size=n1,n2,...`   | 10,100,500,2000      | Operand sizes in base-65536 digits  |
  | `--iter=N`           | 5                    | Iterations per (bench, size) pair   |
  | `--max-mem=N`        | 8                    | Memory limit in GB; skip if exceeded|
  | `--comment=TEXT`     | (empty)              | Free-text note recorded in each row |
  | `--no-header`        | off                  | Suppress CSV header line            |
  | `--header-only`      | off                  | Print header and exit               |

  ## Memory guidance (64-bit system)

  Each base-65536 digit occupies 2 bytes.  The dominant allocation per run:

  - `mul_fft`: two FFT arrays of `next_pow2(2n)` F64 values ≈ 32 × n bytes.
    At 8 GB limit this means n ≤ ~250 million digits.
  - All other benchmarks: O(n) bytes; the time constraint is more binding
    than memory (schoolbook is O(n²), skip large sizes).
  - `mul_schoolbook` is auto-skipped for n > 500 to avoid multi-minute runtimes.
  - `isqrt` squares the operand first (2n-digit input), then runs O(log n)
    Newton steps each using `divrem`; it is auto-skipped for n > 200.

  ## Reproducibility

  Operands are always built from fixed seeds (`0xDEADBEEF` / `0xCAFEBABE`),
  so repeated runs with the same `--size` produce comparable checksums.
  Record the `ponyc` version alongside the log; `build` (debug/release) is
  already in every row.
  """


  new create(env: Env) =>
    """
    Parse CLI flags, optionally print the CSV header, then run each selected
    `(benchmark, size)` pair and emit one CSV row per pair.
    """
    let all_bench_names: Array[String] val = [
      "mul"
      "mul_schoolbook"
      "mul_karatsuba"
      "mul_ntt"
      "mul_fft"
      "string"
      "from_string"
      "pow2"
      "divrem"
      "isqrt"
    ]

    var print_header: Bool    = true
    var header_only: Bool     = false
    var all_benches: Bool     = true
    var n_iter:      USize    = 5
    var max_mem_gb:  U64      = 8
    var comment:     String   = ""

    let selected: Array[String] = Array[String]
    var sizes: Array[USize]     = [10; 100; 500; 2000]

    // Parse CLI flags (skip args(0) = binary name)
    let bench_prefix:   String = "--bench="
    let size_prefix:    String = "--size="
    let iter_prefix:    String = "--iter="
    let maxmem_prefix:  String = "--max-mem="
    let comment_prefix: String = "--comment="

    var i: USize = 1
    while i < env.args.size() do
      let arg = try env.args(i)? else "" end
      if arg == "--no-header" then
        print_header = false
      elseif arg == "--header-only" then
        header_only = true
      elseif arg.at(bench_prefix, 0) then
        all_benches = false
        let val_str = arg.substring(bench_prefix.size().isize())
        for name in val_str.split_by(",").values() do
          selected.push(name)
        end
      elseif arg.at(size_prefix, 0) then
        sizes = Array[USize]
        let val_str = arg.substring(size_prefix.size().isize())
        for s in val_str.split_by(",").values() do
          try sizes.push(s.usize()?) end
        end
      elseif arg.at(iter_prefix, 0) then
        try n_iter = arg.substring(iter_prefix.size().isize()).usize()? end
      elseif arg.at(maxmem_prefix, 0) then
        try max_mem_gb = arg.substring(maxmem_prefix.size().isize()).u64()? end
      elseif arg.at(comment_prefix, 0) then
        comment = arg.substring(comment_prefix.size().isize())
      else
        env.err.print("benchmark: unknown flag: " + arg)
        env.err.print("Run with --help to see usage.")
      end
      i = i + 1
    end

    let build_type: String = ifdef debug then "debug" else "release" end
    (let ts_s, _) = Time.now()
    let max_mem_bytes: U64 = max_mem_gb * 1_073_741_824

    if print_header then
      env.out.print(
        "ts_s,build,bench,n_digits,n_iter,elapsed_ns,ns_per_iter,checksum,comment")
    end

    if header_only then
      return
    end

    // Determine the benchmark list to run.
    let benches: Array[String] = if all_benches then
      let b: Array[String] = Array[String]
      for name in all_bench_names.values() do
        b.push(name)
      end
      b
    else
      selected
    end

    // Run each (bench, size) pair.
    for bench in benches.values() do
      for n in sizes.values() do
        // Per-benchmark size and memory guards.
        let skip_msg: String = _skip_reason(bench, n, max_mem_bytes)
        if skip_msg.size() > 0 then
          env.err.print(
            "benchmark: skipping " + bench + " n=" + n.string()
            + " — " + skip_msg)
          continue
        end

        (let elapsed, let checksum) = _run_bench(bench, n, n_iter)
        let ns_per_iter: U64 =
          if n_iter > 0 then elapsed / n_iter.u64() else 0 end

        env.out.print(
          ts_s.string()          + ","
          + build_type           + ","
          + bench                + ","
          + n.string()           + ","
          + n_iter.string()      + ","
          + elapsed.string()     + ","
          + ns_per_iter.string() + ","
          + checksum.string()    + ","
          + _csv_quote(comment))
      end
    end


  fun _skip_reason(bench: String, n: USize, max_mem_bytes: U64): String =>
    """
    Return a non-empty string describing why `(bench, n)` should be skipped,
    or an empty string if it is safe to run.

    Skipping rules:
    - `mul_schoolbook` (alias `mul`) for n > 500: O(n²) time exceeds minutes.
    - `isqrt` for n > 200: squares first (2n digits) then O(n² log n) Newton.
    - Any benchmark where estimated peak memory exceeds `max_mem_bytes`.

    Memory estimates (conservative upper bounds):
    - `mul_fft`: FFT arrays dominate, ≈ 32 × n bytes.
    - All others: operands + temporaries, ≈ 24 × n bytes.
    """
    // Time guards (independent of memory).
    match bench
    | "mul_schoolbook" =>
      // mul_schoolbook is O(n²); cap at the Karatsuba dispatch threshold.
      if n > 500 then
        return "O(n²) schoolbook too slow for n > 500 (use mul_karatsuba or mul_fft)"
      end
    | "isqrt" =>
      if n > 200 then
        return "squares operand first then O(n² log n) Newton; too slow for n > 200"
      end
    end

    // Memory guard.
    // mul_fft allocates two FFT arrays of next_pow2(2n) F64 values ≈ 32n bytes.
    // mul dispatches to mul_fft for n ≥ 512, so use the same factor.
    let bytes_per_digit: U64 = match bench
      | "mul_fft" => 32
      | "mul_ntt" => 32
      | "mul"     => 32
      else 24
    end
    let estimated: U64 = bytes_per_digit * n.u64()
    if estimated > max_mem_bytes then
      return (
        "estimated memory " + (estimated / 1_073_741_824).string()
        + " GB exceeds --max-mem=" + (max_mem_bytes / 1_073_741_824).string()
        + " GB")
    end

    ""


  fun _make_mpint(n: USize, seed: U64): MPInt =>
    """
    Build a reproducible `n`-digit (base-65536) positive `MPInt` from `seed`.

    The most-significant digit is forced non-zero so the number has exactly
    `n` digits.  Construction is O(n²) — it is done once per benchmark run,
    outside the timed hot loop.

    Returns 0 for `n = 0`.
    """
    if n == 0 then
      return MPInt.from_ilong(0)
    end

    var rand = Rand(seed)

    // Top digit: non-zero so the result has exactly n digits.
    let top: ILong = ((rand.u64() % 65535) + 1).ilong()
    var result: MPInt = MPInt.from_ilong(top).digit_shl(n - 1)

    // Fill positions 0 .. n-2 with random digits.
    var k: USize = 0
    while k < (n - 1) do
      let d: ILong = (rand.u64() % 65536).ilong()
      result = result + MPInt.from_ilong(d).digit_shl(k)
      k = k + 1
    end

    result


  fun _run_bench(bench: String, n: USize, n_iter: USize): (U64, U64) =>
    """
    Time `bench` for operands of size `n` digits, repeated `n_iter` times.

    Returns `(elapsed_ns, checksum)` where `elapsed_ns` is the total
    wall-clock nanoseconds for all iterations and `checksum` is a value
    derived from the last result (used to detect silent correctness changes
    across benchmark runs or compiler versions).

    Operands are built from fixed seeds for cross-run comparability.
    Setup work (operand construction, string pre-computation for
    `from_string`) is intentionally excluded from the timed region.
    """
    let a: MPInt = _make_mpint(n, 0xDEAD_BEEF)
    let b: MPInt = _make_mpint(n, 0xCAFE_BABE)

    // Pre-compute the decimal string for from_string (outside the hot loop).
    let a_str: String = if bench == "from_string" then a.string() else "" end

    var result: MPInt = MPInt.from_ilong(0)

    // str_checksum accumulates a checksum for the string benchmark, which
    // does not produce an MPInt result.
    var str_checksum: U64 = 0

    let t0: U64 = Time.nanos()
    var k: USize = 0

    match bench
    | "mul" =>
      while k < n_iter do
        result = a * b
        k = k + 1
      end
    | "mul_schoolbook" =>
      while k < n_iter do
        result = a.mul_schoolbook(b)
        k = k + 1
      end
    | "mul_karatsuba" =>
      while k < n_iter do
        result = a.mul_karatsuba(b)
        k = k + 1
      end
    | "mul_ntt" =>
      while k < n_iter do
        result = a.mul_ntt(b)
        k = k + 1
      end
    | "mul_fft" =>
      while k < n_iter do
        result = a.mul_fft(b)
        k = k + 1
      end
    | "string" =>
      while k < n_iter do
        let s: String = a.string()
        // XOR string length and first byte into checksum to detect regressions.
        str_checksum = str_checksum xor s.size().u64()
        str_checksum = str_checksum xor (try s(0)?.u64() else 0 end)
        k = k + 1
      end
    | "from_string" =>
      while k < n_iter do
        result = try MPInt.from_string(a_str)? else MPInt.from_ilong(0) end
        k = k + 1
      end
    | "pow2" =>
      while k < n_iter do
        result = a * a
        k = k + 1
      end
    | "divrem" =>
      while k < n_iter do
        (result, _) = a.divrem(b)
        k = k + 1
      end
    | "isqrt" =>
      // Square a and add 1 so the input is not a perfect square (exercises
      // the full Newton loop rather than allowing early exit).
      let a2: MPInt = (a * a) + MPInt.from_ilong(1)
      while k < n_iter do
        result = a2.isqrt()
        k = k + 1
      end
    end

    let elapsed: U64 = Time.nanos() - t0

    let checksum: U64 = if bench == "string" then
      str_checksum
    else
      _checksum(result)
    end

    (elapsed, checksum)


  fun _checksum(m: MPInt): U64 =>
    """
    A lightweight XOR checksum over the raw bytes of `m`.

    Used to detect silent correctness changes across benchmark runs
    (e.g., after an optimisation refactor or a compiler upgrade).
    A change in checksum for a fixed operand seed signals that the
    operation no longer produces the same result.
    """
    var acc: U64 = 0
    for byte in m.raw_digits().values() do
      acc = acc xor byte.u64()
    end
    acc


  fun _csv_quote(s: String val): String val =>
    """
    Return `s` wrapped in double quotes if it contains any comma, double-quote,
    carriage-return, or newline character, per RFC 4180.  Internal double-quote
    characters are escaped by doubling them (`"` → `""`).  Returns `s`
    unchanged when no quoting is required.

    This ensures a user comment such as `"added Karatsuba, see PR #42"` is
    stored as a single CSV field rather than splitting across two columns.
    """
    var needs_quote: Bool = false
    for c in s.values() do
      if (c == ',') or (c == '"') or (c == '\n') or (c == '\r') then
        needs_quote = true
        break
      end
    end
    if not needs_quote then
      return s
    end
    let result: String iso = recover
      let r = String(s.size() + 2)
      r.push('"')
      for c in s.values() do
        if c == '"' then
          r.push('"')
        end
        r.push(c)
      end
      r.push('"')
      r
    end
    consume result
