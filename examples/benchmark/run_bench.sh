#!/usr/bin/env bash
# run_bench.sh — Run the MPInt benchmark suite and append results to bench.csv.
#
# Results are written to examples/benchmark/bench.csv (next to this script).
# The file is created with a CSV header the first time it does not exist.
#
# Usage:
#   ./run_bench.sh [--debug] [--comment=TEXT | "comment text"] [bench flags]
#
# Examples:
#   # First baseline — no optimisations applied yet
#   ./run_bench.sh "no optimisations yet"
#
#   # After implementing Karatsuba dispatch:
#   ./run_bench.sh "automatic mul dispatch via karatsuba" --bench=mul --iter=10
#
#   # Run only the multiplication benchmarks at specific sizes, release build:
#   ./run_bench.sh --comment="NTT test" --bench=mul,mul_fft --size=100,500,2000
#
#   # Debug build (asserts enabled, no compiler optimisations):
#   ./run_bench.sh --debug "debug baseline"
#
#   # Print the CSV header and exit (useful to re-initialise the log file):
#   ./run_bench.sh --header-only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CSV="$SCRIPT_DIR/bench.csv"
COMMENT=""
CONFIG="release"
HEADER_ONLY=false
BENCH_ARGS=()

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------

usage() {
  cat <<EOF
Usage: $(basename "$0") [--debug] [--comment=TEXT | TEXT] [bench flags]

Options:
  --debug            Use the debug build (default: release).
  --comment=TEXT     Free-text note recorded in the 'comment' CSV column.
  TEXT               First non-flag positional argument is used as the comment.
  --header-only      Write the CSV header to bench.csv and exit.
  -h, --help         Show this help and exit.

Bench flags (passed through to the benchmark binary):
  --bench=n1,n2,...  Comma-separated benchmark names (default: all).
                     Names: mul  mul_karatsuba  mul_fft  string  from_string
                            pow2  divrem  isqrt
  --size=n1,n2,...   Operand sizes in base-65536 digits (default: 10,100,500,2000).
  --iter=N           Timed iterations per (bench, size) pair (default: 5).
  --max-mem=N        Memory limit in GB; skip pairs exceeding it (default: 8).

Output:
  Results are appended to: $CSV
  Columns: ts_s, build, bench, n_digits, n_iter, elapsed_ns, ns_per_iter,
           checksum, comment
EOF
  exit 0
}

for arg in "$@"; do
  case "$arg" in
    -h|--help)
      usage
      ;;
    --debug)
      CONFIG="debug"
      ;;
    --header-only)
      HEADER_ONLY=true
      ;;
    --comment=*)
      COMMENT="${arg#--comment=}"
      ;;
    --bench=*|--size=*|--iter=*|--max-mem=*)
      BENCH_ARGS+=("$arg")
      ;;
    -*)
      echo "$(basename "$0"): unknown flag: $arg" >&2
      echo "Run with --help for usage." >&2
      exit 1
      ;;
    *)
      # First bare word becomes the comment.
      if [ -z "$COMMENT" ]; then
        COMMENT="$arg"
      else
        echo "$(basename "$0"): unexpected argument: $arg" >&2
        echo "Run with --help for usage." >&2
        exit 1
      fi
      ;;
  esac
done

# ---------------------------------------------------------------------------
# Locate benchmark binary
# ---------------------------------------------------------------------------

BINARY="$SCRIPT_DIR/../../build/$CONFIG/benchmark"

if [ ! -f "$BINARY" ]; then
  echo "$(basename "$0"): binary not found: $BINARY" >&2
  echo "" >&2
  if [ "$CONFIG" = "release" ]; then
    echo "Build with:  make config=release build-examples" >&2
  else
    echo "Build with:  make build-examples" >&2
  fi
  exit 1
fi

# ---------------------------------------------------------------------------
# Handle --header-only: (re-)write the header and exit
# ---------------------------------------------------------------------------

if $HEADER_ONLY; then
  "$BINARY" --header-only > "$CSV"
  echo "Header written to $CSV"
  exit 0
fi

# ---------------------------------------------------------------------------
# Create bench.csv with header if it does not exist yet
# ---------------------------------------------------------------------------

if [ ! -f "$CSV" ]; then
  "$BINARY" --header-only > "$CSV"
  echo "Created $CSV"
fi

# ---------------------------------------------------------------------------
# Run the benchmarks and append results
# ---------------------------------------------------------------------------

COMMENT_ARG=()
if [ -n "$COMMENT" ]; then
  COMMENT_ARG=("--comment=$COMMENT")
fi

"$BINARY" \
  --no-header \
  "${COMMENT_ARG[@]}" \
  "${BENCH_ARGS[@]}" \
  >> "$CSV"

ROWS=$(( $(wc -l < "$CSV") - 1 ))
echo "Appended to $CSV ($ROWS data row(s) total)"
if [ -n "$COMMENT" ]; then
  echo "Comment: $COMMENT"
fi
