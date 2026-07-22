"""
Arbitrary-precision numbers for Pony.

This package provides:

- `MPInt` — arbitrary-precision integers (base-65536 pure-Pony implementation)
- `MPFloat` — arbitrary-precision floating-point (base-256 pure-Pony implementation)
- `BigReal[T]` — interface satisfied by any arbitrary-precision float backend
- `RoundingMode` — IEEE 754 / MPFR rounding modes (defined in `mathx`, re-exported here)

The pure-Pony implementations (`MPInt`, `MPFloat`) are the default and require
no external libraries. An FFI-backed substitute using the GNU MPFR library is
available in the `bignum/gmp` sub-package.

This package depends on `mathx` for `FFT` and `NTT`, which are used internally
by `MPInt` for fast multiplication.
"""
