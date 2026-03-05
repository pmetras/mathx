# GEMINI.md - mathx Project Context

This project is a comprehensive mathematical library for the Pony programming language, providing features ranging from arbitrary-precision arithmetic to complex numbers and Fast Fourier Transforms (FFT).

## Project Overview

- **Purpose:** A high-performance mathematical library for Pony.
- **Key Components:**
  - `MPInt`: Arbitrary-precision integers (implemented using base-65536 `U16` digits for FFT-optimized multiplication).
  - `MPFloat`: Arbitrary-precision floating-point numbers.
  - `FFT`: Fast Fourier Transform implementation.
  - `Primes`: Prime number utilities (Miller-Rabin, factorization).
  - `Complex`: Complex number arithmetic.
  - `Modular`: Modular arithmetic and GCD/LCM.
  - `GMP/MPFR`: Wrappers for GNU Multiple Precision and Multiple Precision Floating-Point Reliable libraries (via FFI).
  - `Limits`: Platform and compiler floating-point limits determination.

## Building and Running

The project uses a `Makefile` for its build system. By default, it builds in `debug` mode.

### Key Commands

- `make unit-tests`: Compiles and runs the unit test suite.
- `make build-examples`: Compiles all example programs in the `examples/` directory.
- `make test`: Runs both unit tests and builds examples.
- `make docs`: Generates API documentation in `build/mathx-docs`.
- `make clean`: Removes the current build directory.
- `make realclean`: Removes all build artifacts and documentation.

### Configuration

- `config=release`: Build in release mode (e.g., `make config=release test`).
- `config=debug`: Build in debug mode (default).

### Binaries Location

- Binaries are placed in `build/debug/` or `build/release/`.
- Test binary: `build/[config]/tests`.
- Example binaries: `build/[config]/[example_name]`.

## Dependencies

- **Pony Compiler (`ponyc`):** Required to build the project.
- **System Libraries:**
  - `libgmp-dev`: Required for GMP-based multi-precision features.
  - `libmpfr-dev`: Required for MPFR-based multi-precision features.
- **Pony Packages:** The project expects sibling directories `../assert` and `../pony_test` for building tests.

## Development Conventions

- **Testing:** Unit tests are located in `tests/` and use the `pony_test` package. Each major component has a corresponding `_tests_[component].pony` file. New features or optimizations must be verified with corresponding test cases.
- **Code Style:** Follows standard Pony idioms. Extensive documentation is provided in docstrings within the `.pony` files.
- **Multi-Precision Strategy:** The pure Pony implementation of `MPInt` uses FFT for multiplication (convolutions), which limits `MPInt` size to approximately 10^6 base-65536 digits due to `F64` precision limits. For higher precision or performance, use the GMP/MPFR wrappers.
- **Examples:** Demonstrations of library features (e.g., fractals, prime counting) are located in `examples/`.
- **Style**: Follows the extended [Pony standard library Style Guide](./STYLE_GUIDE.md).

## Behavior rules

Apply rules from [AGENTS.md](./AGENTS.md).
