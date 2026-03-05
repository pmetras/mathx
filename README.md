# Mathematical library

The library includes various mathematical algorithms.

* `FLimits` class has functions to determines the platform/compiler limits of floating point numbers.
* `Primes` class deals with prime numbers and includes the Miller-Rabin probabilistic primality test.
* `Complex` is a primitive implementing complex numbers arithmetic.
* `FFT` implements the Fast Fourier Transform algorithm.
* `MPInt` and `MPFloat` are multiple precision integer and floating point implementations. A more complete implementation using [GMP](https://gmplib.org/) and [MPFR](https://mpfr.org/) libraries through FFI is also available.

## Setup

You need to install `libgmp-dev` and `libmpfr-dev` if you want to use GMP and MPFR
