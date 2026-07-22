"""
Extended mathematical types for Pony.

* `Complex[F: (Float & FloatingPoint[F]) = F64]`: complex numbers
* `FFT[F: (Float & FloatingPoint[F]) = F64]`: Fourier transforms
* `FLimits[F: FloatingPoint[F] val]`: numerical limits of the platform
* `Modular[A: UnsignedInteger[A] val = USize]`: modular arithmetic
* `NTT[A: (UnsignedInteger[A] & Any val) = USize]`: number-theoretic transforms
* `Prime[A: UnsignedInteger[A] val = USize]`: prime operations
* `UnsignedComp[U: UnsignedInteger[U] val = USize]`: more operations on unsigned integers

Arbitrary-precision numbers are in the `bignum` package (a sibling of `mathx`),
which depends on `mathx` for `FFT` and `NTT`.

TODO:
* `Matrix` and `Vector`
* `Polynom`
* Differentiation and integration
"""