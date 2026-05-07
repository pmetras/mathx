"""
Extended mathematical types for Pony.

* `Complex[F: (Float & FloatingPoint[F]) = F64]`: complex numbers
* `FFT[F: (Float & FloatingPoint[F]) = F64]`: Fourier transforms
* `FLimits[F: FloatingPoint[F] val]`: numerical limits of the platform
* `Modular[A: UnsignedInteger[A] val = USize]`: modular arithmetic
* `MPFContext`: precision and rounding policy for arbitrary-precision floating-point arithmetic.
* `MPFloat`: real numbers with arbitrary precision
* `MPFMathLib`: mathematical functions on `MPFloat`
* `MPFMathLibRep`: mathematical fonctions on `MPFRep`
* `MPFRep`: representation of an arbitrary-precision floating-point number
* `MPInt`: multiple-precision integer
* `NTT[A: (UnsignedInteger[A] & Any val) = USize]`: number theoritic transforms
* `Prime[A: UnsignedInteger[A] val = USize]`: prime operations
* `RoundingMode`: rounding modes for `MPFloat`
* `UnsignedComp[U: UnsignedInteger[U] val = USize]`: more operations on unsigned integers

TODO:
* `Matrix` and `Vector`
* `Polynom`
* Differentiation and integration
"""