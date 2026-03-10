# Floating points platform limits

The class `FLimits` determines the Pony compiler/platform limits for the floating point types (`F32` and `F64`, or eventually `F128` if Pony supports such type in the future). The algorithms used were taken from the MACHAR Fortran source ported to C, and now to Pony. You can read a version of [`machar.c`](machar.c) that compiles with GCC or clang.

The [`platform_limits.pony`](platform_limits.pony) executes the MACHAR algorithms and compare the results with the results of limit functions from `F32` and `F64`.

## Sample output

On linux amd64:

```
Platform limits
===============
Calculated limits for type F32
--------------------------------
On binary hardware (radix == 2), unit = bit

Radix            = 2
Digits           = 24 units
Guard digits     = 0
	0: if floating-point arithmetic rounds, or if it truncates and only `Digit`-base `Radix` digits participate in the post-normalization shift of the floating-point mantissa in multiplication.
	1: if floating-point arithmetic truncates and more than `Digit`-base `Radix` digits participate in the post- normalization shift of the floating-point mantissa in multiplication.
Round style      = 5
	0: if floating-point addition chops
	1: if floating-point addition rounds, but not in the IEEE style
	2: if floating-point addition rounds in the IEEE style
	3: if floating-point addition chops, and there is partial underflow
	4: if floating-point addition rounds, but not in the IEEE style, and there is partial underflow
	5: if floating-point addition rounds in the IEEE style, and there is partial underflow
Machep           = -23
	The exponent for the smallest power of `Radix` (but bounded below by `Digits - 3`) whose sum with 1.0 is greater than 1.0.
Negeps           = -24
	The exponent for the smallest power of `Radix` (but bounded below by `Digits - 3`) whose difference with 1.0 is less than 1.0.
Exponent         = 8 units
Minimal exponent = -126
Maximal exponent = 128
Epsilon          = 1.19209e-07
Negative epsilon = 5.96046e-08
Minimal value    = 1.17549e-38
Maximal value    = 3.40282e+38

Compiler limits for type F32
------------------------------
F.precision2     = 24 bits
F.precision10    = 6
F.min_exp2       = -125
F.min_exp10      = -37
F.max_exp2       = 128
F.max_exp10      = 38
F.epsilon        = 1.19209e-07
F.min_value      = -3.40282e+38
F.max_value      = 3.40282e+38
F.min_normalised = 1.17549e-38

Calculated limits for type F64
--------------------------------
On binary hardware (radix == 2), unit = bit

Radix            = 2
Digits           = 53 units
Guard digits     = 0
	0: if floating-point arithmetic rounds, or if it truncates and only `Digit`-base `Radix` digits participate in the post-normalization shift of the floating-point mantissa in multiplication.
	1: if floating-point arithmetic truncates and more than `Digit`-base `Radix` digits participate in the post- normalization shift of the floating-point mantissa in multiplication.
Round style      = 5
	0: if floating-point addition chops
	1: if floating-point addition rounds, but not in the IEEE style
	2: if floating-point addition rounds in the IEEE style
	3: if floating-point addition chops, and there is partial underflow
	4: if floating-point addition rounds, but not in the IEEE style, and there is partial underflow
	5: if floating-point addition rounds in the IEEE style, and there is partial underflow
Machep           = -52
	The exponent for the smallest power of `Radix` (but bounded below by `Digits - 3`) whose sum with 1.0 is greater than 1.0.
Negeps           = -53
	The exponent for the smallest power of `Radix` (but bounded below by `Digits - 3`) whose difference with 1.0 is less than 1.0.
Exponent         = 11 units
Minimal exponent = -1022
Maximal exponent = 1024
Epsilon          = 2.22045e-16
Negative epsilon = 1.11022e-16
Minimal value    = 2.22507e-308
Maximal value    = 1.79769e+308

Compiler limits for type F64
------------------------------
F.precision2     = 53 bits
F.precision10    = 15
F.min_exp2       = -1021
F.min_exp10      = -307
F.max_exp2       = 1024
F.max_exp10      = 308
F.epsilon        = 2.22045e-16
F.min_value      = -1.79769e+308
F.max_value      = 1.79769e+308
F.min_normalised = 2.22507e-308
```

## C version

To try the C version, compile and run it with:

```sh
# Compile with GCC
$ gcc machar.c -o machar.gcc
# Run
$ ./machar.gcc

# Compile with clang
$ clang machar.c -o machar.clang
# Run
$ ./machar.clang
```
### Sample output

Look at [source code](platform_limits.pony) to get the new names and the explanations of values. On linux amd64:

```
Double  precision MACHAR constants
ibeta  = 2
it     = 53
irnd   = 5
ngrd   = 0
machep = -52
negep  = -53
iexp   = 11
minexp = -1022
maxexp = 1024
eps      2.2204460492503131e-16          0   3CB00000 
epsneg   1.1102230246251565e-16          0   3CA00000 
xmin    2.2250738585072014e-308          0     100000 
xmax    1.7976931348623157e+308   FFFFFFFF   7FEFFFFF 
```
