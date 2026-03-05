# ASCII art Mandelbrot fractal

This example draws the well-known Mandelbrot fractal with ASCII characters. The interesting part is not much in the ASCII ploting but the way the fractal points are calculated in the complex plan.

The Mandelbrot set is defined by the points `c` of the complex plan for which the recurrence suite does not diverge:

  * `z_0 = 0`
  * `z_n+1 = z^2_n + c`

The suite of `z_n` diverges if `z.abs()` increases to infinite. It has been shown that a minimal condition for this is that there exist a `n` so that `z_n.abs() > 2`.

We calculate this `n` value for the points in an interval of the complex plan, and print the speed of divergence wirh various ASCII characters.

