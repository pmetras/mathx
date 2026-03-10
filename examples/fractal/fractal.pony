"""
Draw an ASCII fractal of the Mandelbrot set.

See http://www.hiddendimension.com/FractalMath/Divergent_Fractals_Main.html and
many other resources on the Web for explanations, or read the comments below.
The [https://en.wikipedia.org/wiki/Mandelbrot_set](Wikipedia page) contains
a picture of the [first published picture](https://en.wikipedia.org/wiki/Mandelbrot_set#/media/File:Mandel.png)
of the Mandelbrot, in 1978 set that we try toe reproduce with this example.

Note that we are using `Complex[F64]` because that's a test example and it would
have been enough to use `Complex[F32]` instead for ASCII art... And faster.
You can adapt this example to send the result to [GNUplot](https://www.gnuplot.info/)
for more colourful results.
"""

use "../../mathx"
use "../../assertx"

use "debug"
use "collections"

actor Main
  """
  Main entry point for drawing an ASCII Mandelbrot fractal.
  """
  let _out: OutStream
    """
    The output stream to print the fractal plot.
    """


  new create(env: Env) =>
    """
    Plot the fractal on the standard output as ASCII art.
    """
    _out = env.out
    plot(this~fractal(), 132, -2.5, 1.5, 60, -1.5, 1.5)
    
  
  fun plot(f: {(Complex[F64]): USize},
           x_steps: USize = 132, x_min: F64 = -2.5, x_max: F64 = 1.5,
           y_steps: USize = 60, y_min: F64 = -1.5, y_max: F64 = 1.5) =>
    """
    Plots fractal function `f` on the plan `[x_min, x_max] x [y_min, y_max]`.
    """
    try
      Assert(x_steps != 0, "The number of steps on x axis can't be null.")?
      Assert(y_steps != 0, "The number of steps on y axis can't be null.")?
    end
    
    let x_inc = (x_max - x_min) / F64.from[USize](x_steps)
    let y_inc = (y_max - y_min) / F64.from[USize](y_steps)
    
    for y in Range[F64](y_min, y_max, y_inc) do
      for x in Range[F64](x_min, x_max, x_inc) do
        let c = Complex[F64](x, y)
        let colour = f(c)
        
        match colour
        | if colour < 25 => _out.write(" ")
        | if colour < 50 => _out.write(".")
        | if colour < 75 => _out.write("o")
        | if colour < 100 => _out.write("*")
        | if colour < 150 => _out.write("O")
        | if colour < 200 => _out.write("@")
        else
          _out.write("#")
        end
      end
      _out.print("")
    end
    _out.print("Mandelbrot set on complex plan [" + x_min.string() + "; " +
               x_max.string() + "] x [" + y_min.string() + "; " + y_max.string() + "]")
    
    
  fun fractal(p: Complex[F64]): USize =>
    """
    The equation of the fractal. Returns the speed of divergence of the function
    on a scale `[0, 255]` that is used to determine the colour to plot the point
    `p`.
    
    We use the Mandelbrot equation `z_n+1 = z^2_n + p` and see how fast the
    the modulus of z goes larger than 2.0, meaning that the recurrence will
    diverge. The function returns the value of `n`.
    """
    var n: USize = 0
    var z = Complex[F64](0.0, 0.0)
    var abs = z.abs()
    while (n < 250) and (abs < 2.0) do
      z = (z * z) + p
      n = n + 1
      abs = z.abs()
    end
    n
