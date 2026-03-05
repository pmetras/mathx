// Test GMP/MPFR classes

use "../../mathx/gmp"

actor Main
  new create(env: Env) =>
    printsize(env)
    let a = MPFloat(2.0)
    let pa = a.get_precision()
    env.out.print("a = " + a.string())
    env.out.print("Precision = " + pa.string())
    let b = MPFloat.from_string("3.00000000000001", 50)
    let pb = b.get_precision()
    env.out.print("b = " + b.string())
    env.out.print("Precision = " + pb.string())
    let c = a + b
    env.out.print("c = " + c.string())
    let pc = c.get_precision()
    env.out.print("Precision = " + pc.string())

  fun printsize(env: Env) =>
    env.out.print("USize = " + USize(0).bitwidth().string())
    env.out.print("ULong = " + ULong(0).bitwidth().string())
    env.out.print("ISize = " + ISize(0).bitwidth().string())
    env.out.print("ILong = " + ILong(0).bitwidth().string())
