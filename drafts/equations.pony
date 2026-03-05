primitive QuadraticEquation[F: (FloatingPoint[F] | Complex[F]) & ...]
  """
  Quadratic equation solver trying to limit loss of precision.
  
  The quadratic equations is `a * x^2 + b * x + c = 0`.
  
  Solutions are `x = (-b +/- (b * b - 4 * a * c).sqrt()) / ( 2 * a)`
  """
  
  fun solve(a: F, b: F, c: F): (F, F) =>
    """
    Solves the quadratic equation `a * x^2 + b * x + c = 0` and returns the two
    solutions in a tuple.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    let four = F.from[ISize](4)
    let discrim = (b * b) - (four * a * c)
    let sign =
      match b
      | let _: FloatingPoint[F] => if b > zero then one else -one end
      | let _: Complex[F] => if (b.conjug() * discrim.sqrt()).re() >= zero then one else -one end
    let q = -((b + (sign * discrim.sqrt())) / two)
    
    (q / a, c / q)
    
    
primitive CubicEquation[F: (FloatingPoint[F] | Complex[F]) & ...]
  """
  Cubic equation solver trying to limit loss of precision.
  
  The cubic equation is `x^3 + a * x^2 + b * x + c = 0`.
  """
  
  fun solve(a: F, b:F, c: F): (F, F, F) =>
    """
    Solves the cubic equation is `x^3 + a * x^2 + b * x + c = 0` and returns the
    three solutions in a tuple.
    """
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    let three = F.from[ISize](3)
    let nine =  F.from[ISize](9)
    let twenty_seven = F.from[ISize](27)
    let fifty_four = F.from[ISize](54)
    
    let q = ((a * a) - (three * a * b)) / nine
    let r = (((two * a * a * a) - (nine * a * b)) + (twenty_seven * c)) / fifty_four
    let r2 = r * r
    let q3 = q * q * q
    
    if q.is_real() and r.is_real() and (r2 < q3) then
      let theta = (r / q3.sqrt()).acos()
      let qsqrt = q.sqrt()
      let a3 = a / three
      let pi = F.pi()
      let pi23 = (two * pi) / three
      let theta3 = theta / three
      let x1 = -(two * qsqrt * theta3.cos()) - a3
      let x2 = -(two * qsqrt * (theta3 + pi23).cos()) - a3
      let x2 = -(two * qsqrt * (theta3 - pi23).cos()) - a3
      
      (x1, x2, x3)
    else
      let one_third = one / three
      let r2q3sqrt = (r2 - q3).sqrt()
      let ba = - (r + r2q3sqrt).pow(one_third)
      
      let sign = (r.conjug() * r2q3sqrt).real()
    end
    