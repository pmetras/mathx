// Various integer functions

use "collections"


primitive MPIntMath

    fun _recfact(start: ILong, n: ILong): MPInt =>
        """
        Helper function to calculate factorial recursively.
        """
        if n <= 16 then 
            var r = MPInt.from[ILong](start)
            for i in Range[ILong](start + 1, start + n) do
                r = r * MPInt.from[ILong](i)
            end
            return r
        end

        let i = n / 2
        _recfact(start, i) * _recfact(start + i, n - i)


    fun factorial(n: ILong): MPInt =>
        """
        Calculate factorial of n, i.e. `n!` by split recursion.

        See http://www.luschny.de/math/factorial/conclusions.html
        for a presentation of factorial algorithms, and if you want to
        implement the parallel swing factorial algorithm.
        """
        _recfact(1, n)