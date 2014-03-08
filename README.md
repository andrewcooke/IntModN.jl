# IntModN.jl

A pragmatic (meaning incomplete, and written by someone who needed
this before he fully understood it) library for doing modular
arithmetic.

The aim is not to encapsulate a large amount of theory, or to describe
the relationships between different structures, but to enable
arithmetic on various types, motivated largely by the practical needs
of crypto code.

Currently incomplete; the hope is to support rings and fields (prime
moduli) of integers, polynomials over those, and related fields
(irreducible factor polynomials).  At the moment only rings and fields
of integers, and simple polynomials over those, are supported.

## Examples

See the
[tests](https://github.com/andrewcooke/IntModN.jl/blob/master/src/IntModN.jl).

### Simultaneous Equations in GF(2)

Answering [this
question](http://math.stackexchange.com/questions/169921/how-to-solve-system-of-linear-equations-of-xor-operation):

```
@zfield 2 begin
    A = [1 1 1 0; 
         1 1 0 1;
         1 0 1 1;
         0 1 1 1]
    b = [1, 1, 0, 1]
    x = A\b
    @assert x == [0, 1, 0, 0]
end
```

### Polynomial Arithmetic in GF(2)

```
julia> using IntModN

julia> x = X(GF2)
P(ZField{2,Int64},1,0)

julia> a = x^3 + x^2 + 1
P(ZField{2,Int64},1,1,0,1)

julia> b = x^2 + 1
P(ZField{2,Int64},1,0,1)

julia> p, q = divrem(a, b)
(P(ZField{2,Int64},1,1),P(ZField{2,Int64},1,0))

julia> println(p * b + q)
x^3 + x^2 + 1 mod 2
```

## Licence

MIT Licence; (c) 2014 Andrew Cooke andrew@acooke.org

[![Build
Status](https://travis-ci.org/andrewcooke/IntModN.jl.png)](https://travis-ci.org/andrewcooke/IntModN.jl)
Julia 0.3 (trunk).
 
