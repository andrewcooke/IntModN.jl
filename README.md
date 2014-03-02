# IntModN.jl

A pragmatic (meaning incomplete, and written by someone who needed this before
he fully understood it) library for doing modular arithmetic.

The aim is not to encapsulate a large amount of theory, or to describe the
theoretical relationships between different structures, but to enable
arithmetic on various types, motivated largely by the practical needs of
crypto code.

Currently incomplete; the hope is to support rings and fields (prime moduli)
of integers, rings and factor rings of polynomials over those, and related
fields (irreducible factor polynomials).  At the moment only rings and fields
of integers are supported.

## Examples

See the
[tests](https://github.com/andrewcooke/IntModN.jl/blob/master/src/IntModN.jl).

### Matrix inversion in GF(2)

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

## Licence

MIT Licence; (c) 2014 Andrew Cooke andrew@acooke.org

[![Build Status](https://travis-ci.org/andrewcooke/IntModN.jl.png)](https://travis-ci.org/andrewcooke/IntModN.jl)
 
