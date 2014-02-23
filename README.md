# IntModN.jl

A type, and associated arithmetic operstors, for groups of integers modulo
some value, N.  Only simple integers are supported (no vector spaces or
polynomials).

If N is prime then the group of possible values for the type form a field and
the multiplicative inverse is defined for all non-zero values.  If N is not
prime then an exception may be raised (unless the particular value being
inverted is co-prime with the modulus).

## Types

The base type is `Z{N,I<:Integer}` where N is an `Int` that defines the
modulus and `I` is the type used to store the values.  So `Z{5,Int}(3)`
creates the value `3 mod 5`, stored internally as an `Int`.

Since most uses will want to use `Int` for storage (since it is both fast and
easy to use with literal values), there is a type alias `GF{N}`.  So the above
value is equivalent to `GF{5}(3)`.

Finally, the type alias `GF2` defines the common case.

## Constructors

* Explicitly give all type parameters: `Z{5, Int}(3)`

* Storage type inferred from value: `Z(0x3, 5) == Z{5, UInt8}(0x3)`

* Type alias: `GF{5}(3) == Z{5, Int}(3)`

## Example

To solve [the problem
here](http://math.stackexchange.com/questions/169921/how-to-solve-system-of-linear-equations-of-xor-operation):

```
    l, o = one(GF2), zero(GF2)
    A = [l l l o; 
         l l o l;
         l o l l;
         o l l l]
    b = [l, l, o, l]
    x = A\b
    @assert x == [o, l, o, o]
```

[![Build Status](https://travis-ci.org/andrewcooke/IntModN.jl.png)](https://travis-ci.org/andrewcooke/IntModN.jl)
