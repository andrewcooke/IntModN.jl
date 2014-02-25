# IntModN.jl

A type, and associated arithmetic operators, for rings of integers modulo some
value, N.  Only simple integers are supported (no vector spaces or
polynomials).

If N is prime then the group of possible values for the type form a field and
the multiplicative inverse is defined for all non-zero values.  If N is not
prime then an exception may be raised (unless the particular value being
inverted is co-prime with the modulus).

## Type and Aliases

The base type is `Z{N,I<:Integer}` where N is an `Int` that defines the
modulus and `I` is the type used to store the values.  So `Z{5,Int}(3)`
creates the value `3 mod 5`, stored internally as an `Int`.

Most uses will want to use `Int` for storage (since it is both fast and easy
to use with literal values), so there is a type alias `GF{N}`.  The value
above can be written as `GF{5}(3)`.

Finally, the type alias `GF2` defines the common case of values modulo 2
(where addition and subtraction are XOR and multiplication is AND).

## Constructors

* Explicitly give all type parameters: `Z{5, Int}(3)`;

* Helper that infers from value: `Z(5, 0x3) == Z{5, UInt8}(0x3)`;

* Type alias: `GF{5}(3) == Z{5, Int}(3)`.

The first form requires that the argument be within the range 0-N.  Other
forms normalize as necessary.

## Example - Boolean Algebra

To solve [the problem
here](http://math.stackexchange.com/questions/169921/how-to-solve-system-of-linear-equations-of-xor-operation):

```
    1 = x $ y $ z
    1 = x $ y $ w
    0 = x $ w $ z
    1 = w $ y $ z
```

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

```
    x, y, z, w = 0, 1, 0, 0
```
## Example - DH Key Exchange

From [the example
here](http://en.wikipedia.org/wiki/Diffie%E2%80%93Hellman_key_exchange#Explanation_including_encryption_mathematics):

```
    base = Z(23, 5)
    @assert (base^6).n == 8
    @assert (base^15).n == 19
    @assert ((base^15)^6).n == 2 == ((base^6)^15).n
```

## Licence

MIT Licence; (c) 2014 Andrew Cooke andrew@acooke.org

[![Build Status](https://travis-ci.org/andrewcooke/IntModN.jl.png)](https://travis-ci.org/andrewcooke/IntModN.jl)
 
