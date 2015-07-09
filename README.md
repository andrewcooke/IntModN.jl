[![Build
Status](https://travis-ci.org/andrewcooke/IntModN.jl.png)](https://travis-ci.org/andrewcooke/IntModN.jl)
[![Coverage Status](https://coveralls.io/repos/andrewcooke/IntModN.jl/badge.svg)](https://coveralls.io/r/andrewcooke/IntModN.jl)
[![IntModN](http://pkg.julialang.org/badges/IntModN_release.svg)](http://pkg.julialang.org/?pkg=IntModN&ver=release)

# IntModN.jl

A pragmatic (meaning incomplete, and written by someone who needed
this before he fully understood it) library for doing modular
arithmetic.

The aim is not to encapsulate a large amount of theory, or to describe
the relationships between different structures, but to enable
arithmetic on various types, motivated largely by the practical needs
of crypto code.

* [Examples](#examples)
* [Types](#types)
  * [Integers Modulo N](#integers-modulo-n)
  * [Polynomials](#polynomials)
    * [Polynomials with Integral Coefficients](#polynomials-with-integral-coefficients)
    * [Polynomials over GF(2)](#polynomials-over-gf-2)
	

Incomplete; pull requests welcome.

## Examples

See the
[tests](https://github.com/andrewcooke/IntModN.jl/blob/master/src/Tests.jl).

### Simultaneous Equations

Answering [this
question](http://math.stackexchange.com/questions/169921/how-to-solve-system-of-linear-equations-of-xor-operation) (in GF(2)):

```
julia> using IntModN

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

### Polynomial Arithmetic

```
julia> x = X(GF2)
ZP(ZField{2,Int64},1,0)

julia> a = x^3 + x^2 + 1
ZP(ZField{2,Int64},1,1,0,1)

julia> b = x^2 + 1
ZP(ZField{2,Int64},1,0,1)

julia> p, q = divrem(a, b)
(ZP(ZField{2,Int64},1,1),ZP(ZField{2,Int64},1,0))

julia> println(p * b + q)
x^3 + x^2 + 1 mod 2
```

### Rijndael

The multiplication [described
here](http://en.wikipedia.org/wiki/Finite_field_arithmetic#Rijndael.27s_finite_field):

```
julia> x = GF2X(Uint)
GF2Poly{UInt64}(2)

julia> rijndael = x^8 + x^4 + x^3 + x + 1
GF2Poly{UInt64}(283)

julia> print(ZF(rijndael, x^7 + x^6 + x^3 + x) * ZF(rijndael, x^6 + x^4 + x +
1))
1 mod 2 mod x^8 + x^4 + x^3 + x + 1 mod 2
```

### Fast Polynomials in GF(2)

The examples above could have used any modulus.  I chose GF(2) only
because it is common.

However, the following works only in GF(2) (the trade-off for the lack
of flexibility is speed and compactness - these are encoded as bit
patterns):

```
julia> x = GF2X(Uint8)
GF2Poly{Uint8}(2)

julia> a = x^3 + x^2 + 1
GF2Poly{Uint8}(13)

julia> b = x^2 + 1
GF2Poly{Uint8}(5)

julia> p, q = divrem(a, b)
(GF2Poly{Uint8}(3),GF2Poly{Uint8}(2))

julia> @assert a == p * b + q
```

And you can also use these with factor rings:

```
julia> x = GF2X(Uint)
GF2Poly{UInt64}(2)

julia> rijndael = x^8 + x^4 + x^3 + x + 1
GF2Poly{UInt64}(283)

julia> print(FR(x^7 + x^6 + x^3 + x, rijndael) * FR(x^6 + x^4 + x + 1, rijndael))
1 mod x^8 + x^4 + x^3 + x + 1 mod 2
```

However, note that `rinjdael` here requires 9 bits of storage; there is no
representation with an implicit msb.

## Types

`Residue <: Integer` - abstract superclass for (almost) everything below.
Used to provide some common utilities (like automatic promotion from
integers).

### Integers Modulo N

`ZModN{N,I<:Integer} <: Residue` - abstract superclass for integers modulo
some value, where `N` is the modulus, and so typically an `Int` (yes, that's a
integer as a *type*, not a value), and `I`

This has two concrete subclasses, because when `N` is a prime number we can
define a multiplicative inverse.

`ZRing{N, I<:Integer} <: ZModN{N,I}` - the general case.

`ZField{N, I<:Integer} <: ZModN{N,I}` - assumes that `N` is prime, and so
includes division.

These constructors can be used directly, but do not check that arguments are
consistent with assumptions made in the code (values within range, etc).

The associated functions `ZR()` and `ZF()` are more suitable for "normal" use
(but still do not check primality for fields), and include support for factory
functions:

```julia
julia> ZF(3, 5, UInt8)
ZField{3,UInt8}(2)

julia> ZF(3, 5)
ZField{3,Int64}(2)

julia> GF3 = ZF(3)
(anonymous function)

julia> GF3(5)
ZField{3,Int64}(2)
```

The macros `@zring` and `@zfield` can also be used to convert all integers
with scope:

```julia
julia> @zring 4 begin
          A = [1 2 3 4 5]
       end
1x5 Array{IntModN.ZRing{4,Int64},2}:
 1  2  3  0  1
```

### Polynomials

`Poly <: Residue` - abstract superclass for polynomials.  All share some basic
conventions about accessing coefficients (with `[]`) and iterators.

The types below all form rings, not fields, because polynomials do not have
inverses.

Note: Originally, the code used
[Polynomial.jl](https://github.com/vtjnash/Polynomial.jl), but that had some
weird design decisions so I wrote my own code.  Since then,
[Polynomials.jl](https://github.com/Keno/Polynomials.jl) fixed some of the
issues, so at some point it may make sense to revert to that package.

#### Polynomials With Integral Coefficients

`ZPoly{I<:Integer} <: Poly` - a simple wrapper around an array of integral
coefficients (including `ZModN` subclasses).  The coefficients are in the
"usual" order, so `[i]` gives the ith coefficient, and the leading coefficient
is always non-zero (or the array is empty).

As with integers mod N, the constructor can be used directly, but it is
generally preferable to use `ZP()`, which has various forms.

In addition, there's support for the natural syntax `x^n...` via `X()`:

```julia
julia> x = X(ZF(2))
ZP(IntModN.ZField{2,Int64},1,0)

julia> x^3 + x
ZP(IntModN.ZField{2,Int64},1,0,1,0)
```

#### Polynomials over GF(2)

`GF2Poly{U<:Unsigned} <: Poly` - specialized support for polynomials over
GF(2).  Coefficients can only be 0 or 1, so we can use bit fields (integers)
for their values.

As always, you can use the constructor directly, or the utilities `GF2P()` and
`GF2X()`.

The bit pattern can be displayed with `bits()` and addition is binary xor:

```julia
julia> x = GF2X(Uint8)
GF2Poly{UInt8}(2)

julia> a = x^7 + x^3
GF2Poly{UInt8}(136)

julia> b = x^3 + x^2 + 1
GF2Poly{UInt8}(13)

julia> a+b
GF2Poly{UInt8}(133)

julia> bits(a), bits(b), bits(a+b)
("10001000","00001101","10000101")
```

### Quotient (Factor) Rings

These used to be a spearate type, but can now be handled as `ZRing()` and
`ZField)(` with polynomial arguments.  The latter is appropriate when the
ideal is irreducible (maximal) (I think).

See the [Rijndael](#rijndael) example.



