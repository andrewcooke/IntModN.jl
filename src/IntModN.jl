
# a pragmatic library for doing modular arithmetic.

# the aim here is not to encapsulate a large amount of theory, or to describe
# the theoretical relationships between different structures, but to enable
# arithmetic on various types, motivated largely by the practical needs of
# crypto code.

# currently incomplete; the aim is to support rings and fields (prime moduli)
# of integers, rings and factor rings of polynomials over those, and related
# fields (irreducible factor polynomials).

# much thanks to Andreas Noack Jensen
# https://groups.google.com/d/msg/julia-users/Ui977brdzAU/u4rWiDeJv-MJ
# https://gist.github.com/andreasnoackjensen/9132511
# (but since this changes things significanty, mistakes are mine...)


module IntModN

import Base: show, showcompact, zero, one, inv, real, abs, convert,
       promote_rule

export Z, ZField, ZRing, ZR, ZF, GF2


# --- integers modulo n

#     we need two types here because prime moduli have a faster inverse
#     (via euler's theorem) (note - we do NOT test for primality).

abstract Z{N, I<:Integer} <: Number

# you can construct these directly, but code assumes that the integer is
# already reduced mod N, that capacities are ok, etc.

immutable ZRing{N, I<:Integer} <: Z{N,I}
    i::I
end

immutable ZField{N, I<:Integer} <: Z{N,I}
    i::I
end

# but typically you would use one of these contructors

ZR{I<:Integer}(n::Int, i, ::Type{I}) = ZRing{n,I}(validate(n, convert(I, i)))
ZR(n::Int, i::Integer) = ZRing{n, typeof(i)}(validate(n, i))
ZF{I<:Integer}(n::Int, i, ::Type{I}) = ZField{n,I}(validate(n, convert(I, i)))
ZF(n::Int, i::Integer) = ZField{n, typeof(i)}(validate(n, i))
Z2(n) = ZF(2, n)
GF2 = Z2

zero{N,I}(::Type{ZRing{N,I}}) = ZRing{N,I}(zero(I))
one{N,I}(::Type{ZRing{N,I}}) = ZRing{N,I}(one(I))
zero{N,I}(::Type{ZField{N,I}}) = ZField{N,I}(zero(I))
one{N,I}(::Type{ZField{N,I}}) = ZField{N,I}(one(I))

function validate(n, i)
    @assert n > 0 "modulus ($n) too small"
    @assert n <= typemax(typeof(i)) "modulus ($n) too large for $(typeof(i))"
    mod(i, n)
end

convert{X<:Integer}(::Type{X}, z::ZRing) = convert(X, z.i)
convert{X<:Integer}(::Type{X}, z::ZField) = convert(X, z.i)
modulus{N, I}(::ZRing{N, I}) = N
modulus{N, I}(::ZField{N, I}) = N

showcompact{N,I}(io::IO, z::Z{N,I}) = showcompact(op, convert(I, z))
show{N,I}(io::IO, z::Z{N,I}) = print(io, "$(convert(I, z)) mod $(modulus(z))")

=={N,I}(a::Z{N,I}, b::Z{N,I}) = convert(I, a) == convert(I, b)  # equal N
<{N,I}(a::Z{N,I}, b::Z{N,I}) = convert(I, a) < convert(I, b)

real{N,I}(a::Z{N,I}) = real(convert(I, a))
abs{N,I}(a::Z{N,I}) = abs(convert(I, a))

promote_rule{N, I<:Unsigned}(::Type{ZField{N,I}}, ::Type{Uint}) = Uint
promote_rule{N, I<:Integer}(::Type{ZField{N,I}}, ::Type{Int}) = Int
promote_rule{N, I<:Unsigned}(::Type{ZRing{N,I}}, ::Type{Uint}) = Uint
promote_rule{N, I<:Integer}(::Type{ZRing{N,I}}, ::Type{Int}) = Int

lift(c, n, i, f, a) = c(convert(i, mod(f(convert(i, a)), n)))
lift(c, n, i, f, a, b) = c(convert(i, mod(f(convert(i, a), convert(i, b)), n)))
-{N,I}(a::ZRing{N,I}) = lift(ZRing{N,I}, N, I, -, a)
-{N,I}(a::ZField{N,I}) = lift(ZField{N,I}, N, I, -, a)
+{N,I}(a::ZRing{N,I}, b::ZRing{N,I}) = lift(ZRing{N,I}, N, I, +, a, b)
+{N,I}(a::ZField{N,I}, b::ZField{N,I}) = lift(ZField{N,I}, N, I, +, a, b)
-{N,I}(a::ZRing{N,I}, b::ZRing{N,I}) = lift(ZRing{N,I}, N, I, -, a, b)
-{N,I}(a::ZField{N,I}, b::ZField{N,I}) = lift(ZField{N,I}, N, I, -, a, b)
*{N,I}(a::ZRing{N,I}, b::ZRing{N,I}) = lift(ZRing{N,I}, N, I, *, a, b)
*{N,I}(a::ZField{N,I}, b::ZField{N,I}) = lift(ZField{N,I}, N, I, *, a, b)
^{N,I}(a::ZRing{N,I}, p::Integer) = lift(ZRing{N,I}, N, I, x -> powermod(x, p, N), a)
^{N,I}(a::ZField{N,I}, p::Integer) = lift(ZField{N,I}, N, I, x -> powermod(x, p, N), a)

# http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Modular_integers
function inverse{I<:Integer}(a::I, n::I)
    t::I, newt::I = zero(I), one(I)
    r::I, newr::I = n, a
    while newr != 0
        q::I = div(r, newr)
        t, newt = newt, t - q * newt
        r, newr = newr, r - q * newr
    end
    @assert r <= 1 "$a is not invertible mod $n"
    t < 0 ? t + n : t
end

inv{N,I}(a::ZRing{N,I}) = ZRing{N,I}(inverse(N, convert(I, N)))
inv{N,I}(a::ZField{N,I}) = a^(N-2)
/{N,I}(a::Z{N,I}, b::Z{N,I}) = a * inv(b)


function test_constructor()

    @assert string(ZR(3, 2, Int)) == "2 mod 3"
    @assert string(ZR(3, 2)) == "2 mod 3"
    @assert string(ZF(3, 2, Int)) == "2 mod 3"
    @assert string(ZF(3, 2)) == "2 mod 3"

    @assert GF2(1) == Z2(1) == ZField{2,Int}(1) == ZF(2, 1)

    @assert isa(ZR(4, 0x3).i, Uint8)

    try
        bad = ZR(256, 0x0)
        error("expected failure")
    catch e
        @assert isa(e, ErrorException)
        @assert search(string(e), "too large") != 0:-1 
    end
    println("test_constructor ok")
end

function test_arithmetic()

    @assert GF2(1) + GF2(1) == GF2(0)
    @assert zero(ZRing{5,Int}) - one(ZRing{5,Int}) == ZR(5, 4)
    @assert zero(ZField{5,Int}) - one(ZField{5,Int}) == ZF(5, 4)
    @assert inv(ZF(5, 3)) == ZF(5, 2)
    @assert ZF(5, 3) * ZF(5, 2) == one(ZField{5, Int})
    @assert GF2(0)^0 == GF2(1)

    try
        inv(ZR(6, 2))
        error("expected failure")
    catch e
        @assert isa(e, ErrorException)
        @assert search(string(e), "not invertible") != 0:-1 
    end

    println("test_arithmetic ok")
end

function test_matrix_inverse()

    l, o = GF2(1), GF2(0)

    A = [l l;
         o l]
    b = [l, l]
    x = A\b
    @assert x == [o, l]

    # http://math.stackexchange.com/questions/169921/how-to-solve-system-of-li#near-equations-of-xor-operation
    A = [l l l o; 
         l l o l;
         l o l l;
         o l l l]
    b = [l, l, o, l]
    x = A\b
    @assert x == [o, l, o, o]

    println("test_matrix_inverse ok")
end

function test_power()

    # http://en.wikipedia.org/wiki/Diffie%E2%80%93Hellman_key_exchange#Explanation_including_encryption_mathematics
    base = ZF(23, 5)
    @assert (base^6).i == 8
    @assert (base^15).i == 19
    @assert ((base^15)^6).i == 2
    @assert ((base^6)^15).i == 2

    println("test_power ok")
end

function tests()
    test_constructor()
    test_arithmetic()
    test_matrix_inverse()
    test_power()
end

tests()

end
