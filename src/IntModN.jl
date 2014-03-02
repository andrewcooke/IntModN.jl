
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
# Pkg.clone("https://github.com/astrieanna/TypeCheck.jl.git")
#using TypeCheck

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

# but typically you would use one of these constructors

ZR{I<:Integer}(n::Int, i, ::Type{I}) = ZRing{n,I}(validate(n, convert(I, i)))
ZR(n::Int, i::Integer) = ZRing{n, typeof(i)}(validate(n, i))
ZF{I<:Integer}(n::Int, i, ::Type{I}) = ZField{n,I}(validate(n, convert(I, i)))
ZF(n::Int, i::Integer) = ZField{n, typeof(i)}(validate(n, i))
Z2(n) = ZF(2, n)
GF2 = Z2

# and all the methods that make this a number

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

showcompact{N,I}(io::IO, z::Z{N,I}) = showcompact(io, convert(I, z))
show{N,I}(io::IO, z::Z{N,I}) = print(io, "$(convert(I, z)) mod $(modulus(z))")

=={N,I}(a::Z{N,I}, b::Z{N,I}) = convert(I, a) == convert(I, b)  # equal N
<={N,I}(a::Z{N,I}, b::Z{N,I}) = convert(I, a) <= convert(I, b)
<{N,I}(a::Z{N,I}, b::Z{N,I}) = convert(I, a) < convert(I, b)

real{N,I}(a::Z{N,I}) = a
abs{N,I}(a::Z{N,I}) = a

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


function test_z_constructor()

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
    println("test_z_constructor ok")
end

function test_z_arithmetic()

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

    println("test_z_arithmetic ok")
end

function test_z_matrix_inverse()

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

    println("test_z_matrix_inverse ok")
end

function test_z_power()

    # http://en.wikipedia.org/wiki/Diffie%E2%80%93Hellman_key_exchange#Explanation_including_encryption_mathematics
    base = ZF(23, 5)
    @assert (base^6).i == 8
    @assert (base^15).i == 19
    @assert ((base^15)^6).i == 2
    @assert ((base^6)^15).i == 2

    println("test_z_power ok")
end

function test_z_coverage()
    println("Integer", methodswithdescendants(Integer, lim=20))
    println("ZRing", methodswithdescendants(ZRing, lim=20))
    println("test_z_coverage ok")
end

function tests_z()
    test_z_constructor()
    test_z_arithmetic()
    test_z_matrix_inverse()
    test_z_power()
#    test_z_coverage()
end


# --- polynomials over integers modulo n

abstract P{N, T<:Z}

immutable PRing{N, T<:Z} <:P{N, T}    
    a::Array{T,1}
    function PRing(a)
        @assert length(a) == N "wrong length ($N != length($a))"
        @assert N > 0 "empty polynomial"
        new(a)
    end
end

# constructors

PR{T<:Z}(a::T...) = PRing{length(a), T}([a...])
function PGF2(n::Int, indices::Int...)  # indices are 1-indexed
    a, o = zeros(ZField{2,Int}, n), GF2(1)
    for i in indices
        a[i] = o
    end
    PR(a...)
end

# number methods

zero{N,T<:Z}(::Type{PRing{N,T}}) = PRing{N,T}(zeros(T, N))
one{N,T<:Z}(::Type{PRing{N,T}}) = PRing{N,T}(ones(T, N))

# does not create a new copy
convert{N,T<:Z}(::Type{Array{T,1}}, a::PRing{N,T}) = a.a

showcompact{N,T<:Z}(io::IO, p::PRing{N,T}) = showcompact(io, convert(Array{T,1}, p))
function show{N,T<:Z}(io::IO, p::PRing{N,T})
    for (i, x) in enumerate(p.a)
        if i > 1
            print(io, " + ")
        end
        showcompact(io, x)
        if i == 2
            print(io, " x")
        elseif i > 2
            print(io, " x^$(i-1)")
        end
    end
    print(io, " mod $(modulus(p.a[1]))")  # array length > 0
end

=={N,T<:Z}(a::P{N,T}, b::P{N,T}) = a.a == b.a
<={N,T<:Z}(a::P{N,T}, b::P{N,T}) = a.a <= b.a
<{N,T<:Z}(a::P{N,T}, b::P{N,T}) = a.a < b.a


function test_p_constructor()
    @assert PRing{2,ZField{2,Int}}([GF2(0), GF2(1)]) == PR(GF2(0), GF2(1)) == PGF2(2, 2)
    @assert PR(GF2(0), GF2(0)) == zero(PRing{2, ZField{2,Int}})
    println("test_p_constructor ok")
end

function test_p_type()
    @assert convert(Array{ZField{2,Int},1}, PGF2(2)) == [GF2(0), GF2(0)]
    @assert string(PGF2(3,1)) == "1 + 0 x + 0 x^2 mod 2" 
    println("test_p_type ok")
end

function tests_p()
    test_p_constructor()
    test_p_type()
end


function tests()
    tests_z()
    tests_p()
end

# run by travis (see .travis.yml in root dir)
#tests()

end
