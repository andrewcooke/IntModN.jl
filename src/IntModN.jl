
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

export ZModN, ZField, ZRing, ZR, ZF, GF2, @zring, @zfield


# --- integers modulo n

#     we need two types here because prime moduli have a faster inverse
#     (via euler's theorem) (note - we do NOT test for primality).

abstract ZModN{N, I<:Integer} <: Number

# you can construct these directly, but code assumes that the integer is
# already reduced mod N, that capacities are ok, etc.

immutable ZRing{N, I<:Integer} <: ZModN{N,I}
    i::I
end

immutable ZField{N, I<:Integer} <: ZModN{N,I}
    i::I
end

# but typically you would use one of these constructors

ZR{I<:Integer}(n::Int, i, ::Type{I}) = ZRing{n,I}(validate(n, convert(I, i)))
ZR(n, i::Integer) = ZRing{n, typeof(i)}(validate(n, i))
ZR(n) = i::Integer -> ZR(n, i)  # used to construct constructors
ZF{I<:Integer}(n::Int, i, ::Type{I}) = ZField{n,I}(validate(n, convert(I, i)))
ZF(n, i::Integer) = ZField{n, typeof(i)}(validate(n, i))
ZF(n) = i::Integer -> ZF(n, i)  # used to construct constructors

GF2 = ZF(2)

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
modulus{N, I<:Integer}(::Type{ZModN{N, I}}) = N
modulus{T<:ZModN}(::Type{T}) = modulus(super(T))  # jameson type chain
modulus{T<:ZModN}(::T) = modulus(T)

showcompact{N,I}(io::IO, z::ZModN{N,I}) = showcompact(io, convert(I, z))
show{N,I}(io::IO, z::ZModN{N,I}) = print(io, "$(convert(I, z)) mod $(modulus(z))")

=={N,I}(a::ZModN{N,I}, b::ZModN{N,I}) = convert(I, a) == convert(I, b)  # equal N
<={N,I}(a::ZModN{N,I}, b::ZModN{N,I}) = convert(I, a) <= convert(I, b)
<{N,I}(a::ZModN{N,I}, b::ZModN{N,I}) = convert(I, a) < convert(I, b)

real{N,I}(a::ZModN{N,I}) = a
abs{N,I}(a::ZModN{N,I}) = a

promote_rule{N, I<:Unsigned}(::Type{ZField{N,I}}, ::Type{Uint}) = Uint
promote_rule{N, I<:Integer}(::Type{ZField{N,I}}, ::Type{Int}) = Int
promote_rule{N, I<:Unsigned}(::Type{ZRing{N,I}}, ::Type{Uint}) = Uint
promote_rule{N, I<:Integer}(::Type{ZRing{N,I}}, ::Type{Int}) = Int

liftz(c, n, i, f, a) = c(convert(i, mod(f(convert(i, a)), n)))
liftz(c, n, i, f, a, b) = c(convert(i, mod(f(convert(i, a), convert(i, b)), n)))
-{N,I}(a::ZRing{N,I}) = liftz(ZRing{N,I}, N, I, -, a)
-{N,I}(a::ZField{N,I}) = liftz(ZField{N,I}, N, I, -, a)
+{N,I}(a::ZRing{N,I}, b::ZRing{N,I}) = liftz(ZRing{N,I}, N, I, +, a, b)
+{N,I}(a::ZField{N,I}, b::ZField{N,I}) = liftz(ZField{N,I}, N, I, +, a, b)
-{N,I}(a::ZRing{N,I}, b::ZRing{N,I}) = liftz(ZRing{N,I}, N, I, -, a, b)
-{N,I}(a::ZField{N,I}, b::ZField{N,I}) = liftz(ZField{N,I}, N, I, -, a, b)
*{N,I}(a::ZRing{N,I}, b::ZRing{N,I}) = liftz(ZRing{N,I}, N, I, *, a, b)
*{N,I}(a::ZField{N,I}, b::ZField{N,I}) = liftz(ZField{N,I}, N, I, *, a, b)
^{N,I}(a::ZRing{N,I}, p::Integer) = liftz(ZRing{N,I}, N, I, x -> powermod(x, p, N), a)
^{N,I}(a::ZField{N,I}, p::Integer) = liftz(ZField{N,I}, N, I, x -> powermod(x, p, N), a)

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
/{N,I}(a::ZModN{N,I}, b::ZModN{N,I}) = a * inv(b)

# rewrite literals mod n 
# (these rewrite ALL ints, including args to literal ZF(...), etc)

macro zring(n::Int, code)
    :($(esc(rewrite_ints(n, x -> ZR(n, x), code))))
end

macro zfield(n::Int, code)
    :($(esc(rewrite_ints(n, x -> ZF(n, x), code))))
end

function rewrite_ints(n, c, val)
    if isa(val, Integer)
        :($(c(val)))
    else
        if isa(val, Expr) && val.head != :line
            val.args = [rewrite_ints(n, c, a) for a in val.args]
        end
        val
    end
end

# tests

function test_z_constructor()

    @assert string(ZR(3, 2, Int)) == "2 mod 3"
    @assert string(ZR(3, 2)) == "2 mod 3"
    @assert string(ZF(3, 2, Int)) == "2 mod 3"
    @assert string(ZF(3, 2)) == "2 mod 3"
    
    Z2 = ZF(2)
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

    @zfield 2 begin
        A = [1 1;
             0 1]
        b = [1, 1]
        x = A\b
        @assert x == [0, 1]
    end

    # http://math.stackexchange.com/questions/169921/how-to-solve-system-of-linear-equations-of-xor-operation
    @zfield 2 begin
        A = [1 1 1 0; 
             1 1 0 1;
             1 0 1 1;
             0 1 1 1]
        b = [1, 1, 0, 1]
        x = A\b
        @assert x == [0, 1, 0, 0]
    end

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

function test_z_macros()
    begin
        b = 7
    end
    @assert b == 7 
    @zfield 5  begin
        b = 1 + 2 * 3
        a = b / 3
    end
    @assert a == ZF(5, 4)
    @assert ZR(3, 2) == @zring 3 1 + 4
    println("test_z_macros ok")
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
    test_z_macros()
#    test_z_coverage()
end


# --- polynomials over integers modulo n

immutable ZPoly{T<:ZModN}
    a::Array{T,1}
    ZPoly(a) = length(a) == 0 || a[end] > zero(T) ? new(a) : error("zero leading coeff")
end

# constructors

ZP() = error("provide at least one (zero?) coeff")
ZP(c::Function) = error("provide at least one (zero?) coeff")
function ZP{T<:ZModN}(coeffs::(Int, T)...)
    coeffs = collect(filter(c -> c[2] > zero(T), coeffs))
    m = length(coeffs) == 0 ? 0 : maximum([i for (i, z) in coeffs])
    a = zeros(T, m)
    for (i, z) in coeffs
        a[i] = z
    end
    ZPoly{T}(a)
end
ZP{I<:Integer}(c::Function, coeffs::(Int, I)...) = ZP([(i, c(n)) for (i, n) in coeffs]...)
ZP{T<:ZModN}(z::T...) = ZP(enumerate(z)...)
ZP{I<:Integer}(c::Function, i::I...) = ZP(c, enumerate(i)...)

# number methods

zero{T<:ZModN}(::Type{ZPoly{T}}) = ZPoly{T}(T[])
one{T<:ZModN}(::Type{ZPoly{T}}) = ZPoly{T}([one(T)])

# does not create a new copy
convert{T<:ZModN}(::Type{Array{T,1}}, a::ZPoly{T}) = a.a
modulus{T<:ZModN}(::ZPoly{T}) = modulus(T)
modulus{T<:ZModN}(::Type{ZPoly{T}}) = modulus(T)

showcompact{T<:ZModN}(io::IO, p::ZPoly{T}) = showcompact(io, convert(Array{T,1}, p))
function show{T<:ZModN}(io::IO, p::ZPoly{T})
    n = length(p.a)
    if n == 0
       print(io, "0")
    else
        for i in n:-1:1
            if p.a[i] > zero(T)
                if i < n
                    print(io, " + ")
                end
                if p.a[i] > one(T) || i == 1
                    showcompact(io, p.a[i])
                    if i > 1
                        print(io, " ")
                    end
                end
                if i == 2
                    print(io, "x")
                elseif i > 2
                    print(io, "x^$(i-1)")
                end
            end
        end
    end
    print(io, " mod $(modulus(T))")
end

=={T<:ZModN}(a::ZPoly{T}, b::ZPoly{T}) = a.a == b.a
<={T<:ZModN}(a::ZPoly{T}, b::ZPoly{T}) = a.a <= b.a
<{T<:ZModN}(a::ZPoly{T}, b::ZPoly{T}) = a.a < b.a

liftp(c, f, a) = c(f(a.a))
liftp(c, f, a, b) = c(f(a.a, b.a))
-{T<:ZModN}(a::ZPoly{T}) = liftp(ZPoly{T}, -, a)
+{T<:ZModN}(a::ZPoly{T}, b::ZPoly{T}) = liftp(ZPoly{T}, +, a, b)
-{T<:ZModN}(a::ZPoly{T}, b::ZPoly{T}) = liftp(ZPoly{T}, -, a, b)


function test_p_constructor()
    @assert ZPoly{ZField{2,Int}}([GF2(0), GF2(1)]) == ZP(ZF(2), (2, 1))
    @assert ZP(GF2, 0, 0, 0) == zero(ZPoly{ZField{2,Int}})
    @zfield 3 begin
        @assert ZP(1,2,3,4) == ZP(1,2,0,1)
    end
    println("test_p_constructor ok")
end

function test_p_type()
    @assert convert(Array{ZField{2,Int},1}, ZP(GF2, 0, 1)) == [GF2(0), GF2(1)]
    @assert string(ZP(GF2, 1, 0, 1)) == "x^2 + 1 mod 2"
    println("test_p_type ok")
end

function test_p_arithmetic()
    @assert ZP(GF2, (4,1), (1,1)) + ZP(GF2, (3,1), (1,1)) == ZP(GF2, 0, 0, 1, 1)
    println("test_p_arithmetic ok")
end

function tests_p()
    test_p_constructor()
    test_p_type()
    test_p_arithmetic()
end


function tests()
    tests_z()
    tests_p()
end

# run by travis (see .travis.yml in root dir)
#tests()

end
