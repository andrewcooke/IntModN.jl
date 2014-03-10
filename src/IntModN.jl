
# a pragmatic library for doing modular arithmetic.

# the aim here is not to encapsulate a large amount of theory, or to
# describe the theoretical relationships between different structures,
# but to enable arithmetic on various types, motivated largely by the
# practical needs of crypto code.

# currently incomplete; the aim is to support rings and fields (ie
# prime moduli) of integers, polynomials over those, and fields of
# irreducible factor polynomials.

# much thanks to Andreas Noack Jensen
# https://groups.google.com/d/msg/julia-users/Ui977brdzAU/u4rWiDeJv-MJ
# https://gist.github.com/andreasnoackjensen/9132511 (but since this
# changes things significanty, mistakes are mine...)


module IntModN

import Base: show, showcompact, zero, one, inv, real, abs, convert,
       promote_rule, length, getindex, setindex!, start, done, next,
       rand, rand!, print

# Pkg.clone("https://github.com/astrieanna/TypeCheck.jl.git")
#using TypeCheck
# Pkg.clone("https://github.com/vtjnash/Polynomial.jl")
using Polynomial

export ZModN, ZField, ZRing, ZR, ZF, GF2, @zring, @zfield, Poly, X, P,
       order, modulus, factor, itype,
       FModN, FRing, FR



# --- integers modulo n


# we need two types here because prime moduli have a faster inverse
# (via euler's theorem) (note - we do NOT test for primality).

# Real rather than Number because we want < etc to work
abstract ZModN{N, I<:Integer} <: Real

# you can construct these directly, but code assumes that the integer is
# already reduced mod N, that capacities are ok, etc.

immutable ZRing{N, I<:Integer} <: ZModN{N,I}
    i::I
end

# this assumes N is prime and uses a faster inverse
immutable ZField{N, I<:Integer} <: ZModN{N,I}
    i::I
end

# but typically you would use one of these constructors

function validate{I<:Integer}(n::Int, i::I)
    @assert n > 0 "modulus ($n) too small"
    @assert n <= typemax(typeof(i)) "modulus ($n) too large for $(typeof(i))"
    mod(i, n)
end

ZR{I<:Integer}(n::Int, i, ::Type{I}) = ZRing{n,I}(validate(n, convert(I, i)))
ZR(n, i::Integer) = ZRing{n, typeof(i)}(validate(n, i))
ZR(n) = i::Integer -> ZR(n, i)  # used to construct constructors
ZF{I<:Integer}(n::Int, i, ::Type{I}) = ZField{n,I}(validate(n, convert(I, i)))
ZF(n, i::Integer) = ZField{n, typeof(i)}(validate(n, i))
ZF(n) = i::Integer -> ZF(n, i)  # used to construct constructors

GF2 = ZF(2)

modulus{N, I<:Integer}(::Type{ZModN{N, I}}) = N
modulus{T<:ZModN}(::Type{T}) = modulus(super(T))  # jameson type chain
modulus{T<:ZModN}(::T) = modulus(T)
order{T<:ZModN}(::Type{T}) = modulus(T)
order{T<:ZModN}(::T) = modulus(T)
itype{N, I<:Integer}(::Type{ZModN{N, I}}) = I
itype{T<:ZModN}(::Type{T}) = itype(super(T))  # jameson type chain
itype{T<:ZModN}(::T) = itype(T)
duplicate{Z<:ZModN}(::Type{Z}, i::Integer) = Z(validate(modulus(Z), convert(itype(Z), i)))

# all the methods that make these numbers

zero{N,I}(::Type{ZRing{N,I}}) = ZRing{N,I}(zero(I))
one{N,I}(::Type{ZRing{N,I}}) = ZRing{N,I}(one(I))
zero{N,I}(::Type{ZField{N,I}}) = ZField{N,I}(zero(I))
one{N,I}(::Type{ZField{N,I}}) = ZField{N,I}(one(I))

convert(::Type{Bool}, z::ZRing) = convert(Bool, z.i)
convert(::Type{Bool}, z::ZField) = convert(Bool, z.i)
convert{X<:Integer}(::Type{X}, z::ZRing) = convert(X, z.i)
convert{X<:Integer}(::Type{X}, z::ZField) = convert(X, z.i)
convert{N,I<:Integer}(::Type{ZRing{N,I}}, i::I) = ZR(N, i)
convert{N,I<:Integer}(::Type{ZField{N,I}}, i::I) = ZF(N, i)

showcompact{N,I}(io::IO, z::ZModN{N,I}) = showcompact(io, z.i)
print{N,I}(io::IO, z::ZModN{N,I}) = print(io, "$(z.i) mod $N")
show{N,I<:Integer}(io::IO, z::ZField{N,I}) = print(io, "ZField{$N,$I}($(z.i))")
show{N,I<:Integer}(io::IO, z::ZRing{N,I}) = print(io, "ZRing{$N,$I}($(z.i))")

# the random api for ints is kinda broken in that it doesn't take a generator
rand{N,I<:Integer}(T::Type{ZRing{N,I}}) = ZR(N, rand(one(I):convert(I,N)))
rand{N,I<:Integer}(T::Type{ZField{N,I}}) = ZF(N, rand(one(I):convert(I,N)))
function rand!{T<:ZModN}(::Type{T}, A::AbstractArray)
    for i in 1:length(A)
        A[i] = rand(T)
    end
    A
end
rand(::Type{ZModN}, dims...) = error("use more specific type")
rand{T<:ZModN}(::Type{T}, dims...) = rand!(T, Array(T, dims))

# equal N from types, so we don't need to test explicitly
=={N,I}(a::ZModN{N,I}, b::ZModN{N,I}) = convert(I, a) == convert(I, b)
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
    while newr != zero(I)
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


# macros to rewrite literals mod n 
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



# --- polynomials


# unlike Polynomial.jl we don't try to handle reals that are close to
# zero; this code is intended mainly for integers, but if it works
# for reals for you, great (this code used to use Polynomial.jl, but
# it didn't extend Number and had various ther small issues).

# like Polynomial.jl we store the coeffs "backwards"; unlike that
# package we expect the array to be trnucated so that the first value
# is non-zero.

immutable ZPoly{T<:Number} <: Number
    a::Vector{T}
    ZPoly(a) = length(a) == 0 || a[0] > zero(T) ? new(a) : error("zero leading coeff")
end

# constructors

function validate{N}(coeffs::Vector{N})
    start = 1
    while start < length(coeffs) && coeffs[start] == zero(N)
        start += 1
    end
    start == 1 ? coeffs : coeffs[start:end]
end

ZP() = error("provide at least one (zero?) coeff")
ZP(c::Function) = error("provide at least one (zero?) coeff")
ZP(c::Function, coeffs...) = ZP(map(c, coeffs))
ZP(coeffs...) = ZPoly(validate([coeffs...]))

# number methods

zero{T<:ZModN}(::Type{ZPoly{T}}) = ZPoly{T}(T[])
one{T<:ZModN}(::Type{ZPoly{T}}) = ZPoly{T}([one(T)])

# does not create a new copy
convert{T<:ZModN}(::Type{Vector{T}}, a::ZPoly{T}) = a.a
modulus{T<:ZModN}(::ZPoly{T}) = modulus(T)
modulus{T<:ZModN}(::Type{ZPoly{T}}) = modulus(T)

showcompact{T<:ZModN}(io::IO, p::ZPoly{T}) = showcompact(io, convert(Vector{T}, p))
function show{T<:ZModN}(io::IO, p::ZPoly{T})
    n = length(p.a)
    if n == 0
       print(io, "0")
    else
        for i in 1:n
            if p.a[i] != zero(T)
                if i == 1
                    if p.a[i] < zero(T)
                        print(io, "-")
                    end
                else
                    print(io, p.a[i] < zero(T) ? " - " : " + ")
                end
                if p.a[i] != one(T) || i == n
                    showcompact(io, abs(p.a[i]))
                    if i != n
                        print(io, " ")
                    end
                end
                if i == n - 1
                    print(io, "x")
                elseif i < n - 1
                    print(io, "x^$(n-i)")
                end
            end
        end
    end
    print(io, " mod $(modulus(T))")
end

# TODO - these broken?
=={T}(a::ZPoly{T}, b::ZPoly{T}) = a.a == b.a
<={T}(a::ZPoly{T}, b::ZPoly{T}) = a.a <= b.a
<{T}(a::ZPoly{T}, b::ZPoly{T}) = a.a < b.a
length(a::ZPoly) = length(a.a)

# can be accessed and modified as arrays
getindex(a::ZPoly, i) = getindex(a.a, i)
setindex!{T}(a::ZPoly{T}, v::T, i) = setindex!(a, v, i)
start(a::ZPoly) = 1
done(a::ZPoly{T}, i) = i > length(a)
next(a::ZPoly{T}, i) = (a[i], i+1)

# apply function, discarding zeros from start if shrink
function apply{T}(f, a::Array{T,1}, b::Array{T,1}; shrink=true)
    (big, small) = length(a) > length(b) ? (a, b) : (b, a)
    result = copy(big)
    d = length(big) - length(small)
    shrink = shrink && d == 0
    for (i, z) in enumerate(small)
        x = f(big[i+d], z)
        if x != zero(T) && shrink
            result = result[i:end]  # d == 0
            shrink = false
        end
        if !shrink
            result[i+d] = x
        end
    end
    if shrink
        result = T[]
    end
    result
end    

liftp(ca, f, a) = c(f(a.a))
liftp(c, f, a, b) = c(apply(f, a.a, b.a))
-{T}(a::ZPoly{T}) = liftp(ZPoly{T}, -, a)
+{T}(a::ZPoly{T}, b::ZPoly{T}) = liftp(ZPoly{T}, +, a, b)
-{T}(a::ZPoly{T}, b::ZPoly{T}) = liftp(ZPoly{T}, -, a, b)

function *{T}(a::ZPoly{T}, b::ZPoly{T})
    (big, small) = length(a) > length(b) ? (a, b) : (b, a)
    if length(small) == 0
        small
    elseif length(small) == 1 && small[1] == one(T)
        big
    else
        result = zeros(T, length(big) + length(small) - 1)
        for (i, z) in enumerate(small)
            if z != zero(T)
                for (j, w) in enumerate(big)
                    result[i+j-1] += w * z
                end
            end
        end
        ZP(result)  # truncates
    end
end

function divmod{T<:ZModN}(a::ZPoly{T}, b::ZPoly{T})
    @assert length(b) > 0 "division by zero polynomial"
    if length(a) == 0
        (a, a)
    elseif length(b) == 1 && b[1] == one(T)
        (a, zero(T))
    elseif length(b) > length(a)
        (zero(T), a)
    else 
        # TODO
        result = zeros(T, length(big) + length(small) - 1)
        for (i, z) in enumerate(small)
            if z != zero(T)
                for j in 1:length(big)
                    result[i+j-1] += big[j] * z
                end
            end
        end
        if result[end] == zero(T)
            i = length(result)
            while i > 0 && result[i] == zero(T)
                i -= 1
            end
            result = result[1:i]
        end
        ZPoly{T}(result)
    end
end









# --- polynomials over integers modulo n


# instead of re-inventing the wheel we extend the existing Polynomial module
# (note that these store in an array from high to low)

# nice syntax for creating polynomials.
# use x=X(ZF(2)) or x=X(GF2) or x=X(ZField{2,Int}) and then x^2+3 etc
X(c::Function) = Poly([c(1), c(0)])
X{Z<:ZModN}(::Type{Z}) = Poly([one(Z), zero(Z)])

# a compact way of generating from an array and a type
# (we can then use this in show to give a compact but accurate reprn)
# (note that we always use x (the default) as the variable in Poly)
P(c::Function, a::Vector) = Poly(map(c, a))
P(c::Function, a...) = P(c, [a...])
P{N,I<:Integer}(::Type{ZField{N,I}}, a::Vector{I}) = P(ZF(N), a)
P{N,I<:Integer}(::Type{ZField{N,I}}, a::I...) = P(ZF(N), [a...])
P{N,I<:Integer}(::Type{ZRing{N,I}}, a::Vector{I}) = P(ZR(N), a)
P{N,I<:Integer}(::Type{ZRing{N,I}}, a::I...) = P(ZR(N), [a...])

# the defaults repeat "mod n" all over the place
function print{Z<:ZModN}(io::IO, p::Poly{Z})
    print_no_mod(io, p)
    print(io, " mod $(modulus(Z))")
end

function print_no_mod{Z<:ZModN}(io::IO, p::Poly{Z})
    n = length(p)
    if n <= 0
        print(io,"0")
    else
        for j = 1:n
            pj = p[j]
            if pj != zero(Z)
                if j == 1 
                    pj < 0 && print(io, "-")
                else
                    pj < 0 ? print(io," - ") : print(io," + ")
                end
                magpj = abs(pj)
                if j == n || magpj != one(Z)
                    print(io, magpj.i)
                end
                exp = n-j
                if exp > 0
                    print(io, p.var)
                    if exp > 1
                        print(io, '^', exp)
                    end
                end
            end
        end
    end
end

function show{Z<:ZModN}(io::IO, p::Poly{Z})
    print(io, "P($Z")
    for i in 1:length(p)
        print(io, ",$(p[i].i)")
    end
    if p.var != :x
        print(io, ";var=:$(p.var)")
    end
    print(io, ")")
end

# needed below to squeeze polynomial into type
poly_to_tuple{Z<:ZModN}(p::Poly{Z}) = tuple(map(a -> a.i, p.a)...)
tuple_to_poly{Z<:ZModN}(::Type{Z}, t::Tuple) = Poly([map(i -> duplicate(Z, i), t)...])



# --- factor rings and fields


abstract FModN{Z<:ZModN, F} <: Real

immutable FRing{Z<:ZModN, F<:Tuple} <: FModN{Z,F}
    p::Poly{Z}
end

function validate{Z<:ZModN}(p::Poly{Z}, f::Poly{Z})
    a, b = divrem(p, f)
    b
end

FR{Z<:ZModN}(p::Poly{Z}, f::Poly{Z}) = FRing{Z, poly_to_tuple(f)}(validate(p, f))

factor{Z<:ZModN, F<:Tuple}(::Type{FRing{Z, F}}) = tuple_to_poly(Z, F)
factor{Z<:ZModN, F<:Tuple}(::FRing{Z, F}) = factor(FRing{Z, F})
modulus{Z<:ZModN, F}(::Type{FModN{Z,F}}) = modulus(Z)
modulus{F<:FModN}(::F) = modulus(F)
order{F<:FModN}(::Type{F}) = modulus(F) ^ Polynomial.deg(factor(F))
order{F<:FModN}(::F) = order(F)

zero{Z<:ZModN, F<:Tuple}(::Type{FRing{Z, F}}) = FRing(Poly(zero(Z)), tuple_to_poly(Z, F))
one{Z<:ZModN, F<:Tuple}(::Type{FRing{Z, F}}) = FRing(Poly(one(Z)), tuple_to_poly(Z, F))

function show(io::IO, r::FRing)
    print(io, "FR(")
    show(io, r.p)
    print(io, ",")
    show(io, factor(r))
    print(io, ")")
end

function print(io::IO, r::FRing)
    print_no_mod(io, r.p)
    print(io, " mod ");
    print(io, factor(r))
end

liftf{F<:FRing}(f, a::F, b::F) = FR(f(a.p, b.p), factor(F))
+{F<:FRing}(a::F, b::F) = liftf(+, a, b)
-{F<:FRing}(a::F, b::F) = liftf(-, a, b)
*{F<:FRing}(a::F, b::F) = liftf(*, a, b)



# --- pull in tests (does this need ot be so ugly?)

d = Pkg.dir("IntModN")
d = "$d/src"
push!(LOAD_PATH, d)
import Tests: tests
pop!(LOAD_PATH)

end
