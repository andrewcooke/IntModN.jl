
# a pragmatic library for doing modular arithmetic.

# the aim here is not to encapsulate a large amount of theory, or to
# describe the theoretical relationships between different structures,
# but to enable arithmetic on various types, motivated largely by the
# practical needs of crypto code.

# incomplete; pull requests welcome

# much thanks to Andreas Noack Jensen
# https://groups.google.com/d/msg/julia-users/Ui977brdzAU/u4rWiDeJv-MJ
# https://gist.github.com/andreasnoackjensen/9132511 (but since this
# changes things significanty, mistakes are mine...)


module IntModN

import Base: show, showcompact, zero, one, inv, real, abs, convert,
       promote_rule, length, getindex, setindex!, start, done, next,
       rand, rand!, print, map

# Pkg.clone("https://github.com/astrieanna/TypeCheck.jl.git")
#using TypeCheck
# Pkg.clone("https://github.com/vtjnash/Polynomial.jl")
#using Polynomial

export ZModN, ZField, ZRing, ZR, ZF, GF2, @zring, @zfield, 
       ZPoly, ZP, X, P, order, modulus, factor, itype,
       GF2Poly, GF2P, GF2X,
       FModN, FRing, FR,
       extended_euclidean


# we generally support promotion into tyoes here (from integers) and
# conversion back.  but that's all.  users can add other promotion /
# conversion if required.


# using Integer as the base type here as: (1) some numeric type seems
# to be necessary to access many generic numerical functions; (2) need
# Real to use comparisons; (3) have integer-like accuracy behaviour
# (no need to worry about rounding etc).

# but we need an additional type so that we can intercept things like
# automatic promotion to floats from integers.

abstract IntegerOnly <: Integer

/{I<:Integer,M<:IntegerOnly}(i::I, m::M) = convert(M, i) / m
/{I<:Integer,M<:IntegerOnly}(m::M, i::I) = m / convert(M, i)
inv{M<:IntegerOnly}(m::M) = error("no inverse defined")


# --- integers modulo n


abstract ZModN{N, I<:Integer} <: IntegerOnly

# we need two types here because prime moduli have a faster inverse
# (via euler's theorem) (note - we do NOT test for primality).

# you can construct these directly, but code assumes that the integer is
# already reduced mod N, that capacities are ok, etc.

immutable ZRing{N, I<:Integer} <: ZModN{N,I}
    i::I
end

# this assumes N is prime and uses a faster inverse
immutable ZField{N, I<:Integer} <: ZModN{N,I}
    i::I
end

# but typically you would use one of these constructors (ZF or ZR)

function prepare_z{I<:Integer}(n::Int, i::I)
    @assert n > 0 "modulus ($n) too small"
    @assert n <= typemax(typeof(i)) "modulus ($n) too large for $(typeof(i))"
    mod(i, n)
end

for (f, Z) in ((:ZR, :ZRing), (:ZF, :ZField))
    @eval $f{I<:Integer}(n::Int, i, ::Type{I}) = $Z{n,I}(prepare_z(n, convert(I, i)))
    @eval $f(n, i::Integer) = $Z{n, typeof(i)}(prepare_z(n, i))
    @eval $f(n) = i::Integer -> $f(n, i)  # used to construct constructors
end
GF2 = ZF(2)

modulus{N, I<:Integer}(::Type{ZModN{N, I}}) = N
modulus{T<:ZModN}(::Type{T}) = modulus(super(T))  # jameson type chain
modulus{T<:ZModN}(::T) = modulus(T)
order{T<:ZModN}(::Type{T}) = modulus(T)
order{T<:ZModN}(::T) = modulus(T)
itype{N, I<:Integer}(::Type{ZModN{N, I}}) = I
itype{T<:ZModN}(::Type{T}) = itype(super(T))  # jameson type chain
itype{T<:ZModN}(::T) = itype(T)

for (f, Z) in ((:ZR, :ZRing), (:ZF, ZField))
    # julia's comversion is pretty damn dumb at times
    @eval convert{N,I<:Integer}(::Type{$Z{N,I}}, a::$Z{N,I}) = a
    # this promotion needed because linalg code has direct compariosn with 0
    @eval promote_rule{N, I<:Integer}(::Type{$Z{N,I}}, ::Type{Int}) = $Z{N,I}
    @eval convert{N,I<:Integer}(::Type{$Z{N,I}}, i::Integer) = $f(N, convert(I, i))
end
# conversion to int, but no promotion
for T in (:Int, :Uint)
    @eval convert(::Type{$T}, z::ZModN) = convert($T, z.i)
end

# all the methods that make these numbers

for Z in (:ZRing, :ZField)
    @eval zero{N,I}(::Type{$Z{N,I}}) = $Z{N,I}(zero(I))
    @eval one{N,I}(::Type{$Z{N,I}}) = $Z{N,I}(one(I))
end
zero{Z<:ZModN}(z::Z) = zero(Z)
one{Z<:ZModN}(z::Z) = one(Z)

showcompact{N,I}(io::IO, z::ZModN{N,I}) = showcompact(io, z.i)
print{N,I}(io::IO, z::ZModN{N,I}) = print(io, "$(z.i) mod $N")
show{N,I<:Integer}(io::IO, z::ZField{N,I}) = print(io, "ZField{$N,$I}($(z.i))")
show{N,I<:Integer}(io::IO, z::ZRing{N,I}) = print(io, "ZRing{$N,$I}($(z.i))")

# the random api for ints is kinda broken in that it doesn't take a generator
rand{N,I<:Integer}(T::Type{ZRing{N,I}}) = ZR(N, rand(one(I):convert(I,N)))
rand{N,I<:Integer}(T::Type{ZField{N,I}}) = ZF(N, rand(one(I):convert(I,N)))
function rand!{T<:ZModN}(::Type{T}, a::AbstractArray)
    for i in 1:length(a)
        a[i] = rand(T)
    end
    a
end
rand(::Type{ZModN}, dims...) = error("use more specific type")
rand{T<:ZModN}(::Type{T}, dims...) = rand!(T, Array(T, dims))

# equal N from types, so we don't need to test explicitly
=={N,I}(a::ZModN{N,I}, b::ZModN{N,I}) = a.i == b.i
<={N,I}(a::ZModN{N,I}, b::ZModN{N,I}) = a.i <= b.i
<{N,I}(a::ZModN{N,I}, b::ZModN{N,I}) = a.i < b.i

real{N,I}(a::ZModN{N,I}) = a
abs{N,I}(a::ZModN{N,I}) = a

liftz(c, n, i, f, a) = c(convert(i, mod(f(a.i), n)))
liftz(c, n, i, f, a, b) = c(convert(i, mod(f(a.i, b.i), n)))
for Z in (ZRing, ZField)    
    @eval -{N,I}(a::$Z{N,I}) = liftz($Z{N,I}, N, I, -, a)
    @eval ^{N,I}(a::$Z{N,I}, p::Int) = liftz($Z{N,I}, N, I, x -> powermod(x, p, N), a)
    for op in (:+, :-, :*, :&, :|, :$)
        @eval $op{N,I}(a::$Z{N,I}, b::$Z{N,I}) = liftz($Z{N,I}, N, I, $op, a, b)
    end
end

# http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#IntegerOnly_integers
function extended_euclidean{I<:Integer}(a::I, n::I)
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

inv{N,I}(a::ZRing{N,I}) = ZRing{N,I}(extended_euclidean(N, convert(I, N)))
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



# --- arbitrary (but "integer") coeff polynomials


# started with Polynomial.jl, but that has various issues; this is
# less general (and faster and uses less memory - see PTests.jl).

# IMPORTANT - arrays may be shared between polynomials, so don't mutate
# contents unless you have a new instance whose contents you know to be
# owned.

immutable ZPoly{I<:Integer} <: IntegerOnly
    a::Vector{I}
    ZPoly(a) = length(a) == 0 || a[1] > zero(I) ? new(a) : error("zero leading coeff")
end

# constructors

function prepare_p{N}(coeffs::Vector{N})
    start = 1
    while start <= length(coeffs) && coeffs[start] == zero(N)
        start += 1
    end
    start == 1 ? coeffs : coeffs[start:end]
end

ZP() = error("provide at least one (zero?) coeff")
ZP(c::Function) = error("provide at least one (zero?) coeff")
ZP(t::Type, coeffs::Vector) = ZP(t[map(c -> convert(t, c), coeffs)...])
ZP(t::Type, coeffs...) = ZP(t[map(c -> convert(t, c), coeffs)...])
ZP(c::Function, coeffs::Vector) = ZP(map(c, coeffs))
ZP(c::Function, coeffs...) = ZP([map(c, coeffs)...])
ZP{I<:Integer}(coeffs::Vector{I}) = ZPoly{I}(prepare_p(coeffs))
ZP{I<:Integer}(coeffs::I...) = ZP([coeffs...])

# needed below to squeeze polynomial into type
# not via convert as we don't know type (length of tuple).  
# TODO - this may cause efficiency issues (since return type varies)?
poly_to_tuple{I<:Integer}(p::ZPoly{I}) = tuple(map(a -> convert(Int, a), p.a)...)
tuple_to_poly{I<:Integer}(::Type{I}, t::Tuple) = ZPoly{I}([map(i -> convert(I, i), t)...])

# nice syntax for creating polynomials.
# use x=X(ZF(2)) or x=X(GF2) or x=X(ZField{2,Int}) and then x^2+3 etc
X(c::Function) = ZP([c(1), c(0)])
X{I<:Integer}(::Type{I}) = ZP([one(I), zero(I)])

modulus{T}(::ZPoly{T}) = modulus(T)
modulus{T}(::Type{ZPoly{T}}) = modulus(T)

# this promotion to write x^2+1 etc
convert{I<:Integer}(::Type{ZPoly{I}}, p::ZPoly{I}) = p
promote_rule{I<:Integer}(::Type{ZPoly{I}}, ::Type{Int}) = ZPoly{I}
convert{I<:Integer}(::Type{ZPoly{I}}, i::Integer) = ZP([convert(I, i)])

# can be accessed and modified (but don't do that in public) as arrays
getindex(a::ZPoly, i) = getindex(a.a, i)
setindex!{T}(a::ZPoly{T}, v::T, i) = setindex!(a.a, v, i)
length(a::ZPoly) = length(a.a)
endof(a::ZPoly) = length(a)
start(::ZPoly) = 1
done(a::ZPoly, i) = i > length(a)
next(a::ZPoly, i) = (a[i], i+1)

show_int{N,I<:Integer}(io::IO, z::ZModN{N,I}) = show(io, convert(I, z))
show_int{I<:Integer}(io::IO, i::I) = show(io, convert(I, i))
print_mode{Z<:ZModN}(io::IO, p::ZPoly{Z}) = print(io, " mod $(modulus(Z))")
print_mode{I<:Integer}(io::IO, p::ZPoly{I}) = nothing

function print_no_mod{I<:Integer}(io::IO, p::ZPoly{I})
    n = length(p)
    if n <= 0
        show_int(io, zero(I))
    else
        for j = 1:n
            pj = p[j]
            if pj != zero(I)
                if j == 1 
                    pj < 0 && print(io, "-")
                else
                    pj < 0 ? print(io," - ") : print(io," + ")
                end
                magpj = abs(pj)
                if j == n || magpj != one(I)
                    show_int(io, magpj)
                end
                exp = n-j
                if exp > 0
                    print(io, :x)
                    if exp > 1
                        print(io, '^', exp)
                    end
                end
            end
        end
    end
end

function print{I<:Integer}(io::IO, p::ZPoly{I})
    print_no_mod(io, p)
    print_mode(io, p)
end

function show{I<:Integer}(io::IO, p::ZPoly{I})
    print(io, "ZP($I")
    for i in 1:length(p)
        print(io, ",")
        show_int(io, p[i])
    end
    print(io, ")")
end

showcompact(io::IO, p::ZPoly) = showcompact(io, p.a)

# number methods

zero{T}(::Type{ZPoly{T}}) = ZPoly{T}(T[])
one{T}(::Type{ZPoly{T}}) = ZPoly{T}([one(T)])
zero{Z<:ZPoly}(z::Z) = zero(Z)
one{Z<:ZPoly}(z::Z) = one(Z)

=={T}(a::ZPoly{T}, b::ZPoly{T}) = a.a == b.a
function cmp{T}(a::ZPoly{T}, b::ZPoly{T}, eq)
    if length(a) < length(b)
        true
    elseif length(a) > length(b)
        false
    else
        for i in 1:length(a)
            if a[i] < b[i]
                return true
            elseif a[i] > b[i]
                return false
            end
        end
        eq
    end
end
<={T}(a::ZPoly{T}, b::ZPoly{T}) = cmp(a, b, true)
<{T}(a::ZPoly{T}, b::ZPoly{T}) = cmp(a, b, false)

# some optimisation here to reduce copying of arrays
# big is modified; small is not - see liftp and +,- below
# apply function, discarding zeros from start if shrink
function apply{T}(f, big::Array{T,1}, small::Array{T,1})
    shift = length(big) - length(small)
    @assert shift >= 0
    shrink = shift == 0
    for i in 1:length(small)
        x = f(big[i+shift], small[i])
        if x != zero(T) && shrink
            big = big[i:end]  # shift == 0
            shift = 1-i
            shrink = false
        end
        if !shrink
            big[i+shift] = x
        end
    end
    if shrink
        big = T[]
    end
    big
end    

liftp(c, f, a) = c(f(a.a))
liftp(c, f, aa, ba) = c(apply(f, aa, ba))
-{T}(a::ZPoly{T}) = liftp(ZPoly{T}, -, a)
+{T}(a::ZPoly{T}, b::ZPoly{T}) = length(a) >= length(b) ? liftp(ZPoly{T}, +, copy(a.a), b.a) : liftp(ZPoly{T}, +, copy(b.a), a.a)
-{T}(a::ZPoly{T}, b::ZPoly{T}) = length(a) >= length(b) ? liftp(ZPoly{T}, -, copy(a.a), b.a) : liftp(ZPoly{T}, +, -b.a, a.a)

function *{T}(a::ZPoly{T}, b::ZPoly{T})
    big, small = length(a) > length(b) ? (a, b) : (b, a)
    if length(small) == 0  # no need to test big
        small
    elseif length(small) == 1 && small[1] == one(T)
        big
    else
        result = zeros(T, length(big) + length(small) - 1)
        for i in 1:length(small)
            z = small[i]  # enumerate here uses strange amount of memory
            if z != zero(T)
                for j in 1:length(big)
                    result[i+j-1] += big[j] * z
                end
            end
        end
        ZP(result)
    end
end

function _divrem{T}(a::Vector{T}, b::Vector{T})
    @assert length(b) > 0 "division by zero polynomial"
    if length(a) == 0
        (a, a)
    elseif length(b) == 1 && b[1] == one(T)
        (a, T[])
    elseif length(b) > length(a)
        (T[], a)
    elseif a == b
        ([one(T)], T[])
    else 
        shift = length(a) - length(b)
        rem, div = copy(a), zeros(T, shift + 1)
        for s in 0:shift
            factor = rem[s+1] / b[1]
            div[s+1] = factor
            for i in 1:length(b)
                rem[s+i] -= factor * b[i]
            end
        end
        div, rem
    end
end

function divrem{T}(a::ZPoly{T}, b::ZPoly{T})
    map(ZP, _divrem(a.a, b.a))
end

/{T}(a::ZPoly{T}, b::ZPoly{T}) = ZP(_divrem(a.a, b.a)[1])
%{T}(a::ZPoly{T}, b::ZPoly{T}) = ZP(_divrem(a.a, b.a)[2])
# ^ will come from power_by_squaring in stdlib
degree{T}(a::ZPoly{T}) = length(a) - 1

# defining bitwise ops seems odd, but makes sense for GF2, especially
# when using the compact binary repn below.  for N>2 we apply the
# operation to each coefficient pair in turn.

function same_length{T}(a::Vector{T}, b::Vector{T})
    big, small = length(a) > length(b) ? (a, b) : (b, a)
    shift = length(big) - length(small)
    small = vcat(zeros(T, shift), small)
    big, small
end

map{P<:ZPoly}(f::Function, a::P, b::P) = ZP(map(f, same_length(a.a, b.a)...))
for op in (:&, :|, :$)
    @eval $op{P<:ZPoly}(a::P, b::P) = map($op, a, b)
end
<<{T}(a::ZPoly{T}, n::Int) = a * X(T)^n
>>>{T}(a::ZPoly{T}, n::Int) = a / X(T)^n



# --- GF(2) polynomials encoded as bits


type GF2Poly{I<:Unsigned} <: IntegerOnly
    i::I
end

# constructors
GF2P{U<:Unsigned}(p::U) = GF2Poly{U}(p)
GF2P{S<:Signed}(p::S) = GF2P(unsigned(p))
GF2X{U<:Unsigned}(::Type{U}) = GF2Poly{U}(one(U))

modulus(::Type{GF2Poly}) = 2
modulus(::GF2Poly) = 2

convert{I<:Integer, U<:Unsigned}(::Type{GF2Poly{U}}, i::I) = GF2P(convert(U, i % 2))
promote_rule{I<:Integer, U<:Unsigned}(::Type{GF2Poly{U}}, ::Type{I}) = GF2Poly{U}

# TODO - getindex etc?

# for display we first convert to to ZPoly.  this is not efficient,
# but we probably need theinterop anyway, and who cares if printing is
# efficient?

convert{I<:Integer, U<:Unsigned}(::Type{ZPoly{ZField{2,I}}}, p::GF2Poly{U}) = _convert(ZPoly{ZField{2,I}}, p)
convert{I<:Integer, U<:Unsigned}(::Type{GF2Poly{U}}, p::ZPoly{ZField{2,I}}) = _convert(GF2Poly{U}, p)

function _convert{T, F}(::Type{T}, from::F)
    result, mask, x = zero(T), one(T), one(T) << 1
    none = zero(F)
    while from != none
        println("from ", from)
        println("mask ", mask)
        if from & 1 != none
            result |= mask
        end
        mask <<= 1
        from >>>= 1
    end
    result
end

for (name, op) in ((:+, $), (:-, $), (:|, |), (:&, &))
    @eval $name{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U}) = GF2Poly($op(a.i, b.i))
end
<<{U<:Unsigned}(a::GF2Poly{U}, n::Int) = GF2Poly(a.i << n)
>>>{U<:Unsigned}(a::GF2Poly{U}, n::Int) = GF2Poly(a.i >>> n)

degree{U<:Unsigned}(p::GF2Poly{U}) = 8*sizeof(U) - leading_zeros(p)



# --- factor rings and fields


# this all assumes that the factor poynomial is irreducible

abstract FModN{Z<:ZModN, F} <: IntegerOnly

# assumes already reduced
immutable FRing{Z<:ZModN, F<:Tuple} <: FModN{Z,F}
    p::ZPoly{Z}
end

function prepare_f{Z<:ZModN}(p::ZPoly{Z}, f::ZPoly{Z})
    p % f
end

FR{Z<:ZModN}(p::ZPoly{Z}, f::ZPoly{Z}) = FRing{Z, poly_to_tuple(f)}(prepare_f(p, f))

factor{Z<:ZModN, F<:Tuple}(::Type{FRing{Z, F}}) = tuple_to_poly(Z, F)
factor{Z<:ZModN, F<:Tuple}(::FRing{Z, F}) = factor(FRing{Z, F})
modulus{Z<:ZModN, F}(::Type{FModN{Z,F}}) = modulus(Z)
modulus{F<:FModN}(::F) = modulus(F)
# assuming irreducible factor
order{F<:FModN}(::Type{F}) = modulus(F) ^ degree(factor(F)) - 1
order{F<:FModN}(::F) = order(F)

zero{Z<:ZModN, F<:Tuple}(::Type{FRing{Z, F}}) = FR(ZP(zero(Z)), tuple_to_poly(Z, F))
one{Z<:ZModN, F<:Tuple}(::Type{FRing{Z, F}}) = FR(ZP(one(Z)), tuple_to_poly(Z, F))
zero{F<:FModN}(f::F) = zero(F)
one{F<:FModN}(f::F) = one(F)

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

=={F<:FRing}(a::F, b::F) = a.p == b.p
<={F<:FRing}(a::F, b::F) = a.p <= b.p
<{F<:FRing}(a::F, b::F) = a.p < b.p

liftf{F<:FRing}(f, a::F, b::F) = FR(f(a.p, b.p), factor(F))
+{F<:FRing}(a::F, b::F) = liftf(+, a, b)
-{F<:FRing}(a::F, b::F) = liftf(-, a, b)
*{F<:FRing}(a::F, b::F) = liftf(*, a, b)

inv{F<:FRing}(f::F) = FR(extended_euclidean(f.p, factor(f)), factor(f))
/{F<:FRing}(a::F, b::F) = a * inv(b)



# --- pull in tests (does this need ot be so ugly?)


# TODO - try include (relative path) here

d = Pkg.dir("IntModN")
d = "$d/src"
push!(LOAD_PATH, d)
import Tests: tests
pop!(LOAD_PATH)

end
