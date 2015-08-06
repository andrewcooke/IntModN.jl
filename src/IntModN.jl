

# TODO - overflow
# TODO - fast mult as http://maths-people.anu.edu.au/~brent/pd/rpb232.pdf

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
       rand, rand!, print, map, leading_zeros, divrem, endof, bits,
       div, rem, mod, trailing_zeros, /, ==, <=, <, -, ^, +, *, &,
       |, $, %, <<, >>>, >>
import Base.Random.AbstractRNG

using AutoTypeParameters 

export ZModN, ZField, ZRing, ZR, ZF, GF2, @zring, @zfield, 
       ZPoly, ZP, X, P, order, modulus, factor, itype,
       GF2Poly, GF2P, GF2X,
       FModN, FRing, FR,
       extended_euclidean


# we generally support promotion into types here (from integers) and
# conversion back.  but that's all.  users can add other promotion /
# conversion if required.


# using Integer as the base type here as: (1) some numeric type seems
# to be necessary to access many generic numerical functions; (2) need
# Real to use comparisons; (3) have integer-like accuracy behaviour
# (no need to worry about rounding etc).

# but we need an additional type so that we can intercept things like
# automatic promotion to floats from integers.

abstract Residue <: Integer

/{I<:Integer,M<:Residue}(i::I, m::M) = convert(M, i) / m
/{I<:Integer,M<:Residue}(m::M, i::I) = m / convert(M, i)
inv{M<:Residue}(m::M) = error("no inverse defined")


# --- integers modulo n

# J is the type of N.  this may seem superfluous, but N itself is frozen, so its type is not always available (see AutoTypeParameters)

abstract ZModN{J<:Integer, N, I<:Integer} <: Residue

# we need two types here because prime moduli have a faster inverse
# (via euler's theorem) (note - we do NOT test for primality).

# you can construct these directly, but code assumes that the integer is
# already reduced mod N, that capacities are ok, etc.

immutable ZRing{J<:Integer, N, I<:Integer} <: ZModN{J,N,I}
    i::I
    ZRing(i) = (@assert isa(i, J); new(i))
end

# this assumes N is prime and uses a faster inverse
immutable ZField{J<:Integer, N, I<:Integer} <: ZModN{J,N,I}
    i::I
    ZField(i) = (@assert isa(i, J); new(i))
end

# but typically you would use one of these constructors (ZF or ZR)

function prepare_z{I<:Integer}(n::Int, i::I)
    n > 0 || error("modulus ($n) too small")
    n <= typemax(typeof(i)) || error("modulus ($n) too large for $(typeof(i))")
    mod(i, n)
end

for (f, Z) in ((:ZR, :ZRing), (:ZF, :ZField))
    @eval $f{J<:Integer, I<:Integer}(n, ::Type{J}, i, ::Type{I}) = $Z{J, freeze(convert(J, n)), I}(prepare_z(n, convert(I, i)))
    @eval $f{I<:Integer}(n::Int, i, ::Type{I}) = $Z{typeof(n), freeze(n), I}(prepare_z(n, convert(I, i)))
    @eval $f(n::Int, i::Integer) = $Z{typeof(n), freeze(n), typeof(i)}(prepare_z(n, i))
    @eval $f(n::Integer, i::Integer) = $Z{typeof(n), freeze(n), typeof(i)}(mod(i, n))
    @eval $f(n) = i::Integer -> $f(n, i)  # used to construct constructors
end
GF2 = ZF(2)

modulus{J<:Integer, N, I<:Integer}(::Type{ZModN{J, N, I}}) = convert(J, thaw(eval, N))
modulus{T<:ZModN}(::Type{T}) = modulus(super(T))  # jameson type chaina
modulus{T<:ZModN}(::T) = modulus(T)
order{T<:ZModN}(::Type{T}) = modulus(T)
order{T<:ZModN}(::T) = modulus(T)
itype{J, N, I<:Integer}(::Type{ZModN{J, N, I}}) = I
itype{T<:ZModN}(::Type{T}) = itype(super(T))  # jameson type chain
itype{T<:ZModN}(::T) = itype(T)
ntype{J, N, I<:Integer}(::Type{ZModN{J, N, I}}) = J
ntype{T<:ZModN}(::Type{T}) = ntype(super(T))  # jameson type chain
ntype{T<:ZModN}(::T) = ntype(T)

for (f, Z) in ((:ZR, :ZRing), (:ZF, ZField))
    # julia's comversion is pretty damn dumb at times
    @eval convert{J<:Integer,N,I<:Integer}(::Type{$Z{J,N,I}}, a::$Z{J,N,I}) = a
    # this promotion needed because linalg code has direct compariosn with 0
    @eval promote_rule{J<:Integer, N, I<:Integer}(::Type{$Z{J,N,I}}, ::Type{Int}) = $Z{J,N,I}
    @eval convert{J<:Integer,N,I<:Integer}(::Type{$Z{J,N,I}}, i::Integer) = $f(convert(J, N), convert(I, i))
end
# conversion to int, but no promotion
for T in (:Int, :Uint)
    @eval convert(::Type{$T}, z::ZModN) = convert($T, z.i)
end

# all the methods that make these numbers

for Z in (:ZRing, :ZField)
    @eval zero{J,N,I}(::Type{$Z{J,N,I}}) = $Z{J,N,I}(zero(I))
    @eval one{J,N,I}(::Type{$Z{J,N,I}}) = $Z{J,N,I}(one(I))
end
zero{Z<:ZModN}(z::Z) = zero(Z)
one{Z<:ZModN}(z::Z) = one(Z)

showcompact{J,N,I}(io::IO, z::ZModN{J,N,I}) = showcompact(io, z.i)
print{J,N,I}(io::IO, z::ZModN{J,N,I}) = print(io, "$(z.i) mod $N")
show{J,N,I<:Integer}(io::IO, z::ZField{J,N,I}) = print(io, "ZField{$J,$N,$I}($(z.i))")
show{J,N,I<:Integer}(io::IO, z::ZRing{J,N,I}) = print(io, "ZRing{$J,$N,$I}($(z.i))")

rand_{J,N,I<:Integer}(T::Type{ZRing{J,N,I}}) = ZR(N, J, rand(one(I):convert(I,modulus(T))), I)
rand_{J,N,I<:Integer}(T::Type{ZField{J,N,I}}) = ZF(N, J, rand(one(I):convert(I,modulus(T))), I)
function rand!{T<:ZModN}(::Type{T}, a::AbstractArray)
    for i in 1:length(a)
        a[i] = rand_(T)
    end
    a
end

rand(::Type{ZModN}, dims...) = error("use more specific type")

function rand{T<:ZModN}(::Type{T}, dims...)
    if length(dims) > 0
        rand!(T, Array(T, dims))
    else
        rand_(T)
    end
end

if VERSION >= v"0.4-"
    rand{J<:Integer,N,I<:Integer}(rng::AbstractRNG, T::Type{ZRing{J,N,I}}) = ZR(N, J, rand(rng, one(I):convert(I,N)), I)
    rand{J<:Integer,N,I<:Integer}(rng::AbstractRNG, T::Type{ZField{J,N,I}}) = ZF(N, J, rand(rng, one(I):convert(I,N)), I)
    function rand!{T<:ZModN}(rng::AbstractRNG, ::Type{T}, a::AbstractArray)
        for i in 1:length(a)
            a[i] = rand(rng, T)
        end
        a
    end
    rand(rng::AbstractRNG, ::Type{ZModN}, dims::Int...) = 
    error("use more specific type")
    rand{T<:ZModN}(rng::AbstractRNG, ::Type{T}, dims::Int...) = 
    rand!(rng, T, Array(T, dims))    
end

# equal N from types, so we don't need to test explicitly
=={T,N,I}(a::ZModN{T,N,I}, b::ZModN{T,N,I}) = a.i == b.i
<={T,N,I}(a::ZModN{T,N,I}, b::ZModN{T,N,I}) = a.i <= b.i
<{T,N,I}(a::ZModN{T,N,I}, b::ZModN{T,N,I}) = a.i < b.i

real{T,N,I}(a::ZModN{T,N,I}) = a
abs{T,N,I}(a::ZModN{T,N,I}) = a

liftz(c, i, f, a) = c(convert(i, mod(f(a.i), modulus(a))))
liftz(c, i, f, a, b) = c(convert(i, mod(f(a.i, b.i), modulus(a))))
for Z in (ZRing, ZField)    
    @eval -{J,N,I}(a::$Z{J,N,I}) = liftz($Z{J,N,I}, I, -, a)
    @eval ^{J,N,I}(a::$Z{J,N,I}, p::Int) = liftz($Z{J,N,I}, I, x -> powermod(x, p, modulus(a)), a)
    for op in (:+, :-, :*, :&, :|, :$)
        @eval $op{J,N,I}(a::$Z{J,N,I}, b::$Z{J,N,I}) = liftz($Z{J,N,I}, I, $op, a, b)
    end
end

# http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Residue_integers
function extended_euclidean{I<:Integer}(a::I, n::I)
    t::I, newt::I = zero(I), one(I)
    r::I, newr::I = n, a
    while newr != zero(I)
        q::I = div(r, newr)
        t, newt = newt, t - q * newt
        r, newr = newr, r - q * newr
    end
    r <= 1 || error("$a is not invertible mod $n")
    t < 0 ? t + n : t
end

inv{J,N,I}(a::ZRing{J,N,I}) = ZRing{J,N,I}(extended_euclidean(a.i, convert(I, modulus(a))))
inv{J,N,I}(a::ZField{J,N,I}) = a^(modulus(a)-2)

function /{T,N,I}(a::ZModN{T,N,I}, b::ZModN{T,N,I})
    if a == b
        one(a)
    elseif a == zero(a)
        a
    elseif b == zero(b)
        error("division by zero: $a/$b")
    else
        a * inv(b)
    end
end

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



# --- supertype for all integer polynomial representations

# I is the type that [] returns when indexing a coefficient
abstract Poly{I<:Integer} <: Residue

endof(a::Poly) = length(a)
start(::Poly) = 1
done(a::Poly, i) = i > length(a)
next(a::Poly, i) = (a[i], i+1)

length{P<:Poly}(a::P) = error("define a length() for type $P")
getindex{P<:Poly}(a::P, i) = error("define a getindex() for type $P")

# delegate modulus to the coefficient type

modulus{T}(::Poly{T}) = modulus(T)
modulus{T}(::Type{Poly{T}}) = modulus(T)   # is this useful?

# generic comparison

function cmp{P<:Poly}(a::P, b::P, eq)
    if length(a) < length(b)
        true
    elseif length(a) > length(b)
        false
    else
        for i in length(a):-1:1
            if a[i] < b[i]
                return true
            elseif a[i] > b[i]
                return false
            end
        end
        eq
    end
end
<={P<:Poly}(a::P, b::P) = cmp(a, b, true)
<{P<:Poly}(a::P, b::P) = cmp(a, b, false)

# exotic polynomials may require a suitable show_coeff and/or show_mod to
# display well

show_coeff(io::IO, i::Integer) = show(io, i)
show_coeff{T,N,I<:Integer}(io::IO, z::ZModN{T,N,I}) = show(io, convert(I, z))

show_mod{I<:Integer}(io::IO, p::Poly{I}) = nothing

function print_no_mod{I<:Integer}(io::IO, p::Poly{I})
    n = length(p)
    if n <= 0
        show_coeff(io, zero(I))
    else
        for j = n:-1:1
            pj = p[j]
            if pj != zero(I)
                if j == n
                    pj < 0 && print(io, "-")
                else
                    pj < 0 ? print(io," - ") : print(io," + ")
                end
                magpj = abs(pj)
                if j == 1 || magpj != one(I)
                    show_coeff(io, magpj)
                end
                if j > 1
                    print(io, :x)
                    if j > 2
                        print(io, '^', j-1)
                    end
                end
            end
        end
    end
end

function print{I<:Integer}(io::IO, p::Poly{I})
    print_no_mod(io, p)
    show_mod(io, p)
end

function show{I<:Integer}(io::IO, p::Poly{I})
    print(io, "ZP($I")
    for i in length(p):-1:1
        print(io, ",")
        show_coeff(io, p[i])
    end
    print(io, ")")
end


# --- arbitrary (but "integer") coeff polynomials

immutable ZPoly{I<:Integer} <: Poly{I}
    a::Vector{I}
    ZPoly(a) = length(a) == 0 || a[end] != zero(I) ? new(a) : error("zero leading coeff")
end

immutable FZP{I<:Integer}
    a::Tuple
end

# constructors

function prepare_p{N}(coeffs::Vector{N})
    start = endof(coeffs)
    while start > 0 && coeffs[start] == zero(N)
        start -= 1
    end
    start == endof(coeffs) ? coeffs : coeffs[1:start]
end

ZP() = error("provide at least one (zero?) coeff")
ZP(c::Function) = error("provide at least one (zero?) coeff")
# the vector constructors are in storage order
ZP(t::Type, coeffs::Vector) = ZP(t[map(c -> convert(t, c), coeffs)...])
ZP(c::Function, coeffs::Vector) = ZP(map(c, coeffs))
ZP{I<:Integer}(coeffs::Vector{I}) = ZPoly{I}(prepare_p(coeffs))
# the direct (tuple) constructors are in natural order
ZP(t::Type, coeffs...) = ZP(reverse(t[map(c -> convert(t, c), coeffs)...]))
ZP(c::Function, coeffs...) = ZP(reverse([map(c, coeffs)...]))
ZP{I<:Integer}(coeffs::I...) = ZP(reverse([coeffs...]))

# nice syntax for creating polynomials.
# use x=X(ZF(2)) or x=X(GF2) or x=X(ZField{2,Int}) and then x^2+3 etc
X(c::Function) = ZP([c(0), c(1)])
X{I<:Integer}(::Type{I}) = ZP([zero(I), one(I)])

# this promotion to write x^2+1 etc
convert{I<:Integer}(::Type{ZPoly{I}}, p::ZPoly{I}) = p
promote_rule{I<:Integer}(::Type{ZPoly{I}}, ::Type{Int}) = ZPoly{I}
convert{I<:Integer}(::Type{ZPoly{I}}, i::Integer) = ZP([convert(I, i)])

# can be accessed and modified (but don't do that in public) as arrays
getindex(a::ZPoly, i) = getindex(a.a, i)
setindex!{T,I<:Integer}(a::ZPoly{T}, v::I, i) = setindex!(a.a, convert(T, v), i)
length(a::ZPoly) = length(a.a)

show_mod{Z<:ZModN}(io::IO, p::ZPoly{Z}) = print(io, " mod $(modulus(Z))")
showcompact(io::IO, p::ZPoly) = showcompact(io, p.a)

# number methods

zero{T}(::Type{ZPoly{T}}) = ZPoly{T}(T[])
zero{Z<:ZPoly}(z::Z) = zero(Z)
one{T}(::Type{ZPoly{T}}) = ZPoly{T}([one(T)])
one{Z<:ZPoly}(z::Z) = one(Z)

=={T}(a::ZPoly{T}, b::ZPoly{T}) = a.a == b.a
-{T}(a::ZPoly{T}) = ZPoly{T}(-a.a)

function trailing_zeros{T}(p::ZPoly{T})
    n = 1
    while n < length(p) && p[n] == zero(T)
        n += 1
    end
    n - 1
end

# a and b same length; discard leading zeroes
function _truncate{T}(f, a::Vector{T}, b::Vector{T})
    ZERO, la = zero(T), length(a)
    copied = false
    for i in la:-1:1
        x::T = f(a[i], b[i])
        if x != ZERO && !copied
            result = Array(T, i)
            copied = true
        end
        if copied
            result[i] = x
        end
    end
    if !copied
        result = T[]
    end
    result
end

# big longer than small; copy extra values
function _skip{T}(f, big::Vector{T}, small::Vector{T})
    lb, ls = length(big), length(small)
    result = Array(T, lb)
    for i in lb:-1:(ls+1)
        result[i] = big[i]
    end
    for i in ls:-1:1
        result[i] = f(big[i], small[i])
    end
    result
end

function +{T}(a::ZPoly{T}, b::ZPoly{T})
    la, lb = length(a), length(b)
    if la == 0
        b
    elseif lb == 0
        a
    elseif la == lb
        ZPoly{T}(_truncate(+, a.a, b.a))
    elseif la > lb
        ZPoly{T}(_skip(+, a.a, b.a))
    else
        ZPoly{T}(_skip(+, b.a, a.a))
    end
end

function -{T}(a::ZPoly{T}, b::ZPoly{T})
    la, lb = length(a), length(b)
    if la == lb
        if a == b
            zero(ZPoly{T})
        else
            ZPoly{T}(_truncate(-, a.a, b.a))
        end
    elseif la > lb
        ZPoly{T}(_skip(-, a.a, b.a))
    else
        # in this case, b is strictly bigger than a, and leading
        # coeffs must be negated.  so there's no chance to truncate
        # and we cannot shift past the start.  so do it here.
        result = Array(T, lb)
        for i = lb:-1:(la+1)
            result[i] = -b.a[i]
        end
        for i = la:-1:1
            result[i] = a[i] - b[i]
        end
        ZPoly{T}(result)
    end
end

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
    la, lb = length(a), length(b)
    lb > 0 || error("division by zero polynomial")
    if la == 0
        (a, a)
    elseif lb == 1 && b[1] == one(T)
        (a, T[])
    elseif lb > la
        (T[], a)
    elseif a == b
        ([one(T)], T[])
    else 
        shift = la - lb  # non-negative
        rem, div = copy(a), zeros(T, shift + 1)
        for s in 0:shift
            factor = rem[end-s] / b[end]
            if factor != zero(T)
                div[end-s] = factor
                for i in 1:lb
                    rem[i+shift-s] -= factor * b[i]
                end
            end
        end
        div, rem
    end
end

function divrem{T}(a::ZPoly{T}, b::ZPoly{T})
    map(ZP, _divrem(a.a, b.a))  # ZP discards leading zeroes
end

div{T}(a::ZPoly{T}, b::ZPoly{T}) = ZP(_divrem(a.a, b.a)[1])
/{T}(a::ZPoly{T}, b::ZPoly{T}) = ZP(_divrem(a.a, b.a)[1])
rem{T}(a::ZPoly{T}, b::ZPoly{T}) = ZP(_divrem(a.a, b.a)[2])
mod{T}(a::ZPoly{T}, b::ZPoly{T}) = ZP(_divrem(a.a, b.a)[2])
%{T}(a::ZPoly{T}, b::ZPoly{T}) = ZP(_divrem(a.a, b.a)[2])
# ^ will come from power_by_squaring in stdlib
degree{T}(a::ZPoly{T}) = length(a) - 1

# defining bitwise ops seems odd, but makes sense for GF2, especially
# when using the compact binary repn below.  for N>2 we apply the
# operation to each coefficient pair in turn.

function same_length{T}(a::Vector{T}, b::Vector{T})
    big, small = length(a) > length(b) ? (a, b) : (b, a)
    pad = length(big) - length(small)
    small = vcat(small, zeros(T, pad))
    big, small
end

map{P<:ZPoly}(f::Function, a::P, b::P) = ZP(map(f, same_length(a.a, b.a)...))
for op in (:&, :|, :$)
    @eval $op{P<:ZPoly}(a::P, b::P) = map($op, a, b)
end
<<{T}(a::ZPoly{T}, n::Int) = ZPoly{T}(vcat(zeros(T, n), a.a))
>>>{T}(a::ZPoly{T}, n::Int) = ZPoly{T}(a.a[n+1:end])
>>{T}(a::ZPoly{T}, n::Int) = a >>> n



# --- GF(2) polynomials encoded as bits

# the type of Poly{} is the type returned by []
immutable GF2Poly{U<:Unsigned} <: Poly{ZField{Int,2,Int}}
    i::U
end

# constructors
GF2P{U<:Unsigned}(p::U) = GF2Poly{U}(p)
GF2P{S<:Signed}(p::S) = GF2P(unsigned(p))
GF2X{U<:Unsigned}(::Type{U}) = GF2Poly{U}(one(U) + one(U))
GF2X() = GF2X(Uint)

modulus(::Type{GF2Poly}) = 2
modulus(::GF2Poly) = 2

convert{I<:Integer, U<:Unsigned}(::Type{GF2Poly{U}}, i::I) = GF2P(convert(U, i % 2))
promote_rule{I<:Unsigned, U<:Unsigned}(::Type{GF2Poly{U}}, ::Type{I}) = GF2Poly{U}
promote_rule{U<:Unsigned}(::Type{GF2Poly{U}}, ::Type{Int}) = GF2Poly{U}

# promote to ZPoly mainly for equality tests (convert below)
promote_rule{U<:Unsigned,J<:Integer,I<:Integer}(::Type{GF2Poly{U}}, ::Type{ZPoly{ZField{J,2,I}}}) = ZPoly{ZField{J,2,I}}

# no setindex! as immutable
getindex{U<:Unsigned}(a::GF2Poly{U}, i) = a.i & (one(U) << (i-1)) == 0 ? GF2(0) : GF2(1)
length{U<:Unsigned}(a::GF2Poly{U}) = 8 * sizeof(U) - leading_zeros(a.i)

# for display we first convert to to ZPoly.  this is not efficient,
# but we probably need the interop anyway, and who cares if printing
# is efficient?

convert{J<:Integer, I<:Integer, U<:Unsigned}(::Type{ZPoly{ZField{J,2,I}}}, p::GF2Poly{U}) = _convert(ZPoly{ZField{J,2,I}}, p)
convert{J<:Integer, I<:Integer, U<:Unsigned}(::Type{GF2Poly{U}}, p::ZPoly{ZField{J,2,I}}) = _convert(GF2Poly{U}, p)
convert{U<:Unsigned}(::Type{GF2Poly{U}}, p::GF2Poly{U}) = p

function _convert{T, F}(::Type{T}, from::F)
    result, mask = zero(T), one(T)
    for _ in 1:length(from)
        if from & one(F) != zero(F)
            result += mask
        end
        from >>>= 1
        mask <<= 1
    end
    result
end

#print{U<:Unsigned}(io::IO, p::GF2Poly{U}) = print(io, convert(ZPoly{ZField{2,Int}}, p))
#print_no_mod{U<:Unsigned}(io::IO, p::GF2Poly{U}) = print_no_mod(io, convert(ZPoly{ZField{2,Int}}, p))
show_mod{U<:Unsigned}(io::IO, p::GF2Poly{U}) = print(io, " mod 2")
show{U<:Unsigned}(io::IO, p::GF2Poly{U}) = print(io, "GF2Poly{$U}($(p.i))")

bits(p::GF2Poly) = bits(p.i)

zero{U<:Unsigned}(::Type{GF2Poly{U}}) = GF2P(zero(U))
zero{U<:Unsigned}(::GF2Poly{U}) = GF2P(zero(U))
one{U<:Unsigned}(::Type{GF2Poly{U}}) = GF2P(one(U))
one{U<:Unsigned}(::GF2Poly{U}) = GF2P(one(U))

-{U<:Unsigned}(a::GF2Poly{U}) = a
for (name, op) in ((:+, :$), (:-, :$), (:|, :|), (:&, :&), (:$, :$))
    @eval $name{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U}) = GF2Poly{U}($op(a.i, b.i))
end
for op in (:(==), :<=, :<)
    @eval $op{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U}) = $op(a.i, b.i)
end
<<{U<:Unsigned}(a::GF2Poly{U}, n::Int) = GF2Poly(a.i << n)
>>>{U<:Unsigned}(a::GF2Poly{U}, n::Int) = GF2Poly(a.i >>> n)
>>{U<:Unsigned}(a::GF2Poly{U}, n::Int) = a >>> n
leading_zeros{U<:Unsigned}(a::GF2Poly{U}) = leading_zeros(a.i)

function *{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U})
    big, small = a > b ? (a, b) : (b, a)
    if small == zero(GF2Poly{U})  # no need to test big
        small
    elseif small == one(GF2Poly{U})
        big
    else
        result = zero(GF2Poly{U})
        while small != zero(GF2Poly{U}) && big != zero(GF2Poly{U})
            if small & one(GF2Poly{U}) == one(GF2Poly{U})
                result += big
            end
            small >>>= 1
            big <<= 1
        end
        result
    end
end

function divrem{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U})
    ONE, ZERO = one(GF2Poly{U}), zero(GF2Poly{U})
    b != ZERO || error("division by zero polynomial")
    if a == ZERO
        (a, a)
    elseif b == ONE
        (a, ZERO)
    elseif a == b
        (ONE, ZERO)
    else 
        shift = leading_zeros(b) - leading_zeros(a)
        if shift < 0
            (ZERO, a)
        else 
            rem::GF2Poly{U}, div::GF2Poly{U} = a, ZERO
            factor = ONE << shift
            b = b << shift  # no overflow (b *= factor)
            mask = ONE << (8 * sizeof(U) - leading_zeros(rem) - 1)  # msb(rem)
            for _ in shift:-1:0
                if rem & mask != ZERO
                    div += factor
                    rem -= b
                end
                factor >>>= 1
                b >>>= 1
                mask >>>= 1
            end
            div, rem
        end
    end
end

div{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U}) = divrem(a, b)[1]
/{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U}) = divrem(a, b)[1]
rem{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U}) = divrem(a, b)[2]
mod{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U}) = divrem(a, b)[2]
%{U<:Unsigned}(a::GF2Poly{U}, b::GF2Poly{U}) = divrem(a, b)[2]

degree{U<:Unsigned}(p::GF2Poly{U}) = 8*sizeof(U) - leading_zeros(p)


end



#ZField{ATP ZP(IntModN.ZField{2,Int64},1,0,0,0,1,1,0,1,1),
#       IntModN.ZPoly{IntModN.ZField{2,Int64}}}(x^6 + x^5 + x^3 + x^2 + x + 1 mod 2) 
#ZField{ATP ZP(IntModN.ZField{2,Int64},1,0,0,0,1,1,0,1,1),
#       IntModN.ZPoly{IntModN.ZField{2,Int64}}}(x^6 + x^4 + x + 1 mod 2)
