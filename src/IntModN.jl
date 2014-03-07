
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
       rand, rand!
# Pkg.clone("https://github.com/astrieanna/TypeCheck.jl.git")
#using TypeCheck

export ZModN, ZField, ZRing, ZR, ZF, GF2, @zring, @zfield, ZPoly, ZP



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

immutable ZField{N, I<:Integer} <: ZModN{N,I}
    i::I
end

# but typically you would use one of these constructors

function validate(n, i)
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


# all the methods that make these numbers

zero{N,I}(::Type{ZRing{N,I}}) = ZRing{N,I}(zero(I))
one{N,I}(::Type{ZRing{N,I}}) = ZRing{N,I}(one(I))
zero{N,I}(::Type{ZField{N,I}}) = ZField{N,I}(zero(I))
one{N,I}(::Type{ZField{N,I}}) = ZField{N,I}(one(I))

# needed for disambiguation
convert(::Type{Bool}, z::ZRing) = convert(Bool, z.i)
convert(::Type{Bool}, z::ZField) = convert(Bool, z.i)
convert{X<:Integer}(::Type{X}, z::ZRing) = convert(X, z.i)
convert{X<:Integer}(::Type{X}, z::ZField) = convert(X, z.i)
# TODO - check these don't break other stuff
convert{N,I<:Integer}(::Type{ZRing{N,I}}, i::I) = ZR(N, i)
convert{N,I<:Integer}(::Type{ZField{N,I}}, i::I) = ZF(N, i)
modulus{N, I<:Integer}(::Type{ZModN{N, I}}) = N
modulus{T<:ZModN}(::Type{T}) = modulus(super(T))  # jameson type chain
modulus{T<:ZModN}(::T) = modulus(T)

showcompact{N,I}(io::IO, z::ZModN{N,I}) = showcompact(io, convert(I, z))
show{N,I}(io::IO, z::ZModN{N,I}) = print(io, "$(convert(I, z)) mod $(modulus(z))")

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



# --- polynomials over integers modulo n


immutable ZPoly{T<:ZModN}
    a::Array{T,1}
    ZPoly(a) = length(a) == 0 || a[end] > zero(T) ? new(a) : error("zero leading coeff")
end


# constructors

# note that index pairs (i, z) are 0-indexed, despite julia being 1-indexed
# anything else is just hopelessly counterintuitive when read as a polynomial
ZP() = error("provide at least one (zero?) coeff")
ZP(c::Function) = error("provide at least one (zero?) coeff")
function ZP{T<:ZModN}(coeffs::(Int, T)...)
    coeffs = collect(filter(c -> c[2] > zero(T), coeffs))
    m = length(coeffs) == 0 ? 0 : 1 + maximum([i for (i, z) in coeffs])
    a = zeros(T, m)
    for (i, z) in coeffs
        a[i+1] = z
    end
    ZPoly{T}(a)
end
ZP{I<:Integer}(c::Function, coeffs::(Int, I)...) = ZP([(i, c(n)) for (i, n) in coeffs]...)
ZP{T<:ZModN}(z::T...) = ZP([(i-1,z[i]) for i in 1:length(z)]...)
ZP{I<:Integer}(c::Function, i::I...) = ZP(c, [(j-1,i[j]) for j in 1:length(i)]...)

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
length(a::ZPoly) = length(a.a)

# can be accessed and modified as arrays
getindex{T<:ZModN}(a::ZPoly{T}, i) = getindex(a.a, i)
setindex!{T<:ZModN}(a::ZPoly{T}, v::T, i) = setindex!(a, v, i)
start{T<:ZModN}(a::ZPoly{T}) = 1
done{T<:ZModN}(a::ZPoly{T}, i) = i > length(a)
next{T<:ZModN}(a::ZPoly{T}, i) = (a[i], i+1)

# apply function, discarding zeros from highest index end
function apply{T}(f, a::Array{T,1}, b::Array{T,1}; shrink=true)
    na, nb = length(a), length(b)
    c = na > nb ? copy(a) : copy(b)
    shrink = shrink && na == nb
    for i in min(na, nb):-1:1
        x = f(a[i], b[i])
        if x != zero(T) && shrink
            c = c[1:i]
            shrink = false
        end
        if !shrink
            c[i] = x
        end
    end
    if shrink
        c = T[]
    end
    c
end    

liftp(c, f, a) = c(f(a.a))
liftp(c, f, a, b) = c(apply(f, a.a, b.a))
-{T<:ZModN}(a::ZPoly{T}) = liftp(ZPoly{T}, -, a)
+{T<:ZModN}(a::ZPoly{T}, b::ZPoly{T}) = liftp(ZPoly{T}, +, a, b)
-{T<:ZModN}(a::ZPoly{T}, b::ZPoly{T}) = liftp(ZPoly{T}, -, a, b)

function *{T<:ZModN}(a::ZPoly{T}, b::ZPoly{T})
    (long, short) = length(a) > length(b) ? (a, b) : (b, a)
    if length(short) == 0
        short
    elseif length(short) == 1 && short[1] == one(T)
        long
    else
        result = zeros(T, length(long) + length(short) - 1)
        for (i, z) in enumerate(short)
            if z != zero(T)
                for j in 1:length(long)
                    result[i+j-1] += long[j] * z
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
        result = zeros(T, length(long) + length(short) - 1)
        for (i, z) in enumerate(short)
            if z != zero(T)
                for j in 1:length(long)
                    result[i+j-1] += long[j] * z
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

end
