
module IntModN

import Base: show, zero, one, inv, real, abs, convert, promote_rule

# much thanks to Andreas Noack Jensen
# https://groups.google.com/d/msg/julia-users/Ui977brdzAU/u4rWiDeJv-MJ
# https://gist.github.com/andreasnoackjensen/9132511
# (but since this changes things significanty, mistakes are mine...)

# the underlying integer type can be specified, but typically (via the
# Z function below) will be Int.
type Z{N, I<:Integer}<:Number
    int::I
    function Z(n)
        @assert isa(N, Int) "N ($N) not an Int"
        @assert N > 0 "N ($N) too small"
        @assert N <= typemax(I) "N ($N) too large for $I"
        @assert n >= 0 "n ($n) too small"
        @assert n < N "n ($n) too large (>= $N)"
        new(n)
    end
end

# infer the inderlying type and reduce to canonical form
Z{I<:Integer}(n::I, N::Int) = Z{N, I}(convert(I, mod(n, N)))

# TODO - more aliases?
typealias GF{N} Z{N, Int}
typealias GF2 GF{2}

# TODO - something compact that includes N?
show{N,I}(io::IO, z::Z{N,I}) = print(io, string(z.int))

# there is no conversion between different parameterisations of Z and
# equality is strictly for matching types.
=={N,I}(a::Z{N,I}, b::Z{N,I}) = a.int == b.int
real{N,I}(a::Z{N,I}) = real(a.int)
abs{N,I}(a::Z{N,I}) = a
<{N,I}(a::Z{N,I}, b::Z{N,I}) = a.int < b.int

zero{N,I}(::Type{Z{N,I}}) = Z{N,I}(zero(I))
one{N,I}(::Type{Z{N,I}}) = Z{N,I}(one(I))

# we do allow promotion to Int/Uint (necessary for linalg to work)
promote_rule{N,I<:Unsigned}(::Type{Z{N,I}},::Type{Uint}) = Uint
promote_rule{N,I<:Integer}(::Type{Z{N,I}},::Type{Int}) = Int
convert{N,I<:Unsigned}(::Type{Uint}, a::Z{N,I}) = convert(Uint, a.int)
convert{N,I<:Integer}(::Type{Int}, a::Z{N,I}) = convert(Int, a.int)

+{N,I}(a::Z{N,I}, b::Z{N,I}) = Z{N,I}(convert(I, mod(a.int + b.int, N)))
-{N,I}(a::Z{N,I}, b::Z{N,I}) = Z{N,I}(convert(I, mod(a.int - b.int, N)))
*{N,I}(a::Z{N,I}, b::Z{N,I}) = Z{N,I}(convert(I, mod(a.int * b.int, N)))

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

# TODO - cache values
inv{N,I}(a::Z{N,I}) = Z{N,I}(inverse(a.int, convert(I, N)))
/{N,I}(a::Z{N,I}, b::Z{N,I}) = a * inv(b)


function test_constructor()

    @assert string(Z{3,Int}(2)) == "2"
    @assert string(Z(0x3, 4)) == "3"
    @assert string(GF2(1)) == "1"
    @assert GF2(1) == GF{2}(1) == Z{2,Int}(1) == Z(1, 2)

    @assert isa(Z(0x3, 4).int, Uint8)

    try
        bad = Z{1, Int}(2)
        error("expected failure")
    catch e
        @assert isa(e, ErrorException)
        @assert search(string(e), "too large") != 0:-1 
    end

    println("test_constructor ok")
end

function test_arithmetic()

    @assert GF2(1) + GF2(1) == GF2(0)
    @assert zero(Z{5,Int}) - one(Z{5,Int}) == Z(4, 5)
    @assert inv(GF{5}(3)) == GF{5}(2)
    @assert GF{5}(3) * GF{5}(2) == one(GF{5})

    try
        inv(GF{6}(2))
        error("expected failure")
    catch e
        @assert isa(e, ErrorException)
        @assert search(string(e), "not invertible") != 0:-1 
    end

    println("test_arithmetic ok")
end

function test_matrix_inverse()

    A = [one(GF2) one(GF2);
         zero(GF2) one(GF2)]
    b = [one(GF2), one(GF2)]
    x = A\b
    @assert x == [zero(GF2), one(GF2)]

    println("test_matrix_inverse ok")
end

function tests()
    test_constructor()
    test_arithmetic()
    test_matrix_inverse()
end

tests()

end

