
module IntModN

import Base: show

type Z{N, I<:Integer}
    int::I
    function Z(n)
        @assert isa(N, Int) "N ($N) not an Int"
        @assert n > 0 "n ($n) too small"
        @assert n < N "n ($n) too large (>= $N)"
        new(n)
    end
end
Z{I<:Integer}(n::I, N::Int) = Z{N, I}(n)

typealias GF2 Z{2, Int}

function show{N,I}(io::IO, z::Z{N,I})
    # TODO - something more compact?
    print(io, "$(z.int) mod $N")
end

=={N,I}(a::Z{N,I}, b::Z{N,I}) = a.int == b.int



function test_constructor()

    @assert string(Z{3,Int}(2)) == "2 mod 3"
    @assert string(Z(0x3, 4)) == "3 mod 4"
    @assert string(GF2(1)) == "1 mod 2"
    @assert GF2(1) == Z{2,Int}(1) == Z(1, 2)

    @assert isa(Z(0x3, 4).int, Uint8)

    try
        bad = Z(2, 1)
        error("expected failure")
    catch e
        @assert isa(e, ErrorException)
        @assert search(string(e), "too large") != 0:-1 
    end

    println("test_contructor ok")
end

function tests()
    test_constructor()
end

#tests()

end

