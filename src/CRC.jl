
# calculate crc checksums for data.  this process is closely related
# to the routines in IntModN (particularly GF2Poly), but the
# traditional implementations are optimized to the point where
# implementing them from "first principles" makes little sense
# (although we can use IntModN to check the results here).

module CRC

using IntModN

# the basic idea is to do a polynomial division along the byte stream.
# it's all explained very nicely at
# http://en.wikipedia.org/wiki/Computation_of_cyclic_redundancy_checks
# - here we implement the "byte at a time" approach, without tables.

# we assume that the generator polnomial is encoded as bits, with the
# msb discarded.  so the degree of the polynomial is 8*sizeof(G)+1
function rem_no_table{G<:Unsigned,D<:Unsigned}(generator::G, data::Vector{D})
    word_size = 8 * sizeof(D)
    carry = one(D) << (word_size - 1)
    shift = 8 * sizeof(G) - word_size
    @assert shift < 0 "polynomial smaller than data chunk"
    remainder = zero(G)
    for word in data
        remainder = remainder $ (word << shift)
        for _ in word_size
            if remainder & carry
                remainder = (remainder << 1) $ generator
            else
                remainder <<= 1
            end
        end
    end
    remainder
end



function test_rem_no_table()
    for _ in 1:100
        a = rands(Uint8, 2)
        b = rand(Uint8)
        c = convert(Uint16, a[0]) << 8 + convert(Uint16, a[1])
        x = GF2Poly(c) % GF2Poly(convert(Uint16, b))
        y = rem_no_table(b, a)
        @assert x.i == y
    end
    println("rem_no_table ok")
end

function tests()
    test_rem_no_table()
end

end
