
using Polynomial, IntModN

import Base: promote_rule, convert


function make_random(deg)
    T = ZField{2,Int}
    a = rand!(T, Array(T, rand(0:deg+1)))
    p = ZP(a)
    convert(GF2Poly{Uint}, p), p, Poly(a)
end

function make_randoms(n, deg)
    a = (GF2Poly{Uint}, ZPoly{ZField{2,Int}}, Poly{ZField{2,Int}})[]
    for _ in 1:n
        push!(a, make_random(deg))
    end
    a
end

function test_op(a, idx, op)
    for i in 1:length(a)
        p = a[i][idx]
        for j in 1:length(a)
            q = a[j][idx] 
            try
                r1 = op(p, q)
            catch
                if q != zero(q)
                    println("($p1) $op ($q1)")
                end
                @assert q == zero(q)
            end
        end
    end
end


function do_timing(n, deg)

    # warm up
    a = make_randoms(10, deg)
    for op in (+, -, *, /)
        for idx in 1:3
            test_op(a, idx, op)
        end
    end

    a = make_randoms(n, deg)
    for op in (+, -, *, /)
        println("\n$op")
        for idx in 1:3
            @time test_op(a, idx, op)
        end
    end
end


function timing()
    do_timing(1000, 8)
end
