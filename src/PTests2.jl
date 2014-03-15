
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
    (a, (GF2Poly{Uint}, ZPoly{ZField{2,Int}}, Poly{ZField{2,Int}}))
end

function test_op(a, idx, op, T)
    for i in 1:length(a)
        p::(T[idx]) = a[i][idx]
        for j in 1:length(a)
            q::(T[idx]) = a[j][idx] 
            try
                r1::(T[idx]) = op(p, q)
            catch
                if q != zero(q)
                    println("($p) $op ($q)")
                end
                @assert q == zero(q)
            end
        end
    end
end


function do_timing(n, deg)

    # warm up
    a, T = make_randoms(10, deg)
    for op in (+, -, *, /, %)
        for idx in 1:3
            test_op(a, idx, op, T)
        end
    end

    a, T = make_randoms(n, deg)
#    for op in (+, -, *, /, %)
    for op in (/, %)
        println("\n$op")
        for idx in 1:3
            @time test_op(a, idx, op, T)
        end
    end
end


function timing()
    do_timing(1000, 8)
end
