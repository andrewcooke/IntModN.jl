
# equality and speed tests for ZPoly and Poly
# largely superceded by PTests2.jl

using Polynomial, IntModN

import Base: promote_rule, convert


function random_zp(N, deg)
    T = ZField{N, Int}
    ZP(rand!(T, Array(T, rand(0:deg+1))))
end

function random_p(N, deg)
    T = ZField{N, Int}
    Poly(rand!(T, Array(T, rand(0:deg+1))))
end

function random_ps(N, deg)
    T = ZField{N, Int}
    a = rand!(T, Array(T, rand(0:deg+1)))
    ZP(a), Poly(a)
end

convert{T}(::Type{ZPoly{T}}, p::Poly{T}) = ZP(p.a)
# cannot use promotion with poly as not a Number
=={T}(a::ZPoly{T}, b::Poly{T}) = a == convert(ZPoly{T}, b)
=={T}(a::Poly{T}, b::ZPoly{T}) = b == convert(ZPoly{T}, a)

function random_op()
    [+,-,/,*][rand(1:4)]
end


function test_eq()
    for i in 1:10
        p1, p2 = random_ps(5, 5)
        @assert p1 == p2
    end
    println("test_eq ok")
end

function test_op()
    for i in 1:100
        p1, p2 = random_ps(5, 5)
        q1, q2 = random_ps(5, 5)
        op = random_op()
        try
            r1 = op(p1, q1)
            r2 = op(p2, q2)
            if r1 != r2
                println("($p1) $op ($q1)")
            end
            @assert r1 == r2
        catch
            if length(q1) != 0
                println("($p1) $op ($q1)")
            end
            @assert length(q1) == 0
        end
    end
    println("test_op ok")
end


function do_rand_zp(N, deg; ops=[+,-,*,/])
    try 
        ops[rand(1:length(ops))](random_zp(N, deg), random_zp(N, deg))
    catch
        # ignore /0
    end
end

function do_rand_p(N, deg; ops=[+,-,*,/])
    try
        ops[rand(1:length(ops))](random_p(N, deg), random_p(N, deg))
    catch
        # ignore /0
    end
end

function do_timing(n, N, deg)
    # warm up
    for i in 1:n
        do_rand_zp(N, deg)
        do_rand_p(N, deg)
    end
    println("any operation")
    @time (for i in 1:n; do_rand_zp(N, deg); end)
    @time (for i in 1:n; do_rand_p(N, deg); end)
    println("division")
    @time (for i in 1:n; do_rand_zp(N, deg, ops=[/]); end)
    @time (for i in 1:n; do_rand_p(N, deg, ops=[/]); end)
    println("multiplication")
    @time (for i in 1:n; do_rand_zp(N, deg, ops=[*]); end)
    @time (for i in 1:n; do_rand_p(N, deg, ops=[*]); end)
    println("addition + subtraction")
    @time (for i in 1:n; do_rand_zp(N, deg, ops=[+,-]); end)
    @time (for i in 1:n; do_rand_p(N, deg, ops=[+,-]); end)
end


function tests()
    test_eq()
    test_op()
    do_timing(100000, 5, 10)
end
