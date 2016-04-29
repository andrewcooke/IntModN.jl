

# equality and speed tests for ZPoly, GF2Poly, and Poly

module PTests

using Polynomials, IntModN, Compat, Base.Test

import Base: promote_rule, convert, ==

export tests


convert{T}(::Type{ZPoly{T}}, p::Poly{T}) = ZP(reverse(p.a))
# cannot use promotion with poly as not a Number
=={T}(a::ZPoly{T}, b::Poly{T}) = a == convert(ZPoly{T}, b)
=={T}(a::Poly{T}, b::ZPoly{T}) = convert(ZPoly{T}, a) == b
=={U<:Unsigned,J<:Integer,I<:Integer}(a::GF2Poly{U}, b::Poly{ZField{J,2,I}}) = convert(ZPoly{ZField{J,2,I}}, a) == convert(ZPoly{ZField{J,2,I}}, b)
=={U<:Unsigned,J<:Integer,I<:Integer}(a::Poly{ZField{J,2,I}}, b::GF2Poly{U}) = convert(ZPoly{ZField{J,2,I}}, a) == convert(ZPoly{ZField{J,2,I}}, b)



function make_random(deg, modulus)
    T = ZField{Int,modulus,Int}
    a = rand!(T, Array(T, rand(0:deg+1)))
    p = ZP(a)
    p, Poly(reverse(a))
end

function make_randoms(n, deg, modulus)
    a = @compat Tuple{ZPoly{ZField{Int,modulus,Int}}, Poly{ZField{Int,modulus,Int}}}[]
    for _ in 1:n
        push!(a, make_random(deg, modulus))
    end
    (a, @compat Tuple{ZPoly{ZField{Int,modulus,Int}}, Poly{ZField{Int,modulus,Int}}})
end

function make_random2(deg)
    T = ZField{Int,2,Int}
    a = rand!(T, Array(T, rand(0:deg+1)))
    p = ZP(a)
    convert(GF2Poly{UInt}, p), p, Poly(reverse(a))
end

function make_randoms2(n, deg)
    a = @compat Tuple{GF2Poly{UInt}, ZPoly{ZField{Int,2,Int}}, Poly{ZField{Int,2,Int}}}[]
    for _ in 1:n
        push!(a, make_random2(deg))
    end
    (a, (GF2Poly{UInt}, ZPoly{ZField{Int,2,Int}}, Poly{ZField{Int,2,Int}}))
end

function test_op(a, idx, op, T)
    check_zero = in(op, (rem, div, %, /))
    ZERO = zero(T[1])
    for i in 1:length(a)
        for j in 1:length(a)
            if !check_zero || a[j][1] != ZERO
                op(a[i][idx], a[j][idx])
            end
        end
    end
end

function do_timing(n, deg)

    # warm up
    a, T = make_randoms2(10, deg)
    for op in (+, -, *, /, %)
        for idx in 1:3
            test_op(a, idx, op, T)
        end
    end

    srand(1)
    a, T = make_randoms2(n, deg)
    println("\nfirst poly: ", a[1][1])
    for idx in 1:3
        print("$(idx): $(T[idx])]\n")
    end
    for op in (+, -, *, /, %)
        println("\n$op")
        for idx in 1:3
            gc()
            if VERSION < v"0.4-"
                gc_disable()
            else
                gc_enable(false)
            end
            print("$(idx): ")
            @time test_op(a, idx, op, T)
            if VERSION < v"0.4-"
                gc_enable()
            else
                gc_enable(true)
            end
        end
    end
end

#do_timing(1000, 8)


function test_eq(a, T, op)
    for i1 in 1:length(a)
        for i2 in 1:length(a)
            for t1 in 1:length(a[1])  # all same length
                p1 = a[i1][t1]
                q1 = a[i2][t1]
                for x in op
                    r1 = nothing
                    try
                        r1 = x(p1, q1)
                    catch
                    end
                    for t2 in 1:length(a[1])
#                        println("test $i1 $i2 $(T[t1]) $(T[t2]) $x")
                        p2 = a[i1][t2]
                        q2 = a[i2][t2]
                        r2 = nothing
                        try
                            r2 = x(p2, q2)
                        catch
                        end
                        if r1 != r2
                            println("($p1) $x ($q1)")
                            println("r1  $(typeof(r1))  $(r1)")
                            println("($p2) $x ($q2)")
                            println("r2  $(typeof(r2))  $(r2)")
                        end
                        @test r1 == r2
                    end
                end
            end
        end
    end
end

function do_eq(n, deg)
    print("\ntesting eq (mod 2)...")
    a, T = make_randoms2(n, deg)
    test_eq(a, T, (+, -, *, /, %))
    println("ok")

    print("testing eq (mod 5)...")
    a, T = make_randoms(n, deg, 5)
    test_eq(a, T, (+, -, *, /, %))
    println("ok")
end

#do_eq(30, 8)


function ptests()
    do_timing(1000, 8)
    do_eq(30, 2)
end

end
