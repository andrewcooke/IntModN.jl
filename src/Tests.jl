
module Tests

using IntModN

export tests


# --- integer rings and fields

function test_z_constructor()

    @assert string(ZR(3, 2, Int)) == "2 mod 3"
    @assert string(ZR(3, 2)) == "2 mod 3"
    @assert sprint(show, ZR(3, 2)) == "ZRing{3,Int64}(2)"
    @assert string(ZF(3, 2, Int)) == "2 mod 3"
    @assert string(ZF(3, 2)) == "2 mod 3"
    @assert sprint(show, ZF(3, 2)) == "ZField{3,Int64}(2)"
    
    Z2 = ZF(2)
    @assert GF2(1) == Z2(1) == ZField{2,Int}(1) == ZF(2, 1)

    @assert isa(ZR(4, 0x3).i, Uint8)

    try
        bad = ZR(256, 0x0)
        error("expected failure")
    catch e
        @assert isa(e, ErrorException)
        @assert search(string(e), "too large") != 0:-1 
    end

    @assert one(ZField{2,Int}) == one(GF2(0)) == GF2(1)

    println("test_z_constructor ok")
end

function test_z_random()
    @assert isa(rand(ZRing{3,Uint}), ZRing{3,Uint})
    a = rand(ZField{5,Uint8}, 2, 3)
    @assert size(a) == (2, 3)
    @assert isa(a, Array{ZField{5,Uint8},2})
    println("test_z_random ok")
end

function test_z_arithmetic()

    @assert GF2(1) + GF2(1) == GF2(0)
    @assert zero(ZRing{5,Int}) - one(ZRing{5,Int}) == ZR(5, 4)
    @assert zero(ZField{5,Int}) - one(ZField{5,Int}) == ZF(5, 4)
    @assert inv(ZF(5, 3)) == ZF(5, 2)
    @assert ZF(5, 3) * ZF(5, 2) == one(ZField{5, Int})
    @assert GF2(0)^0 == GF2(1)

    for i in 1:10
        z1, z2 = ZF(5, rand(1:5)), ZF(5, rand(1:5))
        @assert (z1 - z2) + z2 == z1
        @assert z1 - (z2 - z2) == z1
        @assert z2 + (z1 - z2) == z1
    end
    @assert ZF(5, 1) - ZF(5, 4) == ZF(5, 2)

    @assert modulus(GF2(1)) == order(GF2(0)) == 2

    try
        inv(ZR(6, 2))
        error("expected failure")
    catch e
        @assert isa(e, ErrorException)
        @assert search(string(e), "not invertible") != 0:-1 
    end

    println("test_z_arithmetic ok")
end

function test_z_matrix_inverse()

    @zfield 2 begin
        A = [1 1;
             0 1]
        b = [1, 1]
        x = A\b
        @assert x == [0, 1]
    end

    # http://math.stackexchange.com/questions/169921/how-to-solve-system-of-linear-equations-of-xor-operation
    @zfield 2 begin
        A = [1 1 1 0; 
             1 1 0 1;
             1 0 1 1;
             0 1 1 1]
        b = [1, 1, 0, 1]
        x = A\b
        @assert x == [0, 1, 0, 0]
    end

    println("test_z_matrix_inverse ok")
end

function test_z_power()

    # http://en.wikipedia.org/wiki/Diffie%E2%80%93Hellman_key_exchange#Explanation_including_encryption_mathematics
    base = ZF(23, 5)
    @assert (base^6).i == 8
    @assert (base^15).i == 19
    @assert ((base^15)^6).i == 2
    @assert ((base^6)^15).i == 2

    println("test_z_power ok")
end

function test_z_macros()
    begin
        b = 7
    end
    @assert b == 7 
    @zfield 5  begin
        b = 1 + 2 * 3
        a = b / 3
    end
    @assert a == ZF(5, 4)
    @assert ZR(3, 2) == @zring 3 1 + 4
    println("test_z_macros ok")
end

function test_z_coverage()
    println("Integer", methodswithdescendants(Integer, lim=20))
    println("ZRing", methodswithdescendants(ZRing, lim=20))
    println("test_z_coverage ok")
end

function tests_z()
    test_z_constructor()
    test_z_random()
    test_z_arithmetic()
    test_z_matrix_inverse()
    test_z_power()
    test_z_macros()
#    test_z_coverage()
end


# --- polynomials


function test_p_constructor()

    x = ZP([GF2(1), GF2(0)])
    @assert ZP([GF2(1), GF2(0), GF2(1), GF2(0)]) == x^3 + x
    x = X(ZF(5))
    @assert ZP([ZF(5,2), ZF(5,3), ZF(5,4)]) == 2x^2 + 3x + 4
    x = X(ZField{3, Int})
    @assert x + 1 == ZP([ZF(3,1), ZF(3,1)])
    @assert x + 1 == ZP(ZF(3,1), ZF(3,1))
    @assert x + 1 == ZP(ZField{3, Int}, [1, 1])
    @assert x + 1 == ZP(ZField{3, Int}, 1, 1)

    x = X(GF2)
    @assert ZP(GF2, 1, 1, 0) == x^2 + x
    @assert ZP(ZF(2), [1, 1, 0]) == x^2 + x
    @assert ZP(ZField{2,Int}, [1, 1, 0]) == x^2 + x

    p = x^2 + x
    @assert IntModN.encode_factor(p) == (1, 1, 0)
    @assert IntModN.decode_factor(ZPoly{ZField{2,Int}}, (1, 1, 0)) == p

    println("test_p_constructor ok")
end

function test_p_show()

    x = X(Int)
    @assert string(3x^2 + 1) == "3x^2 + 1"
    @assert sprint(show, 3x^2 + 1) == "ZP(Int64,3,0,1)"

    x = X(ZF(5))
    @assert string(3x^2 + 1) == "3x^2 + 1 mod 5"
    @assert sprint(show, 3x^2 + 1) == "ZP(ZField{5,Int64},3,0,1)"
    
    println("test_p_show ok")
end

function test_p_comparison()

    x = X(Int)
    @assert x + 1 == x + 1
    @assert x^2 + 1 != x + 1

    for i in 1:10
        a = ZP(rand!(Array(ZField{5, Int}, rand(0:2))))
        b = ZP(rand!(Array(ZField{5, Int}, rand(0:2))))
        @assert a == b || ((a < b) $ (b < a))
        if a < b || a == b
            @assert a <= b
        end
        if a > b || a == b
            @assert a >= b
        end
        if a == b
            @assert !(a != b)
        else
            @assert a != b
        end
    end

    println("test_p_comparison ok")
end

function test_p_arithmetic()

    x = X(GF2)
    a = x^3 + x^2 + 1
    b = x^2 + 1
    p, q = divrem(a, b)
    @assert string(p * b + q) == "x^3 + x^2 + 1 mod 2"

    x = X(ZF(5))
    p1 =  x^5 + 3x^4 + 4x^3 + 2x^2      + 1 
    p2 = 4x^5        + 2x^3 + 3x^2 +  x + 4
    ex = 2x^5 + 3x^4 + 2x^3 + 4x^2 + 4x + 2
    r = p1 - p2
    @assert r == ex

    println("test_p_arithmetic ok")
end

function test_p_bitwise()
    @assert ZP(3, 2, 1) & ZP(1, 1, 1) == ZP(1, 0, 1)
    @assert ZP(3, 2, 1) | ZP(1, 1) == ZP(3, 3, 1)

    @assert ZP(GF2(0), GF2(0)) $ ZP(GF2(1)) == ZP(GF2(0), GF2(1))    

    println("test_p_bitwise")
end

function test_p_array()

    x = X(GF2)
    p = x^2 + 1
    p[2] = 1
    @assert p == x^2 + x + 1
    @assert p[1] == GF2(1)
    @assert length(p) == 3

    count = 0
    for (i, a) in enumerate(p)
        count = count + 1
        @assert 0 < i < 4
        @assert a == GF2(1)
    end
    @assert count == 3

    println("test_p_array")
end

function tests_p()
    test_p_constructor()
    test_p_show()
    test_p_comparison()
    test_p_arithmetic()
    test_p_bitwise()
    test_p_array()
end


# --- fast polynomials in GF2


function test_p2_constructor()

    x = GF2X(Uint8)
    @assert GF2P(0xa) == x^3 + x

    p = x^2 + x
    @assert IntModN.encode_factor(p) == tuple(6)
    @assert IntModN.decode_factor(GF2Poly{Uint8}, tuple(6)) == p

    println("test_p2_constructor ok")
end

function test_p2_show()

    x = GF2X(Uint)
    @assert string(3x^2 + 1) == "x^2 + 1 mod 2"
    @assert sprint(show, 3x^2 + 1) == "GF2Poly{Uint64}(5)"

    println("test_p2_show ok")
end

function test_p2_comparison()

    x = GF2X(Uint)
    @assert x + 1 == x + 1
    @assert x^2 + 1 != x + 1

    for i in 1:10
        a = GF2P(rand(Uint8))
        b = GF2P(rand(Uint8))
        @assert a == b || ((a < b) $ (b < a))
        if a < b || a == b
            @assert a <= b
        end
        if a > b || a == b
            @assert a >= b
        end
        if a == b
            @assert !(a != b)
        else
            @assert a != b
        end
    end

    println("test_p2_comparison ok")
end

function test_p2_arithmetic()

    x = GF2X(Uint8)
    a = x^3 + x^2 + 1
    b = x^2 + 1
    p, q = divrem(a, b)
    @assert string(p * b + q) == "x^3 + x^2 + 1 mod 2"

    p1 = x^5 + x^4 + x^3 + x^2     + 1 
    p2 = x^5       + x^3 + x^2 + x
    ex =       x^4             + x + 1
    r = p1 - p2
    @assert r == ex

    println("test_p2_arithmetic ok")
end

function test_p2_bitwise()
    @assert GF2P(0x5) & GF2P(0x7) == GF2P(0x5)
    @assert GF2P(0x5) | GF2P(0x3) == GF2P(0x7)
    @assert GF2P(0x0) $ GF2P(0x1) == GF2P(0x1)

    println("test_p2_bitwise")
end

function test_p2_array()

    x = GF2X()
    p = x^2 + 1
    @assert p[1] == GF2(1)
    @assert length(p) == 3

    count = 0
    for (i, a) in enumerate(p)
        count = count + 1
        @assert 0 < i < 4
        @assert (i == 2 && a == GF2(0)) || (i != 2 && a == GF2(1))
    end
    @assert count == 3

    println("test_p2_array")
end

function tests_p2()
    test_p2_constructor()
    test_p2_show()
    test_p2_comparison()
    test_p2_arithmetic()
    test_p2_bitwise()
    test_p2_array()
end


# --- factor rings of polynomials


function test_f_constructor()
    x = X(GF2)
    f = FR(x^3 + x^2, x^2 + 1)
    @assert string(f) == "x + 1 mod x^2 + 1 mod 2"
    @assert sprint(show, f) == "FR(ZP(ZField{2,Int64},1,1),ZP(ZField{2,Int64},1,0,1))"

    println("test_f_constructor ok")
end

function test_f_rijndael()
    x = X(GF2)
    rij = x^8 + x^4 + x^3 + x + 1
    a = FR(x^7 + x^6 + x^3 + x, rij)
    b = FR(x^6 + x^4 + x + 1, rij)
    o = one(a)
    @assert a * b == o
    @assert inv(a) == b
    @assert o / b == a
    @assert o / a == b

    println("test_f_rijndael ok")
end

function test_f_inverse()
    x = X(GF2)
    rij = x^8 + x^4 + x^3 + x + 1
    a = x^7 + x^6 + x^3 + x
    b = x^6 + x^4 + x + 1
    c = extended_euclidean(a, rij)
    @assert b == c

    println("test_f_inverse ok")
end

function tests_f()
    test_f_constructor()
    test_f_rijndael()
    test_f_inverse()
end



function tests()
    tests_z()
    tests_p()
    tests_p2()
    tests_f()
end

end
