
module Tests

using IntModN, Base.Test

export tests


# --- integer rings and fields

function test_z_constructor()

    @test string(ZR(3, 2, Int)) == "2 mod 3"
    @test string(ZR(3, 2)) == "2 mod 3"
    @test sprint(show, ZR(3, 2)) == "ZRing{Int64,3,Int64}(2)"
    @test string(ZF(3, 2, Int)) == "2 mod 3"
    @test string(ZF(3, 2)) == "2 mod 3"
    @test sprint(show, ZF(3, 2)) == "ZField{Int64,3,Int64}(2)"
    
    Z2 = ZF(2)
    @test GF2(1) == Z2(1) == ZField{Int,2,Int}(1) == ZF(2, 1)

    @test isa(ZR(4, 0x3).i, UInt8)

    try
        bad = ZR(256, 0x0)
        error("expected failure")
    catch e
        @test isa(e, ErrorException)
        @test search(string(e), "too large") != 0:-1 
    end

    @test one(ZField{Int,2,Int}) == one(GF2(0)) == GF2(1)

    println("test_z_constructor ok")
end

function test_z_random()
    @test isa(rand(ZRing{Int,3,UInt}), ZRing{Int,3,UInt})
    a = rand(ZField{Int,5,UInt8}, 2, 3)
    @test size(a) == (2, 3)
    @test isa(a, Array{ZField{Int,5,UInt8},2})
    println("test_z_random ok")
end

function test_z_arithmetic()

    @test GF2(1) + GF2(1) == GF2(0)
    @test zero(ZRing{Int,5,Int}) - one(ZRing{Int,5,Int}) == ZR(5, 4)
    @test zero(ZField{Int,5,Int}) - one(ZField{Int,5,Int}) == ZF(5, 4)
    @test inv(ZF(5, 3)) == ZF(5, 2)
    @test ZF(5, 3) * ZF(5, 2) == one(ZField{Int, 5, Int})
    @test GF2(0)^0 == GF2(1)

    for i in 1:10
        z1, z2 = ZF(5, rand(1:5)), ZF(5, rand(1:5))
        @test (z1 - z2) + z2 == z1
        @test z1 - (z2 - z2) == z1
        @test z2 + (z1 - z2) == z1
    end
    @test ZF(5, 1) - ZF(5, 4) == ZF(5, 2)

    @test modulus(GF2(1)) == order(GF2(0)) == 2

    try
        inv(ZR(6, 2))
        error("expected failure")
    catch e
        @test isa(e, ErrorException)
        @test search(string(e), "not invertible") != 0:-1 
    end

    @test ZR(4, 2) / ZR(4, 2) == ZR(4, 1)
    @test ZR(4, 2) / ZR(4, 3) == ZR(4, 2)
    @test_throws ErrorException ZR(4, 3) / ZR(4, 2)

    @test ZF(5, 2) / ZF(5, 2) == ZF(5, 1)
    @test ZF(5, 2) / ZF(5, 3) == ZF(5, 4)
    @test ZF(5, 3) / ZF(5, 2) == ZF(5, 4)

    println("test_z_arithmetic ok")
end

function test_z_matrix_inverse()

    @zfield 2 begin
        A = [1 1;
             0 1]
        b = [1, 1]
        x = A\b
        @test x == [0, 1]
    end

    # http://math.stackexchange.com/questions/169921/how-to-solve-system-of-linear-equations-of-xor-operation
    @zfield 2 begin
        A = [1 1 1 0; 
             1 1 0 1;
             1 0 1 1;
             0 1 1 1]
        b = [1, 1, 0, 1]
        x = A\b
        @test x == [0, 1, 0, 0]
    end

    println("test_z_matrix_inverse ok")
end

function test_z_power()

    # http://en.wikipedia.org/wiki/Diffie%E2%80%93Hellman_key_exchange#Explanation_including_encryption_mathematics
    base = ZF(23, 5)
    @test (base^6).i == 8
    @test (base^15).i == 19
    @test ((base^15)^6).i == 2
    @test ((base^6)^15).i == 2

    println("test_z_power ok")
end

function test_z_macros()
    begin
        b = 7
    end
    @test b == 7 
    @zfield 5  begin
        b = 1 + 2 * 3
        a = b / 3
    end
    @test a == ZF(5, 4)
    @test ZR(3, 2) == @zring 3 1 + 4
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

    x = ZP(GF2(1), GF2(0))
    @test ZP(GF2(1), GF2(0), GF2(1), GF2(0)) == x^3 + x
    x = X(ZF(5))
    @test ZP(ZF(5,2), ZF(5,3), ZF(5,4)) == 2x^2 + 3x + 4
    x = X(ZField{Int, 3, Int})
    @test x + 1 == ZP([ZF(3,1), ZF(3,1)])
    @test x + 1 == ZP(ZF(3,1), ZF(3,1))
    @test x + 1 == ZP(ZField{Int, 3, Int}, [1, 1])
    @test x + 1 == ZP(ZField{Int, 3, Int}, 1, 1)

    x = X(GF2)
    @test ZP(GF2, 1, 1, 0) == x^2 + x
    @test ZP(ZF(2), [0, 1, 1]) == x^2 + x
    @test ZP(ZField{Int,2,Int}, [0, 1, 1]) == x^2 + x

    println("test_p_constructor ok")
end

function test_p_show()

    x = X(Int)
    @test string(3x^2 + 1) == "3x^2 + 1"
    @test sprint(show, 3x^2 + 1) == "ZP(Int64,3,0,1)"

    x = X(ZF(5))
    @test string(3x^2 + 1) == "3x^2 + 1 mod 5"
    @test sprint(show, 3x^2 + 1) in 
    ("ZP(ZField{Int64,5,Int64},3,0,1)", "ZP(IntModN.ZField{Int64,5,Int64},3,0,1)")
    
    println("test_p_show ok")
end

function test_p_comparison()

    x = X(Int)
    @test x + 1 == x + 1
    @test x^2 + 1 != x + 1

    for i in 1:10
        a = ZP(rand!(Array(ZField{Int, 5, Int}, rand(0:2))))
        b = ZP(rand!(Array(ZField{Int, 5, Int}, rand(0:2))))
        @test a == b || ((a < b) $ (b < a))
        if a < b || a == b
            @test a <= b
        end
        if a > b || a == b
            @test a >= b
        end
        if a == b
            @test !(a != b)
        else
            @test a != b
        end
    end

    println("test_p_comparison ok")
end

function test_p_arithmetic()

    x = X(GF2)
    a = x^3 + x^2 + 1
    b = x^2 + 1
    p, q = divrem(a, b)
    @test string(p * b + q) == "x^3 + x^2 + 1 mod 2"

    @test divrem(x,one(x)) == (x, zero(x))
    @test divrem(x,x) == (one(x), zero(x))
    @test_throws ErrorException divrem(x,zero(x))

    x = X(ZF(5))
    p1 =  x^5 + 3x^4 + 4x^3 + 2x^2      + 1 
    p2 = 4x^5        + 2x^3 + 3x^2 +  x + 4
    ex = 2x^5 + 3x^4 + 2x^3 + 4x^2 + 4x + 2
    r = p1 - p2
    @test r == ex

    @test divrem(x,one(x)) == (x, zero(x))
    @test divrem(x,x) == (one(x), zero(x))
    @test_throws ErrorException divrem(x,zero(x))

    # http://www.math.umn.edu/~garrett/coding/Overheads/08_crcs.pdf
    # page 5
    x = X(GF2)
    @test (x^3 + x^2 + 1) + (x^3 + x + 1) == x^2 + x
    # page 6
    x = X(Int)
    @test (2x^3 + 3x^2 + x - 3) * (x^2 - 2x + 1) == (2x^5 - x^4 - 3x^3 - 2x^2 + 7x - 3)
    # page 7
    x = X(GF2)
    @test (x^3 + x + 1) * (x^2 + x + 1) == (x^5 + x^4 + 1)
    # page 8
    x = X(Int)
    @test divrem(x^8 + x^7 + x^4 + x^3 + x + 1, x^5 + x^3 + x + 1) == (x^3 + x^2 - x - 1, x^4 + 3x + 2)
    # page 9
    x = X(GF2)
    @test divrem(x^8 + x^7 + x^4 + x^3 + x + 1, x^5 + x^3 + x + 1) == (x^3 + x^2 + x + 1, x^4 + x)

    println("test_p_arithmetic ok")
end

function test_p_bitwise()
    @test ZP(3, 2, 1) & ZP(1, 1, 1) == ZP(1, 0, 1)
    @test ZP(3, 2, 1) | ZP(1, 1) == ZP(3, 3, 1)
    @test ZP(GF2(0), GF2(0)) $ ZP(GF2(1)) == ZP(GF2(0), GF2(1))    
    @test ZP(1, 0) << 2 == ZP(1, 0, 0, 0)
    @test ZP(1, 0) >>> 1 == ZP(1)

    println("test_p_bitwise")
end

function test_p_array()

    x = X(GF2)
    p = x^2 + x
    @test p[2] == GF2(1)
    @test p[1] == GF2(0)
    @test length(p) == 3

    count = 0
    for (i, a) in enumerate(p)
        count = count + 1
        @test 0 < i < 4
        @test (i == 1 && a == GF2(0)) || (i != 1 && a == GF2(1))
    end
    @test count == 3

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

    x = GF2X(UInt8)
    @test GF2P(0xa) == x^3 + x

    println("test_p2_constructor ok")
end

function test_p2_show()

    x = GF2X(UInt)
    @test string(3x^2 + 1) == "x^2 + 1 mod 2"
    @test sprint(show, 3x^2 + 1) in 
    ("GF2Poly{UInt64}(5)", "GF2Poly{UInt64}(5)")

    println("test_p2_show ok")
end

function test_p2_comparison()

    x = GF2X(UInt)
    @test x + 1 == x + 1
    @test x^2 + 1 != x + 1

    for i in 1:10
        a = GF2P(rand(UInt8))
        b = GF2P(rand(UInt8))
        @test a == b || ((a < b) $ (b < a))
        if a < b || a == b
            @test a <= b
        end
        if a > b || a == b
            @test a >= b
        end
        if a == b
            @test !(a != b)
        else
            @test a != b
        end
    end

    println("test_p2_comparison ok")
end

function test_p2_arithmetic()

    x = GF2X(UInt8)
    a = x^3 + x^2 + 1
    b = x^2 + 1
    p, q = divrem(a, b)
    @test p*b+q == a
    @test string(p * b + q) == "x^3 + x^2 + 1 mod 2"

    p1 = x^5 + x^4 + x^3 + x^2     + 1 
    p2 = x^5       + x^3 + x^2 + x
    ex =       x^4             + x + 1
    r = p1 - p2
    @test r == ex

    @test divrem(x,one(x)) == (x, zero(x))
    @test divrem(x,x) == (one(x), zero(x))
    @test_throws ErrorException divrem(x,zero(x))

    @test divrem(x^2, x) == (x, zero(x))
    @test divrem(x^2 + 1, x) == (x, one(x))
    @test divrem((x + 1) * (x^2 + 1) + 1, x + 1) == (x^2 + 1, one(x))
    @test divrem((x + 1) * (x^3 + 1) + x, x + 1) == (x^3, one(x))

    x = GF2X(UInt)  # x^8 ahead!
    # http://www.math.umn.edu/~garrett/coding/Overheads/08_crcs.pdf
    # page 5
    @test (x^3 + x^2 + 1) + (x^3 + x + 1) == x^2 + x
    # page 7
    @test (x^3 + x + 1) * (x^2 + x + 1) == (x^5 + x^4 + 1)
    # page 9
    @test divrem(x^8 + x^7 + x^4 + x^3 + x + 1, x^5 + x^3 + x + 1) == (x^3 + x^2 + x + 1, x^4 + x)

    println("test_p2_arithmetic ok")
end

function test_p2_bitwise()
    @test GF2P(0x5) & GF2P(0x7) == GF2P(0x5)
    @test GF2P(0x5) | GF2P(0x3) == GF2P(0x7)
    @test GF2P(0x0) $ GF2P(0x1) == GF2P(0x1)
    @test GF2P(0x2) << 2 == GF2P(0x8)
    @test GF2P(0x2) >>> 1 == GF2P(0x1)

    println("test_p2_bitwise")
end

function test_p2_array()

    x = GF2X()
    p = x^2 + x
    @test p[2] == GF2(1)
    @test p[1] == GF2(0)
    @test length(p) == 3

    count = 0
    for (i, a) in enumerate(p)
        count = count + 1
        @test 0 < i < 4
        @test (i == 1 && a == GF2(0)) || (i != 1 && a == GF2(1))
    end
    @test count == 3

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


function test_f_rijndael()
    x = X(GF2)
    rij = x^8 + x^4 + x^3 + x + 1
    a = ZF(rij, x^7 + x^6 + x^3 + x)
    b = ZF(rij, x^6 + x^4 + x + 1)
    o = one(a)
    @test a * b == o
    # none of these work
    @test o / a == b
    @test o / b == a
    @test inv(a) == b

    println("test_f_rijndael ok")
end

function test_f_rijndael2()
    x = GF2X(UInt64)
    rij = x^8 + x^4 + x^3 + x + 1
    a = ZF(rij, x^7 + x^6 + x^3 + x)
    b = ZF(rij, x^6 + x^4 + x + 1)
    o = one(a)
    @test a * b == o
    @test inv(a) == b
    @test o / b == a
    @test o / a == b

    println("test_f_rijndael2 ok")
end

function test_f_rijndael3()
    x = X(GF2)
    rij = x^8 + x^4 + x^3 + x + 1
    a = ZR(rij, x^7 + x^6 + x^3 + x)
    b = ZR(rij, x^6 + x^4 + x + 1)
    o = one(a)
    @test a * b == o
    # none of these work
    @test o / a == b
    @test o / b == a
    @test inv(a) == b

    println("test_f_rijndael3 ok")
end

function test_f_inverse()
    x = X(GF2)
    rij = x^8 + x^4 + x^3 + x + 1
    a = x^7 + x^6 + x^3 + x
    b = x^6 + x^4 + x + 1
    c = extended_euclidean(a, rij)
    @test b == c

    println("test_f_inverse ok")
end

function tests_f()
    test_f_inverse()
    test_f_rijndael3()
# these require correct handling of inverses
# http://www.dtic.mil/dtic/tr/fulltext/u2/a218148.pdf
#    test_f_rijndael()
#    test_f_rijndael2()
end



function tests()
    tests_z()
    tests_p()
    tests_p2()
    tests_f()
end

end
