
module Tests

using IntModN

export tests


# --- integer rings and fields

function test_z_constructor()

    @assert string(ZR(3, 2, Int)) == "2 mod 3"
    @assert string(ZR(3, 2)) == "2 mod 3"
    @assert string(ZF(3, 2, Int)) == "2 mod 3"
    @assert string(ZF(3, 2)) == "2 mod 3"
    
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
    println("test_z_constructor ok")
end

function test_z_random()
    @assert isa(rand(ZRing{3,Uint}), ZRing{3,Uint})
    a = rand(ZField{4,Uint8}, 2, 3)
    @assert size(a) == (2, 3)
    @assert isa(a, Array{ZField{4,Uint8},2})
    println("test_z_random ok")
end

function test_z_arithmetic()

    @assert GF2(1) + GF2(1) == GF2(0)
    @assert zero(ZRing{5,Int}) - one(ZRing{5,Int}) == ZR(5, 4)
    @assert zero(ZField{5,Int}) - one(ZField{5,Int}) == ZF(5, 4)
    @assert inv(ZF(5, 3)) == ZF(5, 2)
    @assert ZF(5, 3) * ZF(5, 2) == one(ZField{5, Int})
    @assert GF2(0)^0 == GF2(1)

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
    @assert ZPoly{ZField{2,Int}}([GF2(0), GF2(1)]) == ZP(ZF(2), (1, 1))
    @assert ZP(GF2, 0, 0, 0) == zero(ZPoly{ZField{2,Int}})
    @zfield 3 begin
        @assert ZP(1,2,3,4) == ZP(1,2,0,1)
    end
    println("test_p_constructor ok")
end

function test_p_type()
    @assert convert(Array{ZField{2,Int},1}, ZP(GF2, 0, 1)) == [GF2(0), GF2(1)]
    @assert string(ZP(GF2, 1, 0, 1)) == "x^2 + 1 mod 2"
    println("test_p_type ok")
end

function test_p_arithmetic()
    @assert ZP(GF2, (3,1), (0,1)) + ZP(GF2, (2,1), (0,1)) == ZP(GF2, 0, 0, 1, 1)
    # http://en.wikipedia.org/wiki/Finite_field_arithmetic#Rijndael.27s_finite_field
    a = ZP(GF2, (6,1), (4,1), (1,1), (0,1))
    b = ZP(GF2, (7,1), (6,1), (3,1), (1,1))
    @assert a * b == ZP(GF2, (13,1), (12,1), (11,1), (10,1), (9,1), (8,1), (6,1), (5,1), (4,1), (3,1), (2,1), (1,1))
    println("test_p_arithmetic ok")
end

function tests_p()
    test_p_constructor()
    test_p_type()
    test_p_arithmetic()
end








function tests()
    tests_z()
    tests_p()
end

end
