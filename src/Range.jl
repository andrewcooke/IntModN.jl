
# this tests the idea that if we push the iterator into the type
# system we'll get fast code.  we don't.  :o(

module RangeX
import Base: first, step, last, next, done, length, isempty, getindex, start

# this is used to select "out" colon method and so trigger the new
# behaviour
immutable MyInt <: Integer
    i::Int
end

immutable TypedRange{S,I,E} <: Ranges{Int}
end

function colon(start::MyInt, inc::MyInt, end_::MyInt)
    TypedRange{start.i, inc.i, end_.i}()
end

start{S,I,E}(r::TypedRange{S,I,E}) = S - I
next{S,I,E}(r::TypedRange{S,I,E}, i::Int) = (i + I, i + I)
#done{S,I,E}(r::TypedRange{S,I,E}, i::Int) = I > 0 ? i >= E : i <= E
# even if we omit the check (so only support count down, it's slow)
done{S,I,E}(r::TypedRange{S,I,E}, i::Int) = i <= E

# these don't seem to be needed for this simple demo
#first{S,I,E}(r::TypedRange{S,I,E}) = S
#step{S,I,E}(r::TypedRange{S,I,E}) = I
#last{S,I,E}(r::TypedRange{S,I,E}) = E
#length{S,I,E}(r::TypedRange{S,I,E}) = 1 + div(E - S, I)
#isempty{S,I,E}(r::TypedRange{S,I,E}) = (I > 0 && S > E) || (I < 0 && S < E)
#getindex{S,I,E}(::TypedRange{S,I,E}, i::Int) = S + I * (i - 1)

function count_type(r, k)
    n = 0
    for i in r
        n += i % 2 == 1 ? 2k : -k
    end
    return n
end

n = 1 << 20

println(count_type(n:-1:1, 3))
println(typeof(n:-1:1))

@time count_type(n:-1:1, 3)
@time count_type(MyInt(n):MyInt(-1):MyInt(1), 3)

#code_native(count_type, (Range{Int}, Int))
#code_native(count_type, (TypedRange{n,-1,1}, Int))
code_llvm(count_type, (Range{Int}, Int))
code_llvm(count_type, (TypedRange{n,-1,1}, Int))

end
