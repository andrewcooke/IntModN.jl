
# this tests the idea that if we push the iterator into the type
# system we'll get fast code.  we don't.  :o(

module Range
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
done{S,I,E}(r::TypedRange{S,I,E}, i::Int) = I > 0 ? i >= E : i <= E
# even if we omit the check (so only support count down, it's slow)
#done{S,I,E}(r::TypedRange{S,I,E}, i::Int) = i <= E

# these don't seem to be needed for this simple demo
#first{S,I,E}(r::TypedRange{S,I,E}) = S
#step{S,I,E}(r::TypedRange{S,I,E}) = I
#last{S,I,E}(r::TypedRange{S,I,E}) = E
#length{S,I,E}(r::TypedRange{S,I,E}) = 1 + div(E - S, I)
#isempty{S,I,E}(r::TypedRange{S,I,E}) = (I > 0 && S > E) || (I < 0 && S < E)
#getindex{S,I,E}(::TypedRange{S,I,E}, i::Int) = S + I * (i - 1)

function count(s, i, e; display=true)
    for i in s:i:e
        if display
            println(i)
        end
    end
end

# demo
count(MyInt(10), MyInt(-1), MyInt(1))
count(10, -1, 1)

# warm up the jit
count(MyInt(10), MyInt(-1), MyInt(1); display=false)
count(10, -1, 1; display=false)

n = 1 << 20
@time count(MyInt(n), MyInt(-1), MyInt(1); display=false)
@time count(n, -1, 1; display=false)

# typical output
# elapsed time: 0.104878996 seconds (67153608 bytes allocated)
# elapsed time: 0.000323335 seconds (128 bytes allocated)

# so this idea results in a 300x slowdown.

end
