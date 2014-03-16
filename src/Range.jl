
# incomplete experiment for faster range using type to store step

immutable MyRange{N} <: Ranges{Int}
    start::Int
    len::Int
end

function mycolon(start::Int, step::Int, stop::Int)
    MyRange{step}(start, (stop - start) // step)
end

next{N}(r::MyRange{N}, i) = i + N::Int

function slow1(n)
    for _ in n:-1:1
    end
end

function fast1(n)
    for _ in mycolon(n, -1, 1)
    end
end

code_native(slow1, (Int,))
code_native(fast1, (Int,))

