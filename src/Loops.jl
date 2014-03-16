
# https://github.com/JuliaLang/julia/issues/6178

function repeat(op, n)
    for i in 1:n
        try
            op(0)
        catch
        end
    end
end

function fast(b)
    @assert b != 0
    for _ in 0:1
    end
end

function slow(b)
    @assert b != 0
    for _ in 1:-1:0
    end
end

for op in (fast, slow)
    repeat(op, 10)
end

for op in (fast, slow, fast, slow)
    println(op)
    @time repeat(op, 200000)
end

