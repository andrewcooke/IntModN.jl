
function slow1(n)
    for _ in n:-1:1
    end
end

function fast1(n)
    for _ in 1:n
    end
end

code_native(slow1, (Int,))
code_native(fast1, (Int,))
