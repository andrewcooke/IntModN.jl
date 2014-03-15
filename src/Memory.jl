
immutable Fast{U<:Unsigned}
    i::U
end
+{U<:Unsigned}(a::Fast{U}, b::Fast{U}) = Fast{U}(a.i $ b.i)

immutable Slow{I<:Unsigned}
    i::I
end
for (name, op) in ((:+, :$),)
    @eval $name{U<:Unsigned}(a::Slow{U}, b::Slow{U}) = Slow{U}(($op)(a.i, b.i))
end
