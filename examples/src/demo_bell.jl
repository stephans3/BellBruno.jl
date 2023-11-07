using BellBruno
using LinearAlgebra

N = 10;
T = diagm(ones(Int64,N))
for n=2:N
    for k=1:n-1
        T[n,k] = sum(T[n-k,1:k]);
    end
end


bp = bell_poly(N)
for (i, el) in enumerate(bp)
    length(el) == i-1
end

for (i, e1) in enumerate(bp)
    if i > 1
        for (k,e2) in enumerate(e1)
            if typeof(e2) == Vector{Int8}
                display(1 == T[i-1,k]) 
            elseif typeof(e2) == Matrix{Int8}
                display(size(e2)[2] == T[i-1,k]) 
            end
        end
    end
end

(typeof(bp[4][1]) == Vector{Int8} ? 1 : size(bp[4][1])) == T[3][1]

bp[2][1] == [1]
bp[3]
bp[4]
