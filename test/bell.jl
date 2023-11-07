
@testset "Length" begin
    N = 10;
    bp = bell_poly(N)
    for (i, el) in enumerate(bp)
        @test length(el) == i-1
    end

    T = diagm(ones(Int64,N))
    for n=2:N
        for k=1:n-1
            T[n,k] = sum(T[n-k,1:k]);
        end
    end

    for (i, e1) in enumerate(bp)
        if i > 1
            for (k,e2) in enumerate(e1)
                if typeof(e2) == Vector{Int8}
                    @test 1 == T[i-1,k]
                elseif typeof(e2) == Matrix{Int8}
                    @test size(e2)[2] == T[i-1,k] 
                end
            end
        end
    end
end

