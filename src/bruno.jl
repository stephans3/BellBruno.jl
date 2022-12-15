# c1 + c2*t + c3*t^2 + c4*t^3 + ...
function power_series(t, c)
    N = length(c)
    tvec = t*ones(N)
    tvec = map(^, tvec, collect(0:N-1))
    return mapreduce(*,+,tvec,c)
end

function power_series_der(t, c, n :: Int64)
    N = length(c)
    c_new = zeros(BigFloat, N-n)

    for i in 1:(N-n)
        c_new[i] = factorial(big(i-1+n))/factorial(big(i-1)) * c[i+n]
    end

    return power_series(t, c_new)
end


# c*t^p
function simple_monomial(t, c, p)
    return c*t^p
end

# t^(p-n) * c * π(p-j) for j=0 to n-1
function simple_monomial_der(t, c, p, n :: Int64)
    arg_prod = BigFloat.(p .- collect(0:n-1))
    c_new = c*reduce(*, arg_prod)
    return c_new * t^(p-n)
end
