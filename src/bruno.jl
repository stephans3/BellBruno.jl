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


"""
    faa_di_bruno(f_der, g, g_der, tgrid,bp,bc)

Faà di Bruno's algorithm to find derivatives of function composition `f(g(t))`

# Arguments

`f_der`: derivative of outer function

`g`: inner function

`g_der`: derivative of inner function

`tgrid`: sampling points (e.g. time grid)

`bp`: Bell polynomial data

`bc`: Bell coefficients data
"""
function faa_di_bruno(f_der, g, g_der, tgrid,bp,bc)
    nt = length(tgrid)
    N = length(bp)-1
    ngrid = 1:N
    g_data = g.(tgrid)
    f_der_data = f_der.(g_data,ngrid')
    g_der_data = g_der.(tgrid,ngrid')

    q = zeros(nt, N);
    for ν=1 : N
        for k=1 : ν
            fi = firstindex(bp[ν+1][k][1,:])
            li = lastindex(bp[ν+1][k][1,:])
            sol_prod = zeros(BigFloat,nt)   # Solution of the product π

            for μ = fi : li
                sol_prod_temp = zeros(BigFloat,nt)
                a = bc[ν+1][k][μ]   # Coefficients
                for (idx, _) in enumerate(tgrid)
                    @views x = g_der_data[idx,:]
                    sol_prod_temp[idx] = a * mapreduce(^, *, x, bp[ν+1][k][:,μ])
                end
                sol_prod = sol_prod + sol_prod_temp
            end

            q[:,ν] = q[:,ν] + big.(f_der_data[:,k]).*sol_prod
        end
    end 
    return q
end