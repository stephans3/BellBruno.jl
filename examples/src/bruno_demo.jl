

"""
    sin_der(x, N_der)

Computes the N_der derivatives of `sin(x)` at positions `x`.
"""
function sin_der(x, N_der)
    data = zeros(length(x),N_der+1)
    for n=0:N_der
        @. data[:,n+1] = sin(x + n*π/2) 
    end
    return data
end

"""
    inv_der(x, N_der)

Computes the N_der derivatives of `1/(x)` at positions `x`.
"""
function inv_der(x, N_der)
    data = zeros(length(x),N_der+1)
    for n=0:N_der
       @. data[:,n+1] = (-1)^n * factorial(n)*inv(x)^(1+n)
    end
    return data
end




function bruno(bp,bc, inner, outer)
    nt, N = size(inner)
    q = zeros(nt, N-1);
    for ν=1 : N-1
        for k=1 : ν
            fi = firstindex(bp[ν+1][k][1,:])
            li = lastindex(bp[ν+1][k][1,:])
            sol_prod = zeros(BigFloat,nt)   # Solution of the product π
    
            for μ = fi : li
                sol_prod_temp = zeros(BigFloat,nt)
                a = bc[ν+1][k][μ]   # Coefficients
                for idx in axes(inner,1)
                    @views x = inner[idx,:]
                    sol_prod_temp[idx] = a * mapreduce(^, *, x[2:end], bp[ν+1][k][:,μ])
                end
                sol_prod = sol_prod + sol_prod_temp
            end
    
            q[:,ν] = q[:,ν] + outer[:,k+1].*sol_prod
        end
    end
    return q
end

using BellBruno
N=10;
bp = bell_poly(N);
bc = bell_coeff(bp);     # Compute bell coefficients


x = 0:0.01:4;
data_in = sin_der(x, N)
data_out = inv_der(x+4,N)



function _bell_coeff12(mat :: Matrix{Int8},n)
    res = Int64[];
    for col in eachcol(mat)
        res = vcat(res, _bell_coeff(col[:],n))
    end
    return res;
end

function _bell_coeff12(v :: Vector{Int8},n)
    a1 = mapreduce(i-> factorial(big(v[i]))*factorial(big(i))^v[i],*, 1:n-1)
    a2 = factorial(big(n-1))*inv(a1)
    return round(BigInt, a2)
end



function _bell_coeff(mat :: Matrix{Int8},n)
    res = Int64[];
    for col in eachcol(mat)
        res = vcat(res, _bell_coeff(col[:],n))
    end
    return res;
end

function _bell_coeff(v :: Vector{Int8},n)
    a1 = mapreduce(i-> factorial(big(v[i]))*factorial(big(i))^v[i],*, 1:n-1)
    a2 = factorial(big(n-1))*inv(a1)
    return round(BigInt, a2)
end


n = 5
_bell_coeff(bp[n][1],n)
_bell_coeff(bp[n][2],n)
_bell_coeff(bp[n][3],n)




nn = 5
d = bp[nn][2][:,2]

_bell_coeff(bp[nn][2],nn)

a1 = mapreduce(i-> factorial(d[i])*factorial(i)^d[i],*, 1:nn-1)
a2 = factorial(nn-1)*inv(a1)
res = round(BigInt, a2)
bc[nn]

bruno(bp,bc,data_in, data_out)




