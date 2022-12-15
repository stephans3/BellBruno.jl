using BellBruno

N_der = 10;              # Number of derivatives
bp = bell_poly(N_der);   # Create bell polynomials
bc = bell_coeff(bp);     # Compute bell coefficients

#=
    f(z) = -1 z^(-2)
    d^n/dz^n f(z) = -1 * z^(-2-n) π(-2-j) for j=0:1:n
=#
function outer_fun!(y, z) 
    c = -1;
    p = -2;

    y[1] = simple_monomial(z, c, p)

    fi = firstindex(y)
    li = lastindex(y)

    for idx in fi:li-1
        y[idx+1] = simple_monomial_der(z, c, p, idx)
    end
end

#=
    g(t)  = t/T - (t/T)^2
    g'(t) = 1/T - 2*t/T^2
    g''(t) = -2/T^2
    g^(n)(t) = 0 for n>2
=#
function inner_fun!(x, t :: Float64; T = 1.0 :: Real)
    c₁ = 0;
    c₂ = 1/T;
    c₃ = -1/T^2;

    x[1] = c₁ + c₂*t + c₃*t^2;  # g(t)  = t/T - (t/T)^2
    x[2] = c₂ + 2*c₃*t;         # g'(t) = 1/T - (2/T^2)*t
    x[3] = 2*c₃;                # g''(t)= -2/T^2
    x[4:end] .= 0;              # g^(n)(t) = 0 for n>2
end

T = 10;     # Final simulation time
dt = 1e-3   # Sampling time
tgrid = dt : dt : T-dt; # Time grid
nt = length(tgrid)      # Number of time steps

# Outer derivatives
g̃ = zeros(nt, N_der+1); # g̃_n(t) := d^n/dt^n g(t)
f̃ = zeros(nt, N_der+1); # f̃_n(t) := d^n/dy^n f(z)

for (idx, elem) in enumerate(tgrid)
    @views x = g̃[idx,:]
    @views y = f̃[idx,:]
    inner_fun!(x, elem, T=T)
    outer_fun!(y, x[1]) # Test if equals to second_fun!(y, g̃[idx,1])
end


#=
    1. Summation to find d^n/dt^n f(g(t)) =: q_n(t)
=#
q = zeros(nt, N_der);
for ν=1 : N_der
    for k=1 : ν
        fi = firstindex(bp[ν+1][k][1,:])
        li = lastindex(bp[ν+1][k][1,:])
        sol_prod = zeros(BigFloat,nt)   # Solution of the product π

        for μ = fi : li
            sol_prod_temp = zeros(BigFloat,nt)
            a = bc[ν+1][k][μ]   # Coefficients
            for (idx, _) in enumerate(tgrid)
                @views x = g̃[idx,:]
                sol_prod_temp[idx] = a * mapreduce(^, *, x[2:end], bp[ν+1][k][:,μ])
            end
            sol_prod = sol_prod + sol_prod_temp
        end

        q[:,ν] = q[:,ν] + f̃[:,k+1].*sol_prod
    end
end

#=
    2. Summation to find d^n/dt^n h(t) = d^n/dt^n exp(f(g(t))
=#
h = zeros(BigFloat,nt, N_der+1);
h[:,1] = exp.(big.(f̃[:,1]))
for ν=1 : N_der
    for k=1 : ν
        fi = firstindex(bp[ν+1][k][1,:])
        li = lastindex(bp[ν+1][k][1,:])
        sol_prod = zeros(BigFloat,nt)

        for μ = fi : li
            sol_prod_temp = zeros(BigFloat,nt)
            a = bc[ν+1][k][μ]   # Coefficients
            for (idx, _) in enumerate(tgrid)
                @views x = q[idx,:]
                sol_prod_temp[idx] = a * mapreduce(^, *, x, bp[ν+1][k][:,μ])
            end
            sol_prod = sol_prod + sol_prod_temp
        end
        h[:,ν+1] = h[:,ν+1] + exp.(big.(f̃[:,1])).*sol_prod
    end
end

using DelimitedFiles
path_q = "q_results.txt"
path_h = "h_results.txt"

open(path_q, "w") do io
    writedlm(io, q)
end

open(path_h, "w") do io
    writedlm(io, h)
end
