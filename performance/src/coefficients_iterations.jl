using BellBruno, LinearAlgebra
N_bell = 40;
bp = bell_poly(N_bell);
bc = bell_coeff(bp);

T = diagm(ones(Int64,N_bell))
for n=2:N_bell
    for k=1:n-1
        T[n,k] = sum(T[n-k,1:k]);
    end
end

Tsum = sum(T, dims=2)[:]
Tsum[10:5:40]
