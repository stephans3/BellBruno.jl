using BellBruno

N_der = 25;              # Number of derivatives
bp = bell_poly(N_der);   # Create bell polynomials
bc = bell_coeff(bp);     # Compute bell coefficients

bp_num_elem = zeros(Int64,N_der+1)
for (idx, elem) in enumerate(bp)
    bp_num_elem[idx] = size(hcat(elem...))[2]
end

coeffs = vcat(bc[N_der+1]...);

using Plots
scatter(bp_num_elem, legend=false, title="Number of monomials per Bell polynomial", xlabel="Index of Bell polynomial")

coeff_title = string("Coefficient value of Bell polynomial ", N_der)
bar(coeffs, legend=false, title=coeff_title, xlabel="Monomial index")