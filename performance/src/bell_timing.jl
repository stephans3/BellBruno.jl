#=
    Computing time for bell_poly
=#

using BellBruno
N_bell = 40;
N_rep = 20;
ctimes = zeros(N_rep, N_bell)

for nr=1:N_rep
    for nb=1:N_bell
        ctimes[nr, nb] = @elapsed bell_poly(nb);
    end
end

ngrid = 1:N_bell
ctimes_mean = sum(ctimes, dims=1)[:]/N_rep
#ctimes_std = sqrt.(sum((ctimes .- ctimes_mean').^2, dims=1) / N_rep )[:]

ctimes_max = maximum(ctimes, dims=1)[:]
ctimes_min = minimum(ctimes, dims=1)[:]

using CairoMakie
fig = Figure(fontsize=22)
ax1 = Axis(fig[1, 1], xlabel = "Number of derivatives n", ylabel="Computing time in [s]", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

ax1.limits=((0, 41), (-0.1,1.5)) 
ax1.xticks = vcat(1,5:5:40);
ax1.yticks= 0:0.2:1.4;
lines!(ngrid, ctimes_mean; linestyle = :dash, linewidth = 4, color = :navyblue, label = "Mean")
band!(ngrid, ctimes_min, ctimes_max; color = (:red, 0.2), label = "Deviation")
axislegend(ax1, position = :lt)
fig


basepath_im = "performance/results/images/"
path2fig = basepath_im*"performance_bell_polynomials.pdf"
save(path2fig, fig, pt_per_unit = 1)