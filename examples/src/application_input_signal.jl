using DelimitedFiles
basepath = "examples/results/data/"
basepath_im = "examples/results/images/"
path_h = basepath * "h_results.txt"


hdata = readdlm(path_h, '\t', BigFloat)
T = 10;     # Final simulation time
dt = 1e-3   # Sampling time
tgrid = dt : dt : T-dt; # Time grid

Ω_int = sum(hdata[:,1])*dt

u_input = zeros(size(hdata)[1])

i = 0;
for col in eachcol(hdata)
    u_input += col / factorial(big(2i+1))
    i += 1
end

u_input = u_input / Ω_int

using CairoMakie
fig = Figure(fontsize=22)
ax1 = Axis(fig[1, 1], xlabel = "Time t in [s]", ylabel="Input signal", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

# xlims!(ax1, -0.5, T+0.5)
# ylims!(ax1, -3e-7, 40e-8)
ax1.limits=((-0.5, T+0.5), (-0.05, 0.75)) 
ax1.xticks = 0:2:T;
#ax1.yticks=-2e-7:1e-7:2e-7;
#ax1.yticks=0.0:2e-8:13e-8;
ax1.yticks=0.0:0.1:0.7;
lines!(tgrid, u_input; linewidth = 4)
fig
path2fig = basepath_im*"input_signal.pdf"
save(path2fig, fig, pt_per_unit = 1)