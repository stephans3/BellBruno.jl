using DelimitedFiles
basepath = "examples/results/data/"
basepath_im = "examples/results/images/"
path_q = basepath * "q_results.txt"
path_h = basepath * "h_results.txt"


hdata = readdlm(path_h, '\t', BigFloat)
T = 10;     # Final simulation time
dt = 1e-3   # Sampling time
tgrid = dt : dt : T-dt; # Time grid

using CairoMakie
# fig = Figure(resolution = (600, 400), fonts = (; regular= "CMU Serif")) ## probably you need to install this font in your system
# size_inches = (4,3)
# size_pt = 72 .* size_inches
# fig = Figure(resolution=size_pt, fontsize=12) ## probably you need to install this font in your system
fig = Figure(fontsize=22)
ax1 = Axis(fig[1, 1], xlabel = "Time t in [s]", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

# xlims!(ax1, -0.5, T+0.5)
# ylims!(ax1, -3e-7, 40e-8)
ax1.limits=((-0.5, T+0.5), (-3.0e-7, 2.0e-7)) 
ax1.xticks = 0:2:T;
ax1.yticks=-2e-7:1e-7:2e-7;

lines!(tgrid, hdata[:,1]; linewidth = 4, label = L"\Omega(t)")
lines!(tgrid, hdata[:,2]; linestyle = :dash, linewidth = 4, label = L"\dot{\Omega}(t)")
lines!(tgrid, hdata[:,3]; linestyle = :dot, linewidth = 4,label = L"\ddot{\Omega}(t)")

# lines!(tgrid, hdata[:,1]; linewidth = 3, label = L"\Omega(t)=\exp(f~(g~(t~))")
# lines!(tgrid, hdata[:,2]; linestyle = :dash, linewidth = 3, label = L"\dot{\Omega}(t)=\frac{d}{d t} \left[\exp(f~(g~(t~))\right]")
# lines!(tgrid, hdata[:,3]; linestyle = :dot, linewidth = 3,label = L"\ddot{\Omega}(t)=\frac{d^{2}}{d t^{2}} \left[\exp(f~(g~(t~))\right]")

axislegend(; position = :lb, bgcolor = (:grey90, 0.1));
fig
path2fig = basepath_im*"example_1.pdf"
save(path2fig, fig, pt_per_unit = 1)


using CairoMakie
# fig = Figure(resolution = (600, 400), fonts = (; regular= "CMU Serif")) ## probably you need to install this font in your system
fig2 = Figure(fontsize=22)
ax2 = Axis(fig2[1, 1], xlabel = "Time t in [s]", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)
# xlims!(ax1, -0.5, T+0.5)
# ylims!(ax1, -3e-7, 40e-8)
ax2.limits=((-0.5, T+0.5), (-5.2e-3, 4.5e-3)) 
ax2.xticks = 0:2:T;
ax2.yticks=-4e-3:1e-3:4e-3;

lines!(tgrid, hdata[:,9]; linewidth = 3, label = L"\Omega^{(8)}(t)")
lines!(tgrid, hdata[:,10]; linestyle = :dash, linewidth = 4, label = L"\Omega^{(9)}(t)")
lines!(tgrid, hdata[:,11]; linestyle = :dot, linewidth = 4,label = L"\Omega^{(10)}(t)")

# lines!(tgrid, hdata[:,9]; linewidth = 3, label = L"\Omega^{(8)}(t) = \frac{d^{8}}{d t^{8}} \Omega(t)")
# lines!(tgrid, hdata[:,10]; linestyle = :dash, linewidth = 3, label = L"\Omega^{(9)}(t) = \frac{d^{9}}{d t^{9}} \Omega(t)")
# lines!(tgrid, hdata[:,11]; linestyle = :dot, linewidth = 3,label = L"\Omega^{(10)}(t) = \frac{d^{10}}{d t^{10}} \Omega(t)")
axislegend(; position = :lb, bgcolor = (:grey90, 0.1));
fig2
path2fig = basepath_im*"example_2.pdf"
save(path2fig, fig2, pt_per_unit = 1)


#axislegend(L"f(x)"; position = :rt, bgcolor = (:grey90, 0.25));
