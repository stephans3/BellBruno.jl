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

begin
    fig = Figure(fontsize=22)
    ax1 = Axis(fig[1, 1], xlabel = "Time t in [s]", ylabelsize = 22,
        xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
        xtickalign = 1., xticksize = 10, 
        xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
        yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
        ytickalign = 1, yticksize = 10, xlabelpadding = 0)
    
    ax1.limits=((-0.5, T+0.5), (-3.0e-7, 2.0e-7)) 
    ax1.xticks = 0:2:T;
    ax1.yticks=-2e-7:1e-7:2e-7;
    
    lines!(tgrid, hdata[:,1]; linewidth = 4, label = L"\Omega(t)")
    lines!(tgrid, hdata[:,2]; linestyle = :dash, linewidth = 4, label = L"\dot{\Omega}(t)")
    lines!(tgrid, hdata[:,3]; linestyle = :dot, linewidth = 4,label = L"\ddot{\Omega}(t)")
    
    axislegend(; position = :lb, backgroundcolor = (:grey90, 0.1));
    fig
    path2fig = basepath_im*"example_1.pdf"
    save(path2fig, fig, pt_per_unit = 1) 
end


begin
    fig2 = Figure(fontsize=22)
ax2 = Axis(fig2[1, 1], xlabel = "Time t in [s]", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

    ax2.limits=((-0.5, T+0.5), (-5.2e-3, 4.5e-3)) 
    ax2.xticks = 0:2:T;
    ax2.yticks=-4e-3:1e-3:4e-3;

    lines!(tgrid, hdata[:,9]; linewidth = 3, label = L"\Omega^{(8)}(t)")
    lines!(tgrid, hdata[:,10]; linestyle = :dash, linewidth = 4, label = L"\Omega^{(9)}(t)")
    lines!(tgrid, hdata[:,11]; linestyle = :dot, linewidth = 4,label = L"\Omega^{(10)}(t)")

    axislegend(; position = :lb, backgroundcolor = (:grey90, 0.1));
    fig2
    path2fig = basepath_im*"example_2.pdf"
    save(path2fig, fig2, pt_per_unit = 1)
end


