# Tutorial for BellBruno

using BellBruno

# Derivative of outer function
p=0.1;
f_der(x,n) = p^n * exp(p*x) 

# Inner function
g(t) = sin(t)

# Derivative of inner function
g_der(t,n) = sin(t+n*π/2) 

# Create Bell polynomial data
bp = bell_poly(10)

# Create Bell coefficient data
bc = bell_coeff(bp)

# Sampling points (time grid)
tgrid = -π:0.01:2π

# Compute derivatives
diff_data = faa_di_bruno(f_der, g, g_der, tgrid,bp, bc)

# Plotting 
using CairoMakie

basepath_im = "examples/results/images/"


begin
    fig = Figure(fontsize=22)
    ax1 = Axis(fig[1, 1], xlabel = "Sampling t", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

    ax1.limits=((-3.5, 7), (-0.25, 0.25)) 
    ax1.xticks = -3:1:7;
    ax1.yticks= -0.2:0.05:0.2;

    lines!(tgrid, diff_data[:,1]; linewidth = 4, label = L"1st")
    lines!(tgrid, diff_data[:,2]; linestyle = :dot, linewidth = 4, label = L"2nd")
    lines!(tgrid, diff_data[:,3]; linestyle = :dashdot, linewidth = 4,label = L"3rd")
    lines!(tgrid, diff_data[:,4]; linestyle = :dashdotdot, linewidth = 4,label = L"4th")
    lines!(tgrid, diff_data[:,5]; linestyle = :dash, linewidth = 4,label = L"5th")

    axislegend(; position = :lb, backgroundcolor = (:grey90, 0.1));
    fig
    # path2fig = basepath_im*"tutorial_1.png"
    # save(path2fig, fig, pt_per_unit = 1) 
end


begin
    fig = Figure(fontsize=22)
    ax1 = Axis(fig[1, 1], xlabel = "Sampling t", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)

    ax1.limits=((-3.5, 7), (-6.0, 4.5)) 
    ax1.xticks = -3:1:7;
    ax1.yticks= -6.:1.0:4.;

    lines!(tgrid, diff_data[:,6]; linewidth = 4, label = L"6th")
    lines!(tgrid, diff_data[:,7]; linestyle = :dot, linewidth = 4, label = L"7th")
    lines!(tgrid, diff_data[:,8]; linestyle = :dashdot, linewidth = 4,label = L"8th")
    lines!(tgrid, diff_data[:,9]; linestyle = :dashdotdot, linewidth = 4,label = L"9th")
    lines!(tgrid, diff_data[:,10]; linestyle = :dash, linewidth = 4,label = L"10th")

    axislegend(; position = :lb, backgroundcolor = (:grey90, 0.1));
    fig
    # path2fig = basepath_im*"tutorial_2.png"
    # save(path2fig, fig, pt_per_unit = 1) 
end

