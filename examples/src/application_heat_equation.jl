
L = 1.; # Length of 1D rod

# Material
λ = 1; # Thermal conductivity
α = 1; # Diffusivity

diff_ref = 1; # (y_f - y_0) = 100 Kelvin

using DelimitedFiles
basepath = "examples/results/data/"
basepath_im = "examples/results/images/"
path_h = basepath * "h_results.txt"


hdata = readdlm(path_h, '\t', BigFloat)
T = 10;     # Final simulation time
dt = 1e-3   # Sampling time
tgrid = dt : dt : T-dt; # Time grid

Ω_int = sum(hdata[:,1])*dt
u_data = zeros(size(hdata)[1])

i = 0;
for col in eachcol(hdata)
    u_data += col / factorial(big(2i+1))
    i += 1
end

u_data = u_data / Ω_int
# Simulation
u_data = vcat(0, u_data, 0)

function input_signal(t,dt)
    if t <= 0
        return u_data[1]
    elseif t >= Tf
        return u_data[end]
    end
    τ = t/dt + 1
    t0 = floor(Int, τ)
    t1 = t0 + 1;

    u0 = u_data[t0]
    u1 = u_data[t1]

    a = u1-u0;
    b = u0 - a*t0

    return a*τ + b;
end


# Diffusion: x-direction
function diffusion_x!(dx,x,Nx, Ny, Δx) # in-place
    
    for iy in 1 : Ny
        for ix in 2 : Nx-1
            i = (iy-1)*Nx + ix
            dx[i] =  (x[i-1] - 2*x[i] + x[i+1])/Δx^2
        end
        i1 = (iy-1)*Nx + 1      # West
        i2 = (iy-1)*Nx + Nx     # East
        dx[i1] = (-2*x[i1] + 2*x[i1+1])/Δx^2
        dx[i2] = (2*x[i2-1] - 2*x[i2])/Δx^2
    end

    nothing 
end


# 1D heat equation
function heat_eq!(dx,x,p,t)       
    # time = t/Tf;
    #u = input_signal(time, p)
    
    u = input_signal(t,ts)
    
    diffusion_x!(dx,x,Nx,1,Δx)
  
    dx .= α * dx
    dx[1] = dx[1] + 2α/(λ * Δx) * u
end


# Discretization  
const Nx = 101;     # Number of elements x-direction
const Δx = L/(Nx-1) # Spatial sampling
const Tf = T;  # Final time
const ts = T / (length(u_data)-1) # 1.0;     # Time step width

# Simulation without optimization
using OrdinaryDiffEq

x0 = zeros(Nx) # 300 * ones(Nx) # Intial values
tspan =  (0.0, Tf)   # Time span

prob = ODEProblem(heat_eq!,x0,tspan)
#alg = KenCarp4()    # Numerical integrator
#sol = solve(prob,alg, saveat = dt)

alg = Euler()    # Numerical integrator
sol = solve(prob,alg,dt=0.2*Δx^2, saveat = ts)


using CairoMakie
tgrid = sol.t;
fig1 = Figure(fontsize=22)
ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel="Temperature", ylabelsize = 22,
    xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
    xtickalign = 1., xticksize = 10, 
    xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
    yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
    ytickalign = 1, yticksize = 10, xlabelpadding = 0)
ax1.limits=((-0.5, T+0.5), (-0.05, 1.05)) 
ax1.xticks = 0:2:T;
ax1.yticks = 0 : 0.2 : 1.0;
lines!(tgrid, sol[end,:];linewidth = 4, label = "x=1.0 m (right side)")
#axislegend(; position = :lt, bgcolor = (:grey90, 0.1));
fig1
path2fig = basepath_im*"heat_equation.pdf"
save(path2fig, fig1, pt_per_unit = 1)



