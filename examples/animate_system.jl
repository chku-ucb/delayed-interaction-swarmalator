# Adjust the path to environment and utils.jl as needed
script_dir = @__DIR__
code_root = joinpath(script_dir, "..") |> normpath
env_dir = joinpath(code_root, "env")
# -------------------------------------------------
using Pkg
Pkg.activate(env_dir)
include(joinpath(code_root,"utils.jl"))
using Random
using Plots
using DifferentialEquations
using LaTeXStrings
using JLD2


# Check the number of threads available
# If you want to speed up the simulation, you can set the number of threads to use
# 'export JULIA_NUM_THREADS=4' in your terminal before running this script
# number 4 here is just an example, you can set it to any number you want

#simple simulation with phase color labeling
function main()
    N = 100 # Number of swarmalators
    K = -0.75 # Coupling strength for the phase dynamics
    J = 1.0 # Coupling strength for the spatial dynamics
    τ = 20.0 # Delay in the system
    t_end = 10.0 # End time for the simulation

    Random.seed!(1234) # Set a random seed for reproducibility
    time, x, y, θ = SimulationSwarmalator(N, K, J, τ, t_end)

    anim = @animate for i in 1:length(time)
        p = plot(xticks = false, 
                yticks =false, 
                axis = false, 
                legend=false, 
                dpi=200, 
                title="time: $(round(time[i], digits=2))")
        scatter!(p, 
            x[i,:], 
            y[i,:], 
            aspect_ratio =:equal, 
            markersize = 1, 
            markerstrokewidth = 0.0,
            zcolor = θ[i,:] .% 2π,
            colorbar=false,
            cmap=cgrad(:lightrainbow, rev=true))
    end
    mp4(anim, "swarmalator_animation.mp4", fps=30)
end

main()