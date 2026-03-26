# Adjust the path to environment and utils.jl as needed
script_dir = @__DIR__
code_root = joinpath(script_dir, "..") |> normpath
env_dir = joinpath(code_root, "env")
# -------------------------------------------------
using Pkg
Pkg.activate(env_dir)
include(joinpath(code_root,"utils.jl"))

using Plots, DifferentialEquations, Distributions, Random, JLD2, LsqFit, DelaunayTriangulation
using Base.Threads

"""
This script serves as an entry point for generating and saving the raw simulation data.
It executes the Swarmalator model simulation for a predetermined set of parameters 
(such as time delay τ) and saves the resulting position and phase dynamics (time, x, y, θ) 
into a JLD2 file for further downstream processing or defect tracking.
"""


println("Num Threads: ", Threads.nthreads())

function main(τ)
    N = 1000 # Number of particles
    K = 0.5 # Phase coupling strength
    J = -0.5 # Spatial Couplign strength
    p = (K, J, τ, N)
    t_end = 200000.0 #Time

    id = 1234
    Random.seed!(id)
    time, x, y, θ = SimulationSwarmalator(N, K, J, τ, t_end)
    @save "sim_τ_$(τ)_id_$(id).jld2" time x y θ 
    
end

println("Running simulations for different τ values...")
println("τ = 20.0")
main(20.0)



