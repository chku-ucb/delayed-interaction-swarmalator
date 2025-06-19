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
using DelaunayTriangulation
using LinearAlgebra
using Statistics
using Distributions
using Base.Threads
using LaTeXStrings
using JLD2


# Simple simulation to find the defect inside of defined bulk radius
# Also compare with the q-factor



function main()
    N = 500 # Number of swarmalators
    K = -0.75 # Coupling strength for the phase dynamics
    J = 1.0 # Coupling strength for the spatial dynamics
    τ = 20.0 # Delay in the system
    t_end = 2000.0 # End time for the simulation

    Random.seed!(1234) # Set a random seed for reproducibility
    time, x, y, θ = SimulationSwarmalator(N, K, J, τ, t_end)

    # Find the bulk radius

    R = sqrt.((x .- mean(x, dims=2)) .^2 .+ (y .- mean(y, dims=2)) .^2)
    avgR = mean(R, dims=2)[:]
    rbulk = R_bulk_long_tau(x, y) # Calculate the bulk radius
    ndef, n6p, nbulk = count_defect(x,y,time,0.85*rbulk) # using 85% of the bulk radius 
    fracs = ndef ./ nbulk
    itx = findfirst(x -> x>t_end/2, time)
    avg = mean(fracs[itx:end])
    stdf = std(fracs[itx:end])
    idxf = findfirst(x -> x< avg + stdf, fracs)
    idxT, params = find_γ_ω(avgR, time)
    println("Fraction of defects: $(avg) ± $(stdf)")
    println("Defect fraction at time $(time[idxf]): $(fracs[idxf])")

    plot(time, fracs, 
        xlabel = "Time", 
        ylabel = "Fraction of defects",
        title = "Defect fraction over time",
        label = "Defect fraction",
        legend = :topright,
        dpi = 200,
        xlims = (0, t_end),
        ylims = (0, 1),
        color = :blue)

    savefig("defect_fraction.png")

end

main()