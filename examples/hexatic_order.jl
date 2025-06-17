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


# Simple simulation with hexatic order parameter color labeling along with
# the phase of order parameter


function main()
    N = 100 # Number of swarmalators
    K = -0.75 # Coupling strength for the phase dynamics
    J = 1.0 # Coupling strength for the spatial dynamics
    τ = 20.0 # Delay in the system
    t_end = 500.0 # End time for the simulation

    Random.seed!(1234) # Set a random seed for reproducibility
    time, x, y, θ = SimulationSwarmalator(N, K, J, τ, t_end)

    # Calculate the velocity 
    x0 = x .- mean(x, dims=2)
    y0 = y .- mean(y, dims=2)

    dx = diff(x0, dims=1)
    dy = diff(y0, dims=1)
    dt = diff(time)
    
    vx = dx ./ dt
    vy = dy ./ dt

    # calculate the hexatic order parameter
    ϕs = zeros(length(time),N)
    psi6 = zeros(length(time),N)

    for t in 1:length(time)
        pts = vcat(x0[t,:]', y0[t,:]')
        tri = triangulate(pts)
        for i in 1:N
            ψ₆ = orientation_order(tri, i, x0[t,:], y0[t,:])
            ϕs[t,i] = angle(ψ₆)
            psi6[t,i] = abs(ψ₆)
        end
    end

    anim = @animate for i in 1:length(time)-1
        p = plot(xticks=false, yticks=false, axis = false, legend=false, dpi=200,grid=false, title = "idx = $(i)")
        scatter!(p,x0[i,:], y0[i,:], aspect_ratio=:equal, makersize = 1, markerstrokewidth=0,zcolor=psi6[i,:],colorbar=true,clim=(0.0, 1.0), cmap=cgrad(:coolwarm, rev=true))
        nvx = 2 .* normalize(vx[i,:])
        nvy = 2 .* normalize(vy[i,:])
        quiver!(p, x0[i,:], y0[i,:], quiver=(nvx, nvy), color=:black, lw=0.5)
        p2 = plot(xticks=false, yticks=false, axis = false, legend=false, dpi=200,grid=false, title = "idx = $(i)", aspect_ratio=:equal)
        nvx = 0.25 *cos.(ϕs[i,:])
        nvy = 0.25 *sin.(ϕs[i,:])
        quiver!(p2, x0[i,:], y0[i,:], quiver=(nvx, nvy), color=:black, lw=0.5)
        plot(p, p2, layout=(1,2), margin = 10Plots.mm, size = (1000, 500))
    end
    mp4(anim, "phase_orientation.mp4", fps = 15)
end

main()

