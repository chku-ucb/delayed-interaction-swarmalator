# Adjust the path to environment and utils.jl as needed
script_dir = @__DIR__
code_root = joinpath(script_dir, "..") |> normpath
env_dir = joinpath(code_root, "env")
# -------------------------------------------------
using Pkg
Pkg.activate(env_dir)
include(joinpath(code_root,"utils.jl"))

using Plots, JLD2, Statistics, LaTeXStrings, LinearAlgebra, Random, DelaunayTriangulation
using Base.Threads

"""
This script processes the raw JLD2 simulation files (containing x, y, θ, t) to 
generate a temporal sequence of images (plots) showcasing structural defects. 
These sequences are specifically pre-processed to be compatible with TrackMate 
in Fiji/ImageJ for converting into TIF tracking stacks. The script also evaluates 
and exports the defect fractions over time.

INPUT: A text file ('name_files.txt') containing names of simulation JLD2 files.
OUTPUT: 1. Sequence of defect figure PNGs grouped in directories.
        2. Defect fraction count data stored in corresponding JLD2 files.
"""

# load name from file: name_20_files.txt
names = readlines("name_files.txt") #list of name of files .jld2 that save the simulation results

for n in 1:length(names)

    data = jldopen(names[n])
    time = data["time"]
    x = data["x"]
    y = data["y"]
    
    id = n
    N = 1000
    x0 = x .- mean(x, dims=2)
    y0 = y .- mean(y, dims=2)

    # rbulk =  R_bulk_long_tau(x0,y0)
    rbulk = R_solid_core(x0,y0)

    colors = Dict(5 => :tomato, 7 => :blue, 6 => :lightcyan3) 
    tinit= findfirst(x -> x>1000, time)
    tend = findfirst(x -> x>20000, time)
    Ndef = []
    T = []
    println("length step: ", (tend - tinit) / 20)
    cp = ["black" for i in 1:N]
    idxt = tinit:tend-2

    for t in 1:20:length(idxt)-2
        pts = vcat(x0[idxt[t],:]', y0[idxt[t],:]')
        tri = triangulate(pts)
        nn = [length(get_neighbours(tri)[i]) for i in 1:N]
        xp = []
        yp = []
        ndff = 0
        for i in 1:N
            
            if nn[i] != 6
                push!(xp, x0[idxt[t],i])
                push!(yp, y0[idxt[t],i])
                if radii(x0, y0, idxt[t], i) < 0.85*rbulk
                    ndff += 1
                end
            end
        end
        push!(Ndef, ndff)
        push!(T, time[idxt[t]])
        
        p = plot(aspect_ratio=:equal, dpi=100, legend=false,xlim=(-rbulk*1.2, rbulk*1.2), ylim=(-rbulk*1.2, rbulk*1.2))
        scatter!(p, xp, yp, c = cp, aspect_ratio=:equal, makersize = 1, markerstrokewidth=0,xticks=false, yticks=false, axis = false)
        savefig(p, "track_55_$(n)/defects_$(round(Int, time[idxt[t]])).png")
        @save "tracks_low/defect_count_tau_55_$(id).jld2" Ndef T
    end
    
end

