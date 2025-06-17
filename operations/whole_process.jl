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


const N = 1000
const K = -0.75
const J = 1.0

function sim_(τ, T, id, N, K, J)
    """
    Calculate average fraction of defective agents over the last 100 time steps.
    Inputs:
        τ - time delay
        T - total simulation time
        N - number of agents (default 1000)
        K - coupling constant (default -0.75)
        J - another model parameter (default 1.0)
    Output:
        Average fraction of defective agents
    """
    Random.seed!(id)
    time, x, y, θ = SimulationSwarmalator(N, K, J, τ, T)
    return (time, x, y, θ)
end


function cal_def2(time, x, y, θ, tau, id)

    R = sqrt.((x .- mean(x, dims=2)) .^2 .+ (y .- mean(y, dims=2)) .^2)
    if tau < 9.5
        rbulk =  R_solid_core(x, y)
    else
        rbulk = R_bulk_long_tau(x, y)
    end
    
    p = scatter(x[end,:], y[end,:], dpi=300, aspect_ratio=:equal)
    plot_circle(mean(x[end,:]), mean(y[end,:]), rbulk, p)
    plot_circle(mean(x[end,:]), mean(y[end,:]), rbulkl, p, c=:blue)
    savefig("cir_tau_$(tau)_id_$(id).pdf")

    ndef, n6p, nbulk = count_defect(x,y,time,0.85*rbulk)

    return (ndef ./ nbulk , time)
end




function main()
    """
    This is the experiment to collect the defect fraction and 
    """

    # taus = [5.0,5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0 ,15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0]; 
    T = 70000.0; 
    Rd = 12; # Number of random seeds to use for averaging
    frac_all = zeros(Rd, length(taus));
    qfac_all = zeros(Rd, length(taus));
    for ir in 1:Rd
        id = rand(1:10000)
        for i in 1:length(taus)
            tau = taus[i]
            time, x, y, θ = sim_(tau, T, id, N, K, J)
            R = sqrt.((x .- mean(x, dims=2)) .^2 .+ (y .- mean(y, dims=2)) .^2)
            avgR = mean(R, dims=2)[:]
            fracs, times = cal_def2(time, x, y, θ, tau)
            itx = findfirst(x -> x>T/2, times)
            avg = mean(fracs[itx:end])
            stdf = std(fracs[itx:end])
            idxf = findfirst(x -> x< avg + stdf, fracs)

            idxT, params = find_γ_ω(avgR, time)
            println("Q-fac: ", abs(params[3])/(2*params[2]))
            try
                frac_all[ir, i] = mean(fracs[idxf:idxf+500])
            catch 
                frac_all[ir, i] = mean(fracs[itx:itx+500])
            end
            println("DefFrac: ", frac_all[ir,i])
            qfac_all[ir,i] = abs(params[3])/(2*params[2])
        end
    end
    ID = rand(1:500000)
    @save "qfac_frac_lowtau_batches_$(ID).jld2" taus frac_all qfac_all
end

main()
