# Adjust the path to environment and utils.jl as needed
script_dir = @__DIR__
code_root = joinpath(script_dir, "..") |> normpath
env_dir = joinpath(code_root, "env")
# -------------------------------------------------
using Pkg
Pkg.activate(env_dir)
include(joinpath(code_root,"utils.jl"))

using Random, Distributions, Statistics, DifferentialEquations, LsqFit, Plots
using Base.Threads

Random.seed!(1234)

# Define the function for the swarmalator model to use in delay differential equations
function swarmalator!(du, u, h, p, t)
    K, J, τ, N = p
    N = Int(N)
    x = u[1:N]
    y = u[N+1:2N]
    θ = u[2N+1:end]
    θd = h(p, t-τ)[2N+1:end]
    @threads for i in 1:N
        dxy(j) = sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2)
        du[i] = sum((1 + J*cos(θd[j] -θ[i]))*(x[j]-x[i])/dxy(j) - (x[j]-x[i])/dxy(j)^2 for j in 1:N if i != j)/N
        du[N+i] = sum((1 + J*cos(θd[j] -θ[i]))*(y[j]-y[i])/dxy(j) - (y[j]-y[i])/dxy(j)^2 for j in 1:N if i != j)/N
        du[2N+i] = K/N*sum(sin(θd[j]-θ[i])/dxy(j) for j in 1:N if i != j)
    end
end

# Define the function for running the simulation
function SimulationSwarmalator(N::Int64, K::Float64, J::Float64, τ::Float64, t_end::Float64)
    p = [K, J, τ, Int(N)]
    x= 2 .* rand(N) .- 1
    y= 2 .* rand(N) .- 1
    θ= 2π .* rand(N)
    u0 = [x; y; θ]
    tspan = (0.0, t_end)
    h(p, t) = u0
    prob = DDEProblem(swarmalator!, u0, h, tspan, p; constant_lags = τ)
    method = MethodOfSteps(Tsit5())
    sol = solve(prob, method, reltol=1e-12, abstol=1e-8)
    time = sol.t
    x = zeros(length(time), N)
    y = zeros(length(time), N)
    θ = zeros(length(time), N)
    for i in 1:length(time)
        x[i,:] = sol.u[i][1:N]
        y[i,:] = sol.u[i][N+1:2N]
        θ[i,:] = sol.u[i][2N+1:end]
    end
    return time, x, y, θ
end


function main()
    τs1 = [1.6, 1.7, 1.8, 1.9, 2.0]
    τs2 = [0.5*i + 2.0 for i in 1:36]
    #τs = vcat(τs1, τs2)
    τs = [10.0, 15.0, 20.0, 25.0]
    N = 100

    ωs = []
    γs = []

    for i in 1:length(τs)

        println("Running simulation for τ = ", τs[i])

        time, x, y, θ = SimulationSwarmalator(N, -0.75, 1.0, τs[i], 1000.0)

        R = sqrt.((x .- mean(x,dims=2)) .^2 .+ (y .- mean(y,dims=2)) .^2);
        avgR = mean(R,dims=2)[:]

        Err = 1000
        tol = 1e-6
        tidx = 0
        t_idx_end = findfirst(x->x>=700, time)

        # Model function
        @. function damped_oscillator(t, p)
            return p[1] * exp(-p[2]*t) *  sin(p[3] * t + p[4])
        end

        ωi = []
        γi = []
        while Err > tol
            # Initial parameters [A, gamma, omega, phi]
            p0 = [0.2, 0.01, 0.1, 0.1]  # Adjust these based on your data

            idxMx = findfirst(x->x>=tidx, time)
            # Your data
            t_data = range(0, t_idx_end-idxMx, length=t_idx_end-idxMx+1)
            x_data = avgR[idxMx:t_idx_end] .- avgR[end]

            # Fit the model
            fit = curve_fit(damped_oscillator, t_data, x_data, p0)

            # Extract fit parameters
            fit_params = fit.param

            pred = damped_oscillator(t_data, fit_params)

            #find mean squared error
            Err = mean((x_data .- pred).^2)
            tidx += 0.1
            push!(ωi, fit_params[3])
            push!(γi, fit_params[2])
        end

        push!(ωs, ωi[end])
        push!(γs, γi[end])
    end

    p = plot(xlabel="τ", ylabel="fit parametes",dpi=300)
    scatter!(p, τs, ωs, label="ω")
    scatter!(p, τs, γs, label="γ")
    savefig(p, "omg_lamd_tau.pdf")
end

main()