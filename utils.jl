using Plots
using DifferentialEquations
using Random
using JLD2
using LsqFit
using DelaunayTriangulation
using LinearAlgebra
using Statistics
using Distributions
using Base.Threads
using LaTeXStrings



function swarmalator!(du, u, h, p, t)
    """
    Function to compute the derivatives of the swarmalator system.
    This function will be used in ODE solvers
    du: array to store the derivatives
    u: current state of the system
    h: time step size
    p: parameters of the system
    t: current time
    """
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

function SimulationSwarmalator(N::Int64, K::Float64, J::Float64, τ::Float64, t_end::Any)
    """
    Simulates the swarmalator system using the DDEProblem from DifferentialEquations.jl.
    N: number of swarmalators
    K: coupling strength for the phase dynamics
    J: coupling strength for the spatial dynamics
    τ: delay in the system
    t_end: end time for the simulation
    Returns the time, x, y, and θ arrays of the swarmalators.
    This function initializes the swarmalators in a random position and phase,
    and then integrates the system using the Method of Steps.
    """
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

function radii(x,y,t,i)
    """
    Computes the radii of the swarmalators at time t for the i-th swarmalator.
    x: array of x-coordinates of the swarmalators
    y: array of y-coordinates of the swarmalators
    t: time index
    i: index of the swarmalator
    Returns the radius of the i-th swarmalator at time t.
    """
    return sqrt.((x[t,i] .- mean(x[t,:])).^2 .+ (y[t,i] .- mean(y[t,:])).^2)
end

function R_bulk_long_tau(x,y)
    """
    Compute the bulk radius of the swarmalators at the end of the simulation.
    x: array of x-coordinates of the swarmalators
    y: array of y-coordinates of the swarmalators
    Returns the bulk radius, which is the radius of the swarmalators at the end of the simulation.
    This function calculates the distance of each swarmalator from the mean position at the end of the simulation,
    sorts these distances, and returns the radius corresponding to the maximum change in distance in the last quarter of the sorted distances.
    x and y should be the last time step of the simulation.
    """
    N = size(x,2)
    rall = sqrt.((x[end,:] .- mean(x[end,:])).^2 .+ (y[end,:] .- mean(y[end,:])).^2)
    dfR = diff(sort(rall))
    rbulk = sort(rall)[argmax(dfR[Int(round(3N/4)):end])+Int(round(3N/4))]
    return rbulk
end

function R_solid_core(x, y)
    """
    Compute the solid core radius of the swarmalators at the end of the simulation.
    x: array of x-coordinates of the swarmalators
    y: array of y-coordinates of the swarmalators
    Returns the radius of the solid core, which is defined as the radius at which the standard deviation of the distances
    from the mean position is maximized in the last 5000 time steps.
    This function calculates the distance of each swarmalator from the mean position at the end of the simulation,
    sorts these distances, and computes the change in sorted distance for the last 5000 time steps.
    The index of the maximum change in distance is used to determine the solid core radius.
    x and y should be the last time step of the simulation.
    """

    R = sqrt.((x .- mean(x, dims=2)).^2 .+ (y .- mean(y, dims=2)).^2 )
    r_sorted = sort(R[end,:])

    N = 1000
    ncom= 5000
    dr = zeros(ncom, N)
    for i in 1:ncom
        r_sorted = sort(R[end-i,:])
        r_sorted2 = sort(R[end-i-1,:])
        dr[i,:] = r_sorted .- r_sorted2
    end
    std_dr = std(dr, dims=1)[:]

    idx_core = argmax(diff(std_dr[1:6:1000])[1:160])*6
    println("idx_core: ", idx_core)
    return sort(R[end,:])[idx_core]
end

function find_γ_ω(avgR, time, tol=1e-6)
    """
    Function to find the damping coefficient (γ) and frequency (ω) of the damped oscillator model
    avgR: average radius of the swarmalators
    time: time array corresponding to avgR
    tol: tolerance for the fitting error
    Returns the index of the maximum radius and the fitted parameters [A, γ, ω, φ, offset].
    The function fits a damped oscillator model to the average radius data using nonlinear curve fitting.
    The model is defined as:
    A * exp(-γ * t) * sin(ω * t + φ) + offset
    where A is the amplitude, γ is the damping coefficient, ω is the frequency, φ is the phase shift, and offset is a constant.
    The fitting process iteratively adjusts the parameters until the mean squared error is below the specified tolerance.
    """

    Err = 1000
    tidx = 0
    t_idx_end = findfirst(x->x>=700, time)

    # Model function
    @. function damped_oscillator(t, p)
        return p[1] * exp(-p[2]*t) *  sin(p[3] * t + p[4]) + p[5]
    end

    fit_params = [0.0, 0.0, 0.0, 0.0, 0.0]
    idxMx = 0
    while Err > tol
        # Initial parameters [A, gamma, omega, phi]
        p0 = [0.2, 0.01, 0.1, 0.1, 0.0]  # Adjust these based on your data

        idxMx = findfirst(x->x>=tidx, time)
        # Your data
        t_data = cumsum([0 ; diff(time[idxMx:t_idx_end])])
        x_data = avgR[idxMx:t_idx_end]

        # Fit the model
        fit = curve_fit(damped_oscillator, t_data, x_data, p0)

        # Extract fit parameters
        fit_params = fit.param

        pred = damped_oscillator(t_data, fit_params)

        #find mean squared error
        Err = mean((x_data .- pred).^2)
        tidx += 0.1
    end
    return idxMx, fit_params
end


function count_defect(x,y,time,rbulk)
    """
    Count the number of defects, 6-fold, and bulk particles in the swarmalator system.
    x: array of x-coordinates of the swarmalators
    y: array of y-coordinates of the swarmalators
    time: time array corresponding to x and y
    rbulk: radius threshold for bulk particles
    Returns three arrays: ndef, n6p, and nbulk.
    ndef: number of defects at each time step -- particles with 5 or 7 neighbours
    n6p: number of 6-fold particles at each time step
    nbulk: number of bulk particles at each time step
    This function uses the Delaunay triangulation to find the neighbours of each particle,
    and counts the number of neighbours for each particle at each time step.
    The number of neighbours is used to classify particles as defects, 6-fold, or bulk.
    """
    N = size(x,2)
    ndef = zeros(length(time))
    n6p = zeros(length(time))
    nbulk = zeros(length(time))
    for t in 1:length(time)
        pts = vcat(x[t,:]',y[t,:]')
        tri = triangulate(pts)
        nn = [length(get_neighbours(tri)[i]) for i in 1:N if radii(x,y,t,i)<rbulk]
        ndef[t] = count(x->x==5,nn)
        n6p[t] = count(x->x==6,nn)
        ndef[t] += count(x->x==7,nn)
        nbulk[t] = length(nn)
    end
    return ndef, n6p, nbulk
end


function orientation_order(tri, idx, x, y)
    """
    Compute the orientation order parameter for a given particle in the swarmalator system.
    tri: Delaunay triangulation of the swarmalator positions
    idx: index of the particle for which to compute the orientation order
    x: array of x-coordinates of the swarmalators
    y: array of y-coordinates of the swarmalators
    Returns the orientation order parameter for the particle at index idx.
    The orientation order parameter is defined as the average of the complex exponential of the angles
    between the particle and its neighbours, normalized by the number of neighbours.
    """
    N = length(get_neighbours(tri)[idx])
    psi = ComplexF64(0.)
    for i in collect(get_neighbours(tri)[idx])
        if i > 0
            θ = atan(y[i]-y[idx],x[i]-x[idx])
            psi += ℯ^(1im*6*θ)
        end
    end
    return psi/N
end