# Aging in a Two-Dimensional Delayed Swarmalator Crystal

This is a repository for code that use in the article with aforementioned title. This code is written in Julia. 

## Preparations

1. Install Julia from [julialang.org](https://julialang.org/downloads/).
2. Open Julia and install the required packages. After clone the repository, run the following commands in the Julia REPL:
```julia
using Pkg
Pkg.activate("env")
Pkg.update()
Pkg.instantiate()
```
If you want to use the Jupyter notebook, you also need to install the `IJulia` package:
```julia
Pkg.add("IJulia")
```

## Running the Code
To run the code you can use the example scripts provided in the `examples` directory with the following command:
```bash
export JULIA_NUM_THREADS=4
julia examples/example_script.jl
```
This is the recommended way to run the simulation because it allows you to use multiple threads for faster computation. You can change the number of threads by modifying the `JULIA_NUM_THREADS` environment variable.

You can check the number of threads available in your Julia installation by running:
```julia   
using Base.Threads
Threads.nthreads()
```


## Examples

- `animate_system.jl`: this script simulates the system with a given set of parameters and number of particles, and saves the animation of the system to mp4 file.
- `hexatic_order.jl`: this script simulates the system with a given set of parameters and number of particles, and calculates the hexatic order parameter for the system. Consequently, it displays the results in animation and saves the results to a mp4 file.
- `defect_in_bulk.jl`: this script simulates the system, find the bulk radius of the system and count the number of defects in the system. It also displays the results in plot. In addition, you can change the function to estimate the bulk radius to be `R_solid_core` to find the solid core radius of the system with $$\tau \le 9.0$$. 


## Operation 

This section is a brief overview of the process that use to find the defect fraction for each time delays.
