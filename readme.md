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


## 
