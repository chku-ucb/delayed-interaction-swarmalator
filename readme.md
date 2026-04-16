# Aging in a Two-Dimensional Delayed Swarmalator Crystal

This is a repository for code that use in the article with aforementioned title. This code is written in Julia. 

> [!WARNING]
> Please note that all of the code in this repository was run and tested specifically on a macOS ARM M-Series CPU environment. Depending on your operating system or architecture, some configurations (especially the paths for Fiji/ImageJ or multi-threading) may require manual adjustments.

## Preparations

1. Install Julia from [julialang.org](https://julialang.org/downloads/).
2. Open Julia and install the required packages. After cloning the repository, run the following commands in the Julia REPL:
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

3. To perform defect tracking, download the Fiji program from the [Fiji Installation Page](https://imagej.net/software/fiji/downloads).

## Running the Code
To run the code you can use the example scripts provided in the `examples` directory with the following command:
```bash
export JULIA_NUM_THREADS=10
julia examples/example_script.jl
```
This is the recommended way to run the simulation because it allows you to use multiple threads for faster computation. You can change the number of threads by modifying the `JULIA_NUM_THREADS` environment variable.

You can check the number of threads available in your Julia installation by running:
```julia   
using Base.Threads
Threads.nthreads()
```


## Examples

- `animate_system.jl`: this script simulates the system with a given set of parameters and number of particles, and saves the results to an mp4 file.
- `hexatic_order.jl`: this script simulates the system with a given set of parameters and number of particles, and calculates the hexatic order parameter for the system. Consequently, it displays the results in animation and saves the results to a mp4 file.
- `defect_in_bulk.jl`: this script simulates the system, finds the bulk radius of the system and counts the number of defects in the system. For experiments where $$\tau \le \tau_c$$, the script includes an option to switch the radius metric function from 'R_bulk_long_tau' to 'R_solid_core' to capture the solid core instead of bulk radius.


## Operation

This section is a brief overview of the process that are used for finding the defect fraction for different time delays.

The `operations` directory contains a pipeline of scripts that automate the entire workflow:

1. **`gen_data.jl`**: A sample Julia script used to generate and save the raw Swarmalator model simulation data (such as positions and phases) to JLD2 files.
2. **`whole_process.jl`**: A comprehensive simulation script that runs the swarmalator model across various time delays (`tau`), calculates defect fractions over time, and stores the results in JLD2 formats.
3. **`get_figs_all.jl`**: Reads the resulting JLD2 simulation files and generates sequences of images (plots) displaying the defects at different timestamps.
4. **`get_track_file.sh`**: A shell script generating Fiji/ImageJ macros to convert the generated image sequences into 8-bit inverted TIF stacks.
5. **`fijiscript.py`**: A Python script intended to be run via Fiji's TrackMate plugin in headless mode. It processes the TIF stacks to detect and track the spots/defects, outputting tracking coordinates as XML files.
6. **`generic_bash.sh`**: A master bash script that ties the tracking pipeline together—it consecutively runs `get_figs_all.jl`, executes `get_track_file.sh`, runs the Fiji tracker (`fijiscript.py`), and analyzes the track crossings.

### Running the Tracking Pipeline

To extract the track files and calculate the particle crossings over the set radius, you can simply execute the master `generic_bash.sh` script from your terminal:

```bash
cd operations
chmod +x generic_bash.sh get_track_file.sh
./generic_bash.sh
```

This single command will sequentially handle the workflow:
1. Generate the figure sequences.
2. Convert the figure sequences into `.tif` stacks for TrackMate.
3. Perform particle tracking to acquire trajectories inside XML files.
4. Run `count_track_crossings.jl` to evaluate the particle trajectories crossing over your specified set radius.

#### Important: Setting the Processing Paths

Before executing the pipeline, you must ensure your paths are configured accurately to match your local setup. The following scripts have vital path variables:

- `generic_bash.sh` & `get_track_file.sh`: Set `FIJI_PATH` to point to your local ImageJ/Fiji executable (e.g., `~/Fiji/Fiji.app/Contents/MacOS/ImageJ`).
- `get_track_file.sh`, `merge_tif.sh` & `generic_bash.sh`: Verify that `BASE_PATH` points correctly to your input defect plot sequence directories (e.g., `./track_55_`) and `OUTPUT_BASE_DIR` maps to your output directory for `.tif` stacks.
- `fijiscript.py`: Ensure `input_dir` accurately points to your compiled sequence `.tif` stacks and `output_dir` denotes where to construct tracking `.xml` files.

*(Note: These variables exist near the top lines of each respective file and are pre-annotated with inline comments for convenience.)*
