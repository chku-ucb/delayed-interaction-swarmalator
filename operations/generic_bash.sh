#!/bin/bash

# 1. get all figs
julia get_figs_all.jl

# 2. run Fiji
./get_track_file.sh

# 3. Path to your Fiji executable (Double check this path!)
#    On macOS, it is usually inside /Applications/Fiji.app/Contents/MacOS/
FIJI_PATH="~/Fiji/Fiji.app/Contents/MacOS/fiji-macos-arm64" # path to Fiji executable

# # # 4. The Base Folder Path (without the number at the end)
BASE_PATH="./track_55_" # path to directory containing sequence of defect image plots

# # # 5. Where you want to save the output stacks
OUTPUT_BASE_DIR="./tif_stacks" # path to directory where the output TIF stacks will be saved

# # # 6. The Fiji Script
FIJI_SCRIPT="./fijiscript.py" # path to the Python script for Fiji TrackMate

"$FIJI_PATH" --headless --console --run "$FIJI_SCRIPT"

julia count_track_crossings.jl