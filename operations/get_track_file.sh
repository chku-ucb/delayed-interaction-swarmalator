#!/bin/bash

# --- CONFIGURATION ---
# 1. Path to your Fiji executable (Double check this path!)
#    On macOS, it is usually inside /Applications/Fiji.app/Contents/MacOS/
FIJI_PATH="~/Fiji/Fiji.app/Contents/MacOS/fiji-macos-arm64" # path to Fiji executable

# 2. The Base Folder Path (without the number at the end)
BASE_PATH="./track_55_" # path to directory containing sequence of defect image plots

# 3. Where you want to save the output stacks
#    (Currently set to save in the SAME parent folder as create_pics, but you can change this)
OUTPUT_BASE_DIR="./tif_stacks" # path to directory where the output TIF stacks will be saved

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_BASE_DIR"

# --- THE LOOP (From 1 to 20) ---
for i in {1..20}
do
    # 1. DEFINE VARIABLES FOR THIS ITERATION (Don't forget these!)
    INPUT_DIR="${BASE_PATH}${i}"
    OUTPUT_FILE="${OUTPUT_BASE_DIR}/track_${i}_stack.tif"
    MACRO_FILE="temp_macro_${i}.ijm"

    # Optional: Print status to terminal so you know it's working
    echo "Processing set $i..."

    # 2. GENERATE THE MACRO
    cat <<EOT > "$MACRO_FILE"
    setBatchMode(true);

    input = "$INPUT_DIR/";
    output = "$OUTPUT_FILE";

    print("Importing sequence from: " + input);
    File.openSequence(input);

    print("Converting to 8-bit...");
    run("8-bit");

    print("Inverting...");
    // The "stack" argument prevents the popup error
    run("Invert", "stack");

    print("Saving stack...");
    saveAs("Tiff", output);

    close();
    print("Finished.");
    eval("script", "System.exit(0);");
EOT

    # 3. RUN FIJI
    "$FIJI_PATH" --headless -macro "$MACRO_FILE"

    # 4. CLEAN UP
    rm "$MACRO_FILE"
done

echo "========================================"
echo "All batches 9-20 complete."