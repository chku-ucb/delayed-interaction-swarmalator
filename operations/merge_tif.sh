#!/bin/bash

# --- 1. SETUP CLUSTER ENVIRONMENT ---
# module load fiji/1.51

# Define your base paths
BASE_PATH="./track_55_" # path to directory containing sequence of defect image plots
OUTPUT_BASE_DIR="./tif_stacks" # path to directory where the output TIF stacks will be saved

# Create output folder if it doesn't exist
mkdir -p "$OUTPUT_BASE_DIR"

# --- 2. THE LOOP (From 1 to 100) ---
for i in {1..100}
do
    # Define variables for this specific track
    INPUT_DIR="${BASE_PATH}${i}"
    OUTPUT_FILE="${OUTPUT_BASE_DIR}/track_${i}_stack.tif"
    MACRO_FILE="temp_macro_${i}.ijm"

    # Only proceed if the input directory exists
    if [ -d "$INPUT_DIR" ]; then
        echo "------------------------------------------"
        echo "Processing track__$i..."

        # --- 3. GENERATE THE MACRO ---
        # Note: We use File.openSequence(input, " filter=defects_"); to ensure 
        # it only grabs your "defects_" files.
        cat <<EOT > "$MACRO_FILE"
        setBatchMode(true);
        input = "$INPUT_DIR/";
        output = "$OUTPUT_FILE";

        print("Importing sequence from: " + input);
        // "sort" handles the 1, 2, 10 numerical order issue
        run("Image Sequence...", "open=[" + input + "] file=defects_ sort");

        print("Converting to 8-bit...");
        run("8-bit");

        print("Inverting...");
        run("Invert", "stack");

        print("Saving stack...");
        saveAs("Tiff", output);

        close();
        eval("script", "System.exit(0);");
EOT

        # --- 4. RUN FIJI HEADLESS ---
        # 'ImageJ-linux64' is the standard command name after 'module load fiji'
        ImageJ-linux64 --headless -macro "$MACRO_FILE"

        # --- 5. CLEAN UP ---
        rm "$MACRO_FILE"
    else
        echo "Directory $INPUT_DIR not found, skipping..."
    fi
done

echo "Done! All files are in $OUTPUT_BASE_DIR"