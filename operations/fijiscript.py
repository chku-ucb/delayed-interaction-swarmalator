import sys
import os
from java.io import File
from ij import IJ
from fiji.plugin.trackmate import Model, Settings, TrackMate, Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.features import FeatureFilter

# --- CONFIGURATION START ---
# Updated with your specific paths
input_dir = "./tif_stacks" # path to directory containing the input TIF stack files
output_dir = "./xml_tracks" # path to directory where TrackMate XML results will be saved
extension = ".tif"

# Parameters (Matching your manual log)
spot_radius = 5.0
threshold = 6.6136
median_filter = True 
quality_filter_value = 8.5 

linking_max_dist = 50.0
gap_closing_max_dist = 15.0
max_frame_gap = 2

min_max_dist_traveled = 10.0
min_track_duration = 30.0
# --- CONFIGURATION END ---

if not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
    except OSError:
        pass # Directory likely exists or permissions issue

files = [f for f in os.listdir(input_dir) if f.endswith(extension)]

print("--- STARTING BATCH PROCESS ---")
print("Found " + str(len(files)) + " files.")

for filename in files:
    print("\nProcessing: " + filename)
    file_path = os.path.join(input_dir, filename)
    
    imp = IJ.openImage(file_path)
    if imp is None:
        print("Could not open: " + filename)
        continue

    # --- 1. ROBUST DIMENSION FIX ---
    dims = imp.getDimensions() # [width, height, nChannels, nSlices, nFrames]
    nChannels, nSlices, nFrames = dims[2], dims[3], dims[4]

    print("  Original Dims -> Slices (Z): " + str(nSlices) + ", Frames (T): " + str(nFrames))

    # If it looks like a Z-stack (many slices, 1 frame), FORCE it to be Time
    if nFrames == 1 and nSlices > 1:
        print("  -> FIXING: Converting Z-slices to Time-frames...")
        imp.setDimensions(nChannels, 1, nSlices)
        imp.setOpenAsHyperStack(True) 
        nFrames = nSlices 
    
    if imp.getNFrames() <= 1:
        print("  WARNING: Image still has 1 Frame! Tracking will fail.")
    else:
        print("  Confirmed: Image now has " + str(imp.getNFrames()) + " Frames.")

    # --- 2. SETUP TRACKMATE ---
    model = Model()
    model.setLogger(Logger.IJ_LOGGER)
    settings = Settings(imp)

    # Detector
    settings.detectorFactory = LogDetectorFactory()
    settings.detectorSettings = {
        'DO_SUBPIXEL_LOCALIZATION': True,
        'RADIUS': float(spot_radius),
        'TARGET_CHANNEL': 1,
        'THRESHOLD': float(threshold),
        'DO_MEDIAN_FILTERING': median_filter,
    }
    
    # Spot Filter
    settings.addSpotFilter(FeatureFilter('QUALITY', float(quality_filter_value), True))

    # Tracker
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
    settings.trackerSettings['LINKING_MAX_DISTANCE'] = float(linking_max_dist)
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = float(gap_closing_max_dist)
    settings.trackerSettings['MAX_FRAME_GAP'] = int(max_frame_gap)
    settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
    settings.trackerSettings['ALLOW_TRACK_MERGING'] = False

    # Track Filters
    settings.addTrackFilter(FeatureFilter('MAX_DISTANCE_TRAVELED', float(min_max_dist_traveled), True))
    settings.addTrackFilter(FeatureFilter('TRACK_DURATION', float(min_track_duration), True))

    # --- 3. RUN PROCESS ---
    trackmate = TrackMate(model, settings)
    
    if not trackmate.checkInput() or not trackmate.process():
        print("  Error: " + str(trackmate.getErrorMessage()))
        imp.close()
        continue
    
    # --- 4. FORCE FEATURE CALCULATION ---
    trackmate.computeSpotFeatures(True)
    trackmate.computeEdgeFeatures(True)
    trackmate.computeTrackFeatures(True)

    # --- 5. CHECK RESULTS ---
    # Using .getNSpots(True) to fix the previous error
    n_spots = model.getSpots().getNSpots(True) 
    n_tracks = model.getTrackModel().nTracks(True)
    
    print("  -> RESULT: Found " + str(n_spots) + " spots and " + str(n_tracks) + " filtered tracks.")

    if n_tracks == 0:
        print("  WARNING: 0 Tracks found. Check your parameters!")

    # --- 6. SAVE ---
    out_name = os.path.splitext(filename)[0] + ".xml"
    out_path = os.path.join(output_dir, out_name)
    
    java_file = File(out_path)
    writer = TmXmlWriter(java_file)
    writer.appendModel(model)
    writer.appendSettings(settings)
    writer.writeToFile()
    
    print("  Saved to: " + out_path)
    imp.close()

print("Batch processing complete.")
