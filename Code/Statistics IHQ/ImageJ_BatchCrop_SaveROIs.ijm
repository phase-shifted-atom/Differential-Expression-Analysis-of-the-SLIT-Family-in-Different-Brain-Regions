// ============================================================================
// Macro: ImageJ_BatchCrop_SaveROIs.ijm
// ----------------------------------------------------------------------------
// PURPOSE:
// - Batch-crop brightfield IHC images to a *fixed, standard region of interest*
//   (ROI) so that all subsequent quantification uses exactly the same field of view.
// - This ensures that between-slide comparisons (e.g. across subjects or groups)
//   are not confounded by differences in which anatomical area was sampled.
// - The macro also provides a simple way to visually document and reproduce
//   the cropping stage in the image analysis pipeline.
//
// CONTEXT IN THE PIPELINE:
// - Typically used after background harmonisation / pre-cleaning of images.
// - The cropped images are then passed to colour deconvolution and ROI-based
//   measurement macros to obtain quantitative IHC metrics (area, intensity, etc.).
//
// HOW (HIGH-LEVEL):
// 1. Ask the user for an input folder and an output folder.
// 2. Optionally let the user fine-tune the fixed ROI and options at run time.
// 3. Recursively (or non-recursively) scan the input folder for .tif/.tiff files.
// 4. For each image: check that the ROI rectangle is fully inside the image,
//    crop, and save to the output folder with a configurable suffix.
// 5. Maintain a progress counter and print informative messages for skipped files.
// ============================================================================
// === Batch Crop TIF/TIFF with a Fixed ROI ===
// - Prompts for input/output folders
// - Uses your rectangle: x=864, y=1296, w=3456, h=2556
// - Saves out as .tif with a suffix (_crop by default)
// - Skips or overwrites existing files (configurable)
// - Optionally processes subfolders

// Main entry macro: handles input/output folders, ROI configuration dialog,
// and progress bar initialisation before calling the recursive processor.
macro "Batch Crop Tiffs" {

    // ---- CONFIG (edit if you want defaults) ----
    // Fixed ROI (same as your recorded macro):
    // WHY:
    // - Using a *fixed* rectangle ensures that all cropped images show the same
    //   anatomical region, which is crucial for fair comparisons between cases.
    // - The coordinates here reflect the region chosen during macro recording;
    //   they can be adjusted globally or via the dialog that appears at runtime.
    x = 864;   // left (X coordinate of the ROI origin)
    y = 1296;  // top (Y coordinate of the ROI origin)
    w = 3456;  // width  of ROI in pixels
    h = 2556;  // height of ROI in pixels

    // Suffix added to the base filename (before the .tif extension).
    // WHY: This marks the images as "cropped" versions while keeping
    // the original filenames recognizable.
    suffix = "_crop";          // appended to filename (before .tif)

    // Whether to recurse into subfolders of the input directory.
    // WHY: If your dataset is organised in subfolders (e.g. per case),
    // this allows one-click processing of the whole tree.
    // Set to false to keep behaviour simple and predictable (top-level only).
    processSubfolders = false; // set true to recurse into subfolders

    // Overwrite policy for existing outputs.
    // - overwrite = true  -> previous outputs are replaced
    // - overwrite = false -> existing outputs are preserved, and those files are skipped.
    // WHY: This protects against accidental loss of previous batch results when desired.
    overwrite = true;          // true: overwrite existing outputs; false: skip
    // -------------------------------------------

    // Ask user to choose the INPUT and OUTPUT folders.
    // If either selection is cancelled, the macro exits to avoid partial runs.
    inputDir  = getDirectory("Choose the INPUT folder with .tif/.tiff images");
    if (inputDir == "") exit("No input folder selected.");

    outputDir = getDirectory("Choose the OUTPUT folder");
    if (outputDir == "") exit("No output folder selected.");
    File.makeDirectory(outputDir);

    // Optional dialog to adjust ROI & options at run time.
    // WHY:
    // - Gives flexibility for different experiments while keeping the script
    //   itself stable and version-controlled.
    // - The initial values shown in the dialog come from the CONFIG section above.
    Dialog.create("Crop settings (optional)");
    Dialog.addNumber("X", x);
    Dialog.addNumber("Y", y);
    Dialog.addNumber("Width", w);
    Dialog.addNumber("Height", h);
    Dialog.addString("Output filename suffix", suffix);
    Dialog.addCheckbox("Process subfolders", processSubfolders);
    Dialog.addCheckbox("Overwrite existing outputs", overwrite);
    Dialog.show();
    x = Dialog.getNumber();
    y = Dialog.getNumber();
    w = Dialog.getNumber();
    h = Dialog.getNumber();
    suffix = Dialog.getString();
    processSubfolders = Dialog.getCheckbox();
    overwrite = Dialog.getCheckbox();

    // Count how many .tif/.tiff files will actually be processed.
    // WHY:
    // - This allows us to initialise a meaningful progress bar and gives a clear
    //   message if the folder is empty or contains only non-image files.
    total = countEligible(inputDir);
    if (total == 0) exit("No .tif/.tiff files found in the selected folder.");
    processed = 0;

    // Enable batch mode to speed up processing and avoid constant window redraws.
    setBatchMode(true);
    processDir(inputDir);
    setBatchMode(false);
    showProgress(1,1);
    print("✅ Done. Processed " + processed + " / " + total + " image(s).");
}

// ---------------------------------------------------------------------------
// processDir(dir)
// Recursively (or non-recursively) processes all .tif/.tiff files in `dir`.
// - If processSubfolders is true, subdirectories are explored as well.
// - For each eligible file, calls processFile(path, name).
// ---------------------------------------------------------------------------
function processDir(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        name = list[i];
        path = dir + name;

        // If this entry is a subfolder:
        if (endsWith(name, "/")) {
            // Recurse into subfolder only when the user chose to processSubfolders.
            if (processSubfolders) processDir(path);
            continue;
        }

        // Only process .tif/.tiff images; ignore other file types.
        lname = toLowerCase(name);
        if (endsWith(lname, ".tif") || endsWith(lname, ".tiff")) {
            processFile(path, name);
        }
    }
}

// ---------------------------------------------------------------------------
// processFile(path, name)
// - Opens one image, checks that the fixed ROI fits within its dimensions,
//   performs the crop, and saves a new .tif in the output directory.
// - Updates the global `processed` counter and progress bar.
// ---------------------------------------------------------------------------
function processFile(path, name) {
    open(path);

    // Query image dimensions to ensure that the requested ROI is valid for this file.
    wImg = getWidth();
    hImg = getHeight();

    // Safety check: skip images where the fixed ROI would fall outside
    // the image boundaries (to avoid runtime errors or partial crops).
    if (x < 0 || y < 0 || (x + w) > wImg || (y + h) > hImg) {
        print("⚠️ Skipping (ROI out of bounds): " + name + " [" + wImg + "x" + hImg + "]");
        close();
        return;
    }

    // Apply the fixed ROI and crop the image.
    // At this point, the image window now contains only the selected rectangle.
    makeRectangle(x, y, w, h);
    run("Crop");

    // Build output filename by stripping the original extension and appending
    // the configured suffix (e.g. "_crop") before ".tif".
    base = stripExt(name);
    outPath = outputDir + base + suffix + ".tif";

    // Handle existing output files according to the overwrite policy:
    // - If overwrite==false and the file exists, skip and keep the previous result.
    // - If overwrite==true and the file exists, delete the old version first.
    if (!overwrite && File.exists(outPath)) {
        print("⏭️ Skipping (exists): " + outPath);
        close();
        return;
    }
    if (overwrite && File.exists(outPath)) File.delete(outPath);

    saveAs("Tiff", outPath);
    close();

    // Increment global counter and update progress bar so the user has feedback
    // on how many images have been processed.
    processed++;
    showProgress(processed, total);
}

// Utility: remove the extension from a filename (e.g. "img01.tif" -> "img01").
// Used to preserve the original base name when constructing the output path.
function stripExt(name) {
    last = lastIndexOf(name, ".");
    if (last == -1) return name;
    return substring(name, 0, last);
}

// ---------------------------------------------------------------------------
// countEligible(dir)
// Recursively counts all .tif/.tiff files in `dir` (respecting processSubfolders)
// to determine the total number of images that will be processed.
// This is used only for progress bar initialisation and final reporting.
// ---------------------------------------------------------------------------
function countEligible(dir) {
    list = getFileList(dir);
    n = 0;
    for (i = 0; i < list.length; i++) {
        name = list[i];
        path = dir + name;
        if (endsWith(name, "/")) {
            if (processSubfolders) n += countEligible(path);
        } else {
            lname = toLowerCase(name);
            if (endsWith(lname, ".tif") || endsWith(lname, ".tiff")) n++;
        }
    }
    return n;
}