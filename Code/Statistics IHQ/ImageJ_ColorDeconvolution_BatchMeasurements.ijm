// ============================================================================
// Macro: ImageJ_ColorDeconvolution_BatchMeasurements.ijm
// ----------------------------------------------------------------------------
// PURPOSE:
// - Automate colour deconvolution and measurement for a *folder* of IHC images.
// - Extract quantitative measurements only from the stain of interest
//   (here, Colour_2 = typically the DAB channel in the [H DAB] vector),
//   while discarding the other channels and the original RGB image.
// - Append measurements to the ImageJ "Results" table so they can be exported
//   as a tabular file (e.g. for downstream R-based statistics).
//
// CONTEXT IN THE PIPELINE:
// - This macro assumes that images have already been pre-processed
//   (e.g. background harmonisation, cropping to a standard ROI).
// - It focuses on extracting stain-specific intensity/area information that
//   will be summarised per ROI/case in the IHC statistics scripts.
//
// HOW (HIGH-LEVEL):
// 1. Ask the user to select an input folder with IHC images.
// 2. Loop over each image file and run "Colour Deconvolution2" with [H DAB].
// 3. Close the original RGB image + non-target channels (Colour_1, Colour_3).
// 4. Measure only Colour_2 (DAB) for each image and close it.
// 5. Leave all measurements accumulated in the "Results" table.
//
// NOTE:
// - This macro does NOT save new image files; it only populates the Results
//   table. Export from "Results" afterwards (e.g. as .csv or .tsv).
// ============================================================================

// -----------------------------
// Macro: ProcessAllImages_FromFolder.ijm
// -----------------------------

// 1) Helper: index of a value in an array
// WHY:
// - ImageJ macro language does not have a built-in function to get the index
//   of a string in an array (e.g. the list of open image titles).
// - This helper lets us safely check if a particular image window exists
//   before attempting to select/close it, which avoids errors when expected
//   windows are missing (e.g. if Colour_3 is not created).
// HOW:
// - Returns the index (0-based) of the first match, or -1 if the value is
//   not found.
function findIndexOf(array, value) {
    for (jj = 0; jj < array.length; jj++) {
        if (array[jj] == value) return jj;
    }
    return -1; // not found
}

// 2) Ask the user to choose the folder containing the images to process
// WHY:
// - This macro is designed for batch processing; we never operate on a single,
//   manually opened image here.
// - If the user cancels or provides an empty path, we stop early to avoid
//   partial or unintended processing.
inputDir = getDirectory("Choose the input folder with images");
tmp = "" + inputDir; // coerce to string; becomes "null" if user cancels
if (tmp=="null" || tmp=="") exit("No folder selected.");

// Retrieve a list of all entries (files and subfolders) in the chosen folder
// and then filter out subfolders during the loop.
fileList = getFileList(inputDir);

// (Optional) Speed up batch processing by disabling UI updates.
// WHY:
// - setBatchMode(true) avoids repainting image windows on screen for every
//   operation, which significantly accelerates large batch runs and reduces
//   flicker. It does not change the underlying computations.
setBatchMode(true);

// 3) Loop over all entries in the folder
for (i = 0; i < fileList.length; i++) {
    name = fileList[i];
    fullPath = inputDir + name;

    // Skip any subfolders: this macro processes only files in the selected
    // directory (no recursion) to keep behaviour predictable.
    if (File.isDirectory(fullPath)) continue;

    // Process only common image extensions.
    // WHY:
    // - This prevents accidental attempts to "open" non-image files and ensures
    //   that only expected formats are processed (TIF/TIFF/PNG/JPEG/BMP).
    lower = toLowerCase(name);
    if (!(endsWith(lower, ".tif") || endsWith(lower, ".tiff") ||
          endsWith(lower, ".png") || endsWith(lower, ".jpg")  ||
          endsWith(lower, ".jpeg")|| endsWith(lower, ".bmp"))) {
        continue;
    }

    // Open the current image file.
    // NOTE: "originalTitle" below will capture the exact window title that ImageJ
    // assigns, which is then used to construct the expected Colour_* titles.
    open(fullPath);
    originalTitle = getTitle();  // use the actual window title as produced by ImageJ

    // Run Colour Deconvolution2 using the [H DAB] stain vector.
    // WHY:
    // - Colour_1 (H) corresponds to the hematoxylin channel,
    // - Colour_2 (DAB) corresponds to the brown chromogen (target IHC signal),
    // - Colour_3 is the residual.
    // We are interested in quantitative measurements from Colour_2 only.
    run("Colour Deconvolution2", "vectors=[H DAB] output=8bit_Transmittance simulated cross hide");

    // Construct the expected titles for the deconvolved images.
    // NOTE:
    // - Colour Deconvolution2 appends suffixes like "-(Colour_1)" etc. to the
    //   original image title. We recreate those names to find/select the correct
    //   windows later.
    titleCol1 = originalTitle + "-(Colour_1)";
    titleCol2 = originalTitle + "-(Colour_2)";
    titleCol3 = originalTitle + "-(Colour_3)";

    // Close (Colour_3) if it exists.
    // WHY:
    // - Colour_3 is typically the residual channel and is not needed for
    //   quantitative analysis here. Closing it saves memory and avoids clutter.
    if (findIndexOf(getList("image.titles"), titleCol3) != -1) {
        selectImage(titleCol3);
        close();
    }

    // Close (Colour_1) if it exists.
    // WHY:
    // - Colour_1 corresponds to the hematoxylin channel. In this project, we are
    //   not quantifying hematoxylin, so we discard it after deconvolution.
    if (findIndexOf(getList("image.titles"), titleCol1) != -1) {
        selectImage(titleCol1);
        close();
    }

    // Close the original RGB image if it is still open.
    // WHY:
    // - After colour deconvolution, the original is no longer needed for this
    //   macro and keeping it would consume extra memory in large batches.
    if (findIndexOf(getList("image.titles"), originalTitle) != -1) {
        selectImage(originalTitle);
        close();
    }

    // Select Colour_2 (DAB channel) -> measure -> close.
    // WHY:
    // - DAB intensity and area are the main quantitative outputs used for the
    //   downstream IHC statistical analyses.
    // - Measurements are appended to the global "Results" table for all images.
    // - After measurement, we close Colour_2 to keep memory use bounded.
    if (findIndexOf(getList("image.titles"), titleCol2) != -1) {
        selectImage(titleCol2);
        run("Measure");
        close();
    }
}

// End batch mode: re-enable normal UI updates now that all images are done.
setBatchMode(false);

// (Optional) Bring the Results window to the front at the end of the run.
// Uncomment the line below if you want the measurements table to pop up
// automatically after processing.
// selectWindow("Results");