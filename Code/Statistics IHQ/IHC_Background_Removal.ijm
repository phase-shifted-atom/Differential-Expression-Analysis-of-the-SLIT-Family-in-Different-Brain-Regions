// ==== Pale low-stain tissue + keep background pure white ====
// - The goal is to:
//   * force the background to be pure white (S=0, B=255 in HSB),
//   * make low-stain tissue/background more pale by **reducing** saturation and increasing brightness,
//   * thereby suppressing non-specific antibody binding that appears as diffuse, high-intensity background.
// ---- TUNABLES ----
hMinBG=0;   hMaxBG=255;  sMinBG=0;   sMaxBG=30;   bMinBG=220; bMaxBG=255;  // BG selection
hMinLS=0;   hMaxLS=255;  sMinLS=30;   sMaxLS=180;  bMinLS=90; bMaxLS=220;  // LOW_S selection
shrinkPx = 0;            // shrink both selections (0,1,2…)
// Intensity adjustment parameters for LOW_S tissue/background:
// - lowS_S_Mult: multiplies saturation within LOW_S selection.
//   In this macro, the intended use is with a value **< 1** (default 0.03),
//   which strongly reduces saturation and makes low-stain regions more pale.
//   Using values >1 would instead enhance colour in LOW_S, which is *not*
//   the target behaviour here.
// - lowS_B_Add: constant brightness added to LOW_S region.
//   Together with saturation reduction, this pushes low-stain areas towards
//   pure white, helping to remove non-specific antibody background.
lowS_S_Mult = 0.03;      // Saturation multiplier on LOW_S (try 0.30 for a visible test)
lowS_B_Add  = 40;        // Brightness add on LOW_S
// Debug
logCoverage = true;      // print coverage % of masks
logMeans    = true;      // print S/B means before/after in LOW_S
// -------------------

function processPaleLowStain() {
  if (bitDepth!=24) run("RGB Color");

  // Convert to HSB (1=Hue, 2=Saturation, 3=Brightness)
  run("HSB Stack");
  hsbTitle = getTitle();
  roiManager("Reset"); // ensure ROI Manager exists and is empty

  // --- helper: build SB mask on HSB, add ROI to manager, return its index and coverage ---
  function addSBROI_toHSB_withCov(sMin,sMax,bMin,bMax, titleHSB, roiName) {
    // Build S mask
    selectWindow(titleHSB); setSlice(2);
    run("Duplicate...", "title=__S");
    setThreshold(sMin, sMax);
    run("Convert to Mask");

    // Build B mask
    selectWindow(titleHSB); setSlice(3);
    run("Duplicate...", "title=__B");
    setThreshold(bMin, bMax);
    run("Convert to Mask");

    // AND -> combined
    run("Image Calculator...", "image1=__S image2=__B operation=AND create");
    selectWindow("Result of __S"); rename("__SB");

    // Coverage % on binary mask
    run("Select All");
    getStatistics(area, mean);
    cov = (mean/255.0)*100.0;
    if (logCoverage) print(roiName+" coverage: "+d2s(cov,2)+"% (S["+sMin+".."+sMax+"], B["+bMin+".."+bMax+"]) ");

    // Create selection on the mask window
    run("Create Selection");
    if (selectionType==-1) {
      // clean temp windows and signal failure
      selectWindow("__S"); close();
      selectWindow("__B"); close();
      selectWindow("__SB"); close();
      return -1;
    }

    // Add to ROI Manager and rename safely
    roiManager("Add");
    idx = roiManager("count")-1;
    roiManager("Select", idx);
    roiManager("Rename", roiName);

    // Transfer selection to the HSB image (so it exists on that window)
    selectWindow(titleHSB);
    roiManager("Select", idx);

    // Clean temps
    selectWindow("__S"); close();
    selectWindow("__B"); close();
    selectWindow("__SB"); close();

    return idx;
  }

  // --- BACKGROUND: build ROI, apply S=0 and B=255 ---
  idxBG = addSBROI_toHSB_withCov(sMinBG, sMaxBG, bMinBG, bMaxBG, hsbTitle, "BG");
  if (idxBG>=0) {
    selectWindow(hsbTitle);
    roiManager("Select", idxBG);
    if (shrinkPx>0) run("Enlarge...", "enlarge=-"+shrinkPx+" pixel");
    setSlice(2); run("Set...", "value=0 slice");     // Saturation = 0
    setSlice(3); run("Set...", "value=255 slice");   // Brightness = 255
    run("Select None");
  }

  // --- LOW_S: REBUILD ROI AFTER BG EDITS ---
  idxLS = addSBROI_toHSB_withCov(sMinLS, sMaxLS, bMinLS, bMaxLS, hsbTitle, "LOW_S");
  if (idxLS>=0) {
    selectWindow(hsbTitle);
    roiManager("Select", idxLS);
    if (shrinkPx>0) run("Enlarge...", "enlarge=-"+shrinkPx+" pixel");

    // Measure S mean before
    if (logMeans) {
      setSlice(2);
      getStatistics(aSpre, mSpre);
    }

    // WHY: This globally scales saturation inside LOW_S. With the default factor
    // <1, saturation is strongly reduced so that low-stain regions become less
    // coloured and more similar to white background, suppressing non-specific DAB.
    // Apply LOW_S S multiply (forced to current slice)
    setSlice(2); run("Multiply...", "value="+lowS_S_Mult+" slice");

    // Measure S mean after
    if (logMeans) {
      setSlice(2);
      getStatistics(aSpost, mSpost);
      // Compute factor safely (avoid max())
      factor = 0;
      if (mSpre>0) factor = mSpost/mSpre;
      print("LOW_S Saturation mean: before="+d2s(mSpre,2)+" after="+d2s(mSpost,2)+" (x"+d2s(factor,2)+")");
    }

    // Measure B mean before
    if (logMeans) {
      setSlice(3);
      getStatistics(aBpre, mBpre);
    }

    // WHY: Adding a constant brightness makes low-stain regions brighter and,
    // combined with saturation reduction, moves them closer to pure white,
    // effectively removing residual non-specific background signal.
    // Apply LOW_S B add (forced to current slice)
    setSlice(3); run("Add...", "value="+lowS_B_Add+" slice");

    // Measure B mean after
    if (logMeans) {
      setSlice(3);
      getStatistics(aBpost, mBpost);
      print("LOW_S Brightness mean: before="+d2s(mBpre,2)+" after="+d2s(mBpost,2)+" (+"+d2s(mBpost-mBpre,2)+")");
    }

    run("Select None");
  } else {
    if (logCoverage) print("LOW_S selection empty — adjust sMaxLS or widen brightness band.");
  }

  // Recombine to RGB
  run("RGB Color");
}

macro "Pale low-stain + white BG (HSB, verified)" {
  if (nImages==0) {
    input = getDirectory("Choose input folder");
    if (input=="") exit("No input folder chosen.");
    output = getDirectory("Choose output folder");
    if (output=="") exit("No output folder chosen.");
    File.makeDirectory(output);

    list = getFileList(input);
    setBatchMode(true);
    for (i=0; i<list.length; i++) {
      name = list[i];
      if (endsWith(name, "/")) continue; // skip subfolders
      lower = toLowerCase(name);
      if (!(endsWith(lower, ".tif")||endsWith(lower, ".tiff")||endsWith(lower, ".jpg")||endsWith(lower, ".jpeg")||endsWith(lower, ".png")||endsWith(lower, ".bmp"))) continue;

      open(input+name);
      print("\n=== Processing: "+name+" ===");
      processPaleLowStain();

      if (endsWith(lower, ".tif")||endsWith(lower, ".tiff")) {
        saveAs("Tiff", output+name);
      } else if (endsWith(lower, ".jpg")||endsWith(lower, ".jpeg")) {
        saveAs("Jpeg", output+name);
      } else if (endsWith(lower, ".png")) {
        saveAs("PNG", output+name);
      } else if (endsWith(lower, ".bmp")) {
        saveAs("BMP", output+name);
      } else {
        base = File.nameWithoutExtension(name);
        saveAs("Tiff", output+base+".tif");
      }

      close();
      roiManager("Reset");
    }
    setBatchMode(false);
    print("\nDone. Output folder: "+output);
  } else {
    setBatchMode(true);
    processPaleLowStain();
    setBatchMode(false);
  }
}

// HOW:
// 1. Convert RGB image to HSB stack (3 slices: Hue, Saturation, Brightness).
// 2. Select background (BG) by thresholding saturation and brightness slices.
// 3. Set BG saturation to 0 and brightness to 255 (pure white).
// 4. For LOW_S mask: multiply saturation by a factor <1 and add brightness locally
//    so that low-stain tissue/background becomes less saturated and brighter (more white),
//    reducing non-specific DAB haze rather than enhancing weak signal.
// 5. Convert back to RGB.