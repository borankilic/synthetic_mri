# mpm_to_weighted_codes (Refactored)

This project provides a comprehensive and refactored pipeline for generating synthetic weighted MRI images from quantitative multi-parameter maps (R1, R2s, PD). It allows for flexible adjustment of sequence parameters (TE, TR, FA) to simulate various image contrasts. The pipeline also includes functionality to calculate the corresponding k-space data for each coil in an MRI scanner by estimating sensitivity maps.

## Project Overview

The core of this project is to take quantitative MRI maps (R1, R2s, PD) and generate realistic synthetic weighted images that would be produced by an MRI scanner with specific sequence parameters. This is achieved by applying the Ernst equation, which models the MRI signal intensity based on tissue properties and sequence parameters.

## New Folder Structure

```
/mpm_to_weighted_codes_refactored/
|
├── data/                     # Input and output data
│   ├── input/                # For raw input data (MPMs, B1 maps, etc.)
│   └── output/               # For generated data (synthetic images, k-space, etc.)
|
├── src/                      # All MATLAB source code
│   ├── main/                 # Main script and configuration
│   │   ├── main.m
│   │   └── config.m
│   │
│   ├── preprocessing/        # Preprocessing steps
│   │   ├── extract_brain.m
│   │   └── register_maps.m
│   │
│   ├── generation/           # Core data generation
│   │   ├── generate_synthetic_images.m
│   │   ├── compute_ernst_signal.m
│   │   ├── create_smaps.m
│   │   └── nifti_to_kspace.m
│   │
│   ├── utils/                # Utility functions
│   │   ├── file_io/
│   │   └── image_processing/
│   │
│   └── visualization/        # Visualization scripts
|
├── scripts/                  # Shell scripts and external tools
│   ├── run_romeo.sh
│   └── romeo.jl
|
└── README.md
```

## How to Run

1.  **Set up the environment:**
    *   Make sure you have MATLAB installed with the Image Processing Toolbox.
    *   Install FSL and make sure it's in your system's PATH.
    *   Install the BART toolbox and set the `BART_TOOLBOX_PATH` environment variable.
    *   Install Julia and the required packages for ROMEO (`ROMEO`, `MriResearchTools`, `ArgParse`).
2.  **Configure the parameters:**
    *   Open the `src/main/config.m` script and modify the paths and parameters as needed.
3.  **Run the script:**
    *   Run the `src/main/main.m` script in MATLAB.

The script will then generate the synthetic images, k-space data, and B0 maps in the specified output directory.

## Dependencies

-   MATLAB (with Image Processing Toolbox)
-   FSL
-   BART
-   Julia (with `ROMEO`, `MriResearchTools`, `ArgParse` packages)
-   SPM12 (for DARTEL)
-   hMRI toolbox (for DARTEL)