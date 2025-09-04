# mpm_to_weighted_codes

This project provides a comprehensive pipeline for generating synthetic weighted MRI images from quantitative multi-parameter maps (R1, R2s, PD). It allows for flexible adjustment of sequence parameters (TE, TR, FA) to simulate various image contrasts. The pipeline also includes functionality to calculate the corresponding k-space data for each coil in an MRI scanner by estimating sensitivity maps.


## Workflow
The core of this project is to take quantitative MRI maps (R1, R2s, PD) and generate realistic synthetic weighted images that would be produced by an MRI scanner with specific sequence parameters. This is achieved by applying the Ernst equation, which models the MRI signal intensity based on tissue properties and sequence parameters.
The main workflow is controlled by the `main.m` script. Here's a high-level overview of the steps:

1.  **Setup Paths and Parameters:** The `main.m` script starts by setting up the paths to the input data and defining the ranges for the sequence parameters (TE, TR, FA).
2.  **Generate Synthetic Images:** The `generate_synthetic_images.m` function is called to generate the synthetic weighted images. This function performs the following sub-steps:
    *   **Brain Extraction:** If specified, it extracts the brain from the quantitative maps using either `extract_brain_bet.m` or `extract_brain_fsl.m`.
    *   **B1 Map Registration:** It registers the B1 map to the quantitative maps using `register_b1_to_ref.m`
    *   **Signal Calculation:** It calculates the synthetic signal for each combination of TE, TR, and FA using `compute_ernst_signal.m`.
    *   **Save Images:** It saves the generated synthetic images as NIfTI files using `save_nifti_image.m` and writes the corresponding metadata to JSON files using `write_json_metadata.m`.
3.  **Generate B0 Maps:** The `runRomeo.m` function is called to generate B0 maps from multi-echo phase data.
4.  **Generate K-space Data:** The `nifti_to_kspace.m` function is called to convert the synthetic images to k-space data.
    *   **Sensitivity Map Generation:** If not already present, it generates sensitivity maps using `create_smaps.m`.
    *   **Phase evolution:** Phases of synthetic images  are assumed to be caused by phase portion of sensitivity maps and the B0 field inhomogenities under the assumption that phase evolves linearly with echo time.
    *   **Calculation of per-coil images:** Images are calculated for each coil seperately using respective sensitvity maps and stacked.
    *   **Taking FFT and Saving the images:** 3D FFT of each per-coil image is calculated and resulting k-space data is saved as seperate magnitude and phase NIFTI files, each with dimensions `(n_x, n_y, n_z, n_coil)`.
    
## How to Run

1.  **Set up the environment:**
    *   Make sure you have MATLAB installed with the Image Processing Toolbox.
    *   Install FSL, make sure it's in your system's PATH and it's accesible from the terminal via `FSL` command.
    *   Install Julia and the required packages for ROMEO (`ROMEO`, `MriResearchTools`, `ArgParse`).
2.  **Configure the paths:**
    *   Open the `main.m` script and modify the paths in the "Setup paths" section to point to your data.
    *   We adopt the folder structure of hMRI-toolbox for input directory. Your patient directory (`patient_dir`) should be structured as below:
    ```
    /patient_dir/
    |
    ├── Results/                     # MPMs (there should be a corresponding json file for each nifti, containing essential metadata)
    │   ├── <R1_filename>.nii           # Must hace 'R1' in its name.
    │   └── <R2_filename>.nii           # Must have 'R2' in its name.
    │   └── <PD_filename>.nii           # Must have 'PD' in its name.
    │
    │
    ├── B1mapCalc/                   # B1+ (RF transmit sensitivity) maps
    │   ├── <B1map_filename>.nii        # Must have 'B1map' in its name.
    │   └── <B1red_filename>.nii        # Anatomical reference file for coregisteration. mUst have 'B1ref' in its name.
    │   
    ├── MPMCalc/                     # FOlder containing Tissue masks (GM, WM, CSF)
        ├── <gm_mask_filename>.m        # Must hace 'c1' in its name
        └── <wm_mask_filename>.m        # Must hace 'c2' in its name
        └── <csf_mask_filename>.m       # Must hace 'c3' in its name

    ```
    NOTE: If your input files are named differently, you can change file name patterns from [create_data_struct.m](./src/utils/file_io/create_data_struct.m)
    * `kspace_path` should point to raw kspace. Pipeline is configured Siemens` twix file format with extension `.dat`. It is read and mapped by [mapVBVD.m](./src/scripts/mapVBVD.m). If you have the kspace data in a differetn format, consider editing the part in [create_smaps.m](./src/generation/create_smaps.m] where the data is imported. Make sure kspace is data is 3D and contains the ACS region even thoguh it is undersampled.
    * `im_mag_dir` should point to the directory where weighted (T1,PD etc.) magnitude images (which MPMs are created from) are stored.  Directory should contain each echo in a seperate NIFTI file. Make sure length of the list `real_echo_times` provided and number of NIFTI files in `im_mag_dir` match.
    * `im_phase_dir` should point to the directory where weighted (T1,PD etc.) phase images (which MPMs are created from) are stored.  Directory should contain each echo in a seperate NIFTI file. Make sure length of the list `real_echo_times` provided and number of NIFTI files in `im_phase_dir` match.
    * `espririt_path` should point to the directory where the MATLAB library Matlab Library including implementations of [SPIRiT, ESPIRiT, Coil Compression, SAKE](https://people.eecs.berkeley.edu/~mlustig/Software.html). You can download the library from [here](https://people.eecs.berkeley.edu/~mlustig/software/SPIRiT_v0.3.tar.gz).
    * `romeo_script_path` should point to your `romeo.jl file`. Please follow the steps in [here ](https://github.com/korbinian90/ROMEO.jl#usage---command-line) to setup ROMEO.
    * `real_echo_times` be a list of echo times in ms. You can find the specific echo time in your sequence configuration file.
    
3.  **Run the script:**
    *   Run the `main.m` script in MATLAB.

The script will then generate the synthetic images, k-space data, and B0 maps in the specified output directory.

## Function Descriptions

Here's a description of the key functions in the project:

### `main.m`

This is the main script that drives the entire workflow. It sets up the paths, defines the sequence parameters, and calls the other functions to generate the synthetic data.

**Key Parameters:**

-   `patient_dir`: The directory containing the patient's quantitative maps.
-   `kspace_path`: The path to the raw k-space data file.
-   `im_mag_dir`: The directory containing the magnitude images for B0 map generation.
-   `im_phase_dir`: The directory containing the phase images for B0 map generation.
-   `output_dir`: The directory where the generated data will be saved.
-   `TE_values`: A vector of echo times (in ms).
-   `TR_values`: A vector of repetition times (in ms).
-   `FA_values`: A vector of flip angles (in degrees).

### `generate_synthetic_images.m`

This function generates the synthetic weighted images.

**Key Parameters:**

-   `patient_dir`: The directory containing the patient's quantitative maps.
-   `output_dir`: The directory where the generated images will be saved.
-   `kspace_path`: The path to the raw k-space data file.
-   `TE_vals`: A vector of echo times (in ms).
-   `TR_vals`: A vector of repetition times (in ms).
-   `alpha_vals`: A vector of flip angles (in degrees).
-   `signal_constant`: A scaling factor for the signal intensity.
-   `interp_method`: The interpolation method to use for B1 map registration.
-   `extract_brain`: The brain extraction method to use ('bet', 'fsl', or 'no').

### `compute_ernst_signal.m`

This function computes the MRI signal using the Ernst equation.

**Inputs:**

-   `data_struct`: A struct containing the PD, R1, R2s, and B1 maps.
-   `TE`: The echo time (in ms).
-   `TR`: The repetition time (in ms).
-   `alpha`: The flip angle (in degrees).
-   `signal_constant`: A scaling factor for the signal intensity.
-   `interp_method`: The interpolation method to use for B1 map resizing.
-   `extract_brain`: A string indicating whether to use the brain-extracted data.

**Output:**

-   `signal`: The calculated MRI signal.

### `create_smaps.m`

This function computes ESPIRiT sensitivity maps from raw k-space data.

**Inputs:**

-   `twix_filename`: The path to the Siemens TWIX file.
-   `sens_maps_out_file`: The output filename for the sensitivity maps.

### `nifti_to_kspace.m`

This function converts NIfTI images to k-space data using an inverse FFT.

**Inputs:**

-   `input_dir`: The directory containing the NIfTI files.

### `runRomeo.m`

This function runs the ROMEO algorithm to generate B0 maps from multi-echo phase data.

**Inputs:**

-   `magFolder`: The folder containing the magnitude images.
-   `phaseFolder`: The folder containing the phase images.
-   `echoTimes`: A vector of echo times.
-   `output_path`: The output path for the B0 map.
-   `mpm_dir`: The directory containing the MPMs.
-   `romeoScriptPath`: The path to the `romeo.jl` script.

### Other Helper Functions

The project also includes several helper functions for tasks such as:

-   `cfl_to_nifti.m`: Converts CFL files to NIfTI format.
-   `create_data_struct.m`: Creates a struct to hold the paths to the input data.
-   `DARTEL.m`: Registers B1 maps using DARTEL.
-   `extract_brain_bet.m`: Extracts the brain using FSL's BET.
-   `extract_brain_fsl.m`: Extracts the brain using FSL and tissue probability maps.
-   `fixHotPixels.m`: Fixes hot pixels in an image.
-   `get_bart_path2.m`: Gets the path to the BART toolbox.
-   `kspace_to_nifti.m`: Converts k-space data to NIfTI format.
-   `load_mri_data.m`: Loads MRI data from NIfTI files.
-   `nifti_to_cfl.m`: Converts NIfTI files to CFL format.
-   `register_b1_to_ref.m`: Registers B1 maps to a reference image.
-   `resize_b1_map.m`: Resizes B1 maps.
-   `save_nifti_image.m`: Saves NIfTI images.
-   `twix2cfl.m`: Converts TWIX files to CFL format.
-   `validate_mri_data.m`: Validates the input MRI data.
-   `view_brain_slices.m`: Visualizes brain slices.
-   `visualize_ESPIRIT.m`: Visualizes ESPIRiT sensitivity maps.
-   `visualize_smap.m`: Visualizes sensitivity maps.
-   `write_json_metadata.m`: Writes metadata to JSON files.

## How to Run

1.  **Set up the environment:**
    *   Make sure you have MATLAB installed with the Image Processing Toolbox.
    *   Install FSL and make sure it's in your system's PATH.
    *   Install the BART toolbox and set the `BART_TOOLBOX_PATH` environment variable.
    *   Install Julia and the required packages for ROMEO (`ROMEO`, `MriResearchTools`, `ArgParse`).
2.  **Configure the paths:**
    *   Open the `main.m` script and modify the paths in the "Setup paths" section to point to your data.
3.  **Run the script:**
    *   Run the `main.m` script in MATLAB.

The script will then generate the synthetic images, k-space data, and B0 maps in the specified output directory.

## Dependencies

-   MATLAB (with Image Processing Toolbox)
-   FSL
-   BART
-   Julia (with `ROMEO`, `MriResearchTools`, `ArgParse` packages)

