function main()
    % Main function for the synthetic MRI generation pipeline

    % Load configuration
    cfg = config();

    % Add src subfolders to path
    addpath(genpath('../'));

    %% Define MRI sequence parameter ranges
    fprintf('=== MRI Sequence Parameter Ranges ===\n\n');

    % TE (Echo Time)
    fprintf('TE (Echo Time): %d-%d ms (%d values)\n', min(cfg.params.TE_values), max(cfg.params.TE_values), length(cfg.params.TE_values));
    fprintf('  Clinical significance:\n');
    fprintf('  - Controls T2* weighting and susceptibility sensitivity\n');
    fprintf('  - Shorter TE reduces signal decay from T2* effects\n');
    fprintf('  - Typical range for spoiled gradient echo sequences\n\n');

    % TR (Repetition Time)
    fprintf('TR (Repetition Time): [%s] ms (%d values)\n', num2str(cfg.params.TR_values), length(cfg.params.TR_values));
    fprintf('  Clinical significance:\n');
    fprintf('  - Controls T1 weighting and scan speed\n');
    fprintf('  - Short TR values for fast spoiled gradient echo imaging\n');
    fprintf('  - Allows rapid acquisition with minimal T1 recovery\n\n');

    % FA (Flip Angle)
    fprintf('FA (Flip Angle): %d-%d degrees (%d values)\n', min(cfg.params.FA_values), max(cfg.params.FA_values), length(cfg.params.FA_values));
    fprintf('  Clinical significance:\n');
    fprintf('  - Controls signal intensity and T1 contrast\n');
    fprintf('  - Lower angles (5-15°): Better T1 weighting, reduced SAR\n');
    fprintf('  - Medium angles (20-35°): Balanced contrast and SNR\n');
    fprintf('  - Higher angles (40-60°): Maximum signal but increased SAR\n\n');

    % Calculate total combinations
    total_combinations = length(cfg.params.TE_values) * length(cfg.params.TR_values) * length(cfg.params.FA_values);
    fprintf('Total parameter combinations: %d × %d × %d = %d\n\n', ...
        length(cfg.params.TE_values), length(cfg.params.TR_values), length(cfg.params.FA_values), total_combinations);

    %% Run the synthetic image generation
    fprintf('Starting synthetic image generation...\n');
    fprintf('Patient directory: %s\n', cfg.paths.patient_dir);
    fprintf('Output directory: %s\n\n', cfg.paths.output_dir);

    % Record start time
    start_time = tic;

    % Generate synthetic images with the specified parameters
    generate_synthetic_images(cfg);

    % Generate delta B0 maps
    runRomeo(cfg);

    % Convert images to k-space
    nifti_to_kspace(cfg);

    % Calculate and display processing time
    elapsed_time = toc(start_time);
    fprintf('\n=== Processing Complete ===\n');
    fprintf('Total processing time: %.2f minutes\n', elapsed_time/60);
    fprintf('Average time per image: %.2f seconds\n', elapsed_time/total_combinations);
    fprintf('Generated %d synthetic MRI images\n', total_combinations);
    fprintf('Output location: %s\n', cfg.paths.output_dir);

    %% Summary of generated files
    fprintf('\n=== Generated Files Summary ===\n');
    fprintf('Each parameter combination generates:\n');
    fprintf('  1. NIfTI image file (.nii)\n');
    fprintf('  2. JSON metadata file (.json)\n');
    fprintf('  3. NIfTI kspace file- both magnitude and phase-(.nii)\n');
    fprintf('  4. JSON metadata for kspace file (.json)\n');
    fprintf('Total files created: %d\n', total_combinations * 4);

    fprintf('\nFile naming convention:\n');
    fprintf('  Images: [PatientID]_synthim_TE[X]_FA[Y]_TR[Z].nii\n');
    fprintf('  Metadata: [PatientID]_synthim_TE[X]_FA[Y]_TR[Z].json\n');

    fprintf('\nScript execution completed successfully!\n');
end

