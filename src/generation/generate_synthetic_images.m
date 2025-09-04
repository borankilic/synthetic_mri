function generate_synthetic_images(cfg)
    % GENERATE_SYNTHETIC_IMAGES Generate synthetic MRI images using Ernst equation
    %   generate_synthetic_images(cfg) generates synthetic images using parameters
    %   from the cfg struct.

    % Extract parameters from config struct
    patient_dir = cfg.paths.patient_dir;
    output_dir = cfg.paths.output_dir;
    kspace_path = cfg.paths.kspace_path;
    TE_vals = cfg.params.TE_values;
    TR_vals = cfg.params.TR_values;
    alpha_vals = cfg.params.FA_values;
    signal_constant = cfg.options.signal_constant;
    interp_method = cfg.options.interp_method;
    extract_brain_method = cfg.options.extract_brain;

    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Create file structure
    data_struct = create_data_struct(patient_dir);

    % Apply brain mask to quantitative maps
    extracted_mpms_temp = fullfile(output_dir, 'extracted_mpms/');
    if ~exist(extracted_mpms_temp, 'dir'), mkdir(extracted_mpms_temp), end

    if ~strcmp(extract_brain_method, 'no')
        extract_brain(data_struct, extracted_mpms_temp, 'method', extract_brain_method, 'patient_dir', patient_dir);
    end

    % Create sensitivity maps if doesnt exist yet
    sens_map_dir = fullfile(output_dir, 'sens_maps');
    smap_file_path = fullfile(sens_map_dir, replace(data_struct.R1.name,'R1.nii', 'smap.nii.gz'));
    if ~exist(sens_map_dir, 'dir'), mkdir(sens_map_dir), end

    if ~exist(replace(smap_file_path, 'smap', 'smap_mag'), 'file')
        fprintf("Creating sensitivity maps using ESPIRiT algorithm.")
        create_smaps(kspace_path, smap_file_path, ...
        'espirit_path'  , '/data/u_kilicb_software/ESPIRiT/', ...
        'ncalib'      , 20, ...
        'ksize'      , [5 5], ...
        'eigThresh1' , 0.1, ...
        'eigThresh2' , 0.5, ...
        'show_figures' ,  'none' , ....   % options: 'none','essential','detailed'
        'verbose' , true, ...
        'data_struct', data_struct, ...
        'smoothing_radius', 25);          % Smoothing radius of gaussian kernel in mm
    
    else
        fprintf("Skipping smap creation. Sensitivity maps already exists.\n");
    end

    % Register the B1 map and smap to ref maps
    registered_B1_temp = fullfile(output_dir, 'registered_B1/');
    if ~exist(registered_B1_temp, 'dir'), mkdir(registered_B1_temp), end

    register_b1_to_ref(data_struct.B1map.filepath, data_struct.B1ref.filepath, data_struct.R1.filepath, registered_B1_temp, interp_method);

    data_struct = create_data_struct(output_dir, data_struct);

    % Load MRI data
    fprintf('Loading MRI data from: %s\n', output_dir);
    data_struct = load_mri_data(data_struct);
    patient_id = erase(data_struct.R1.name, '_R1.nii');
    if contains(patient_id, '.nii'), erase(patient_dir, '.nii'), end

    % Generate synthetic images for all parameter combinations
    fprintf('Generating synthetic images...\n');
    total_combinations = length(TE_vals) * length(TR_vals) * length(alpha_vals);
    combination_count = 0;

    for TE = TE_vals
        for TR = TR_vals
            for alpha = alpha_vals
                combination_count = combination_count + 1;

                fprintf('Processing combination %d/%d: TE=%d, TR=%d, FA=%d\n', ...
                    combination_count, total_combinations, TE, TR, alpha);

                % Generate synthetic signal
                signal = compute_ernst_signal(data_struct, TE, TR, alpha, ...
                    signal_constant, interp_method, extract_brain_method);

                %Fix hot pixels (if any)
                threshold_std_dev = 1000; % Threshold for exluding pixels from the image in units of tsndart deviation
                thresholding_method = 'mean';
                signal = fixHotPixels(signal, threshold_std_dev, thresholding_method);
                %Normalize the signal range for numerical stability
                mean_val = mean(signal(:));
                signal = (signal ./ mean_val);

                % Save result
                output_filename = sprintf('%s_synthim_TE%d_FA%d_TR%d.nii', ...
                    patient_id, TE, alpha, TR);
                save_dir = fullfile(output_dir, 'synthetic_images');
                if ~exist(save_dir, 'dir'), mkdir(save_dir), end
                output_path = fullfile(save_dir, output_filename);
                
                save_nifti_image(signal, output_path, data_struct.R1.ref_header);
                write_json_metadata(output_path, data_struct.R1.ref_header, TE, TR, alpha);

                fprintf('Saved: %s\n', output_filename);
            end
        end
    end

    fprintf('Finished generating %d synthetic images.\n', total_combinations);
end
