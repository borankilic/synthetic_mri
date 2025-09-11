function data_struct  = generate_synthetic_images(cfg)
    % GENERATE_SYNTHETIC_IMAGES Generate synthetic MRI images using Ernst equation
    %   generate_synthetic_images(cfg) generates synthetic images using parameters
    %   from the cfg struct.

    % Extract parameters from config struct
    patient_dir = cfg.paths.patient_dir;
    output_dir = cfg.paths.output_dir;
    TE_vals = cfg.params.TE_values;
    TR_vals = cfg.params.TR_values;
    TI_vals = cfg.params.TI_values;
    alpha_vals = cfg.params.FA_values;
    signal_constant = cfg.options.signal_constant;
    interp_method = cfg.options.interp_method;
    extract_brain_method = cfg.options.extract_brain;
    sequence = cfg.params.sequence;
    threshold_tmap = cfg.options.threshold_tmap;
    threshold_bet = cfg.options.threshold_bet;

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
        extract_brain(data_struct, extracted_mpms_temp, 'method', extract_brain_method, 'patient_dir', patient_dir,'threshold_tmap', threshold_tmap, 'threshold_bet', threshold_bet);
    end



    % Register the B1 map and smap to ref maps
    registered_B1_temp = fullfile(output_dir, 'registered_B1/');
    if ~exist(registered_B1_temp, 'dir'), mkdir(registered_B1_temp), end

    register_b1_to_ref(data_struct, data_struct.R1.filepath, registered_B1_temp, interp_method);

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
            for TI = TI_vals
                for alpha = alpha_vals
                    combination_count = combination_count + 1;
    
                    fprintf('Processing combination %d/%d: TE=%d, TR=%d, FA=%d\n', ...
                        combination_count, total_combinations, TE, TR, alpha);
    
                    % Generate synthetic signal
                    signal = compute_ernst_signal(data_struct, TE, TR, TI, alpha, ...
                        signal_constant, interp_method, extract_brain_method, sequence);
                    
                    % Add RF_bias
                    if cfg.options.add_RF_bias && isfield(data_struct, 'RF_bias_map')
                        fprintf("Applying RF bias to the image")
                        signal = signal .* data_struct.RF_bias_map;
                    end
                    
                    %Fix hot pixels (if any)
                    threshold_std_dev = 1000; % Threshold for exluding pixels from the image in units of tsndart deviation
                    thresholding_method = 'mean';
                    signal = fixHotPixels(signal, threshold_std_dev, thresholding_method);
                    

                    %% Save result
                    output_filename = generate_filename(patient_id, sequence, TE, alpha, TR, TI);
                    save_dir = fullfile(output_dir, 'synthetic_images');
                    if ~exist(save_dir, 'dir'), mkdir(save_dir), end
                    output_path = fullfile(save_dir, output_filename);
                    
                    save_nifti_image(signal, output_path, data_struct.R1.ref_header);
                    write_json_metadata(output_path, data_struct.R1.ref_header, TE, TR, alpha, cfg);
    
                    fprintf('Saved: %s\n', output_filename);
                end
            end
        end
    end
    fprintf('Finished generating %d synthetic images.\n', total_combinations);
end

function output_filename = generate_filename(patient_id, sequence, TE, alpha, TR, TI)
    switch sequence
        case 'GRE'
            output_filename = sprintf('%s_%s_synthim_TE%d_FA%d_TR%d.nii', ...
                patient_id, sequence, TE, alpha, TR);
        case 'MP2RAGE'
            output_filename = sprintf('%s_%s_synthim_TE%d_FA%d_TR%d_TI1_%d_TI2_%d.nii', ...
                patient_id, sequence, TE, alpha, TR, TI(1), TI(2));
        case 'FSE'
            output_filename = sprintf('%s_%s_synthim_TE%d_TR%d.nii', ...
                patient_id, sequence, TE, TR);
        case 'FLAIR'
            output_filename = sprintf('%s_%s_synthim_TE%d_TR%d_TI%d.nii', ...
                patient_id, sequence, TE, TR, TI);
        otherwise
            output_filename = sprintf('%s_%s_synthim_TE%d_FA%d_TR%d.nii', ...
                patient_id, sequence, TE, alpha, TR);
    end
end