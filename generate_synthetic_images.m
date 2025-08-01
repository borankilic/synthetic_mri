function generate_synthetic_images(patient_dir, output_dir, varargin)
%GENERATE_SYNTHETIC_IMAGES Generate synthetic MRI images using Ernst equation
%   generate_synthetic_images(patient_dir, output_dir) generates synthetic
%   images using default parameters.
%   
%   generate_synthetic_images(patient_dir, output_dir, 'Name', Value) 
%   specifies additional parameters:
%   - 'TE_vals': Echo times in ms (default: [10])
%   - 'TR_vals': Repetition times in ms (default: [500])
%   - 'alpha_vals': Flip angles in degrees (default: [5])
%   - 'signal_constant': Signal scaling constant (default: 2000)
%   - 'interp_method': Interpolation method for B1 map (default: 'linear')

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'patient_dir', @(x) isfolder(x));
    addRequired(p, 'output_dir', @(x) true); % Will be created if doesn't exist
    addParameter(p, 'TE_vals', [10], @(x) isnumeric(x) && all(x > 0));
    addParameter(p, 'TR_vals', [500], @(x) isnumeric(x) && all(x > 0));
    addParameter(p, 'alpha_vals', [5], @(x) isnumeric(x) && all(x > 0));
    addParameter(p, 'signal_constant', 2000, @(x) isnumeric(x) && x > 0);
    addParameter(p, 'interp_method', 'trilinear', @(x) ischar(x) || isstring(x));
    addParameter(p, 'extract_brain', 'fsl', @(x) ischar(x) || isstring(x));
    
    
    parse(p, patient_dir, output_dir, varargin{:});
    
    % Extract parameters
    TE_vals = p.Results.TE_vals;
    TR_vals = p.Results.TR_vals;
    alpha_vals = p.Results.alpha_vals;
    signal_constant = p.Results.signal_constant;
    interp_method = p.Results.interp_method;
    extract_brain = p.Results.extract_brain;
    
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Create file structure
    data_struct = create_data_struct(patient_dir);


    %Apply brain mask to quantitative maps
    extracted_mpms_temp = fullfile(output_dir, 'extracted_mpms/');
    if ~exist(extracted_mpms_temp, 'dir'), mkdir(extracted_mpms_temp), end

    if strcmp(extract_brain, 'bet')
        % Brain extraction using BET
        fprintf("Extracting brain using BET.\n");
        extract_brain_bet({data_struct.R1.filepath, data_struct.R2.filepath, data_struct.PD.filepath}, extracted_mpms_temp); 
    
    elseif strcmp(extract_brain, 'fsl')
        % Brain extraction using FSL
        fprintf("Extracting brain using tissue prob maps and applying erosion, dilation to the mask using fslmaths.\n");
        extract_brain_fsl(data_struct, extracted_mpms_temp);
    
    else
        fprintf("Brain extraction skipped.\n");
    end


    %Register the B1 maps to ref maps
    registered_B1_temp = fullfile(output_dir, 'registered_B1/');
    if ~exist(registered_B1_temp, 'dir'), mkdir(registered_B1_temp), end
    register_b1_to_ref(data_struct.B1map.filepath, data_struct.B1ref.filepath, data_struct.R1.filepath, registered_B1_temp, interp_method);
    register_b1_to_ref(data_struct.smap.filepath, data_struct.smap.filepath, data_struct.PD.filepath, registered_B1_temp, "spline")
    
    % smap_resized  = resize_b1_map(data_struct.smap.filepath, data_struct.PD.filepath, interp_method); 
    % [~, name] = fileparts(data_struct.smap.name);
    % smap_filepath = [registered_B1_temp name '_registered4' '.nii.gz'];
    % niftiwrite(smap_resized, smap_filepath);

    data_struct = create_data_struct(output_dir, data_struct);
    
    % Load MRI data
    fprintf('Loading MRI data from: %s\n', output_dir);
    data_struct = load_mri_data(data_struct);
    patient_id = erase(data_struct.R1.name, '_R1.nii');
    if contains(patient_id, '.nii'), erase(patient_dir, '.nii'), end 

    % Apply Gaussian blur to the smap
    sigma = 2; % Or [sigma_x, sigma_y, sigma_z] for anisotropic blur
    data_struct.smap_registered.data = imgaussfilt3(data_struct.smap_registered.data, sigma);
    
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
                    signal_constant, interp_method, extract_brain);
                
                % Save result
                output_filename = sprintf('%s_synthim_TE%d_FA%d_TR%d.nii', ...
                    patient_id, TE, alpha, TR);
                output_path = fullfile(output_dir, output_filename);
                
                save_nifti_image(signal, output_path, data_struct.R1.ref_header);
                
                fprintf('Saved: %s\n', output_filename);
            end
        end
    end
    
    fprintf('Finished generating %d synthetic images.\n', total_combinations);
end