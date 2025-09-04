function extract_brain(data_struct, output_dir, varargin)
    % EXTRACT_BRAIN Extracts the brain from MRI data using different methods
    %   extract_brain(data_struct, output_dir, 'method', method_name)
    %   - data_struct: struct containing file paths
    %   - output_dir: directory to save extracted brains
    %   - method: 'bet' or 'fsl'

    p = inputParser;
    addRequired(p, 'data_struct', @isstruct);
    addRequired(p, 'output_dir', @ischar);
    addParameter(p, 'method', 'fsl', @(x) any(validatestring(x, {'bet', 'fsl'})));
    addParameter(p, 'patient_dir', @ischar);
    parse(p, data_struct, output_dir, varargin{:});

    method = p.Results.method;
    patient_dir = p.Results.patient_dir;

    if strcmp(method, 'bet')
        fprintf('Extracting brain using BET...\n');
        files_to_process = {data_struct.R1.filepath, data_struct.R2.filepath, data_struct.PD.filepath};
        extract_brain_bet(files_to_process, output_dir, patient_dir);
    elseif strcmp(method, 'fsl')
        fprintf('Extracting brain using FSL...\n');
        extract_brain_fsl(data_struct, output_dir);
    end
end

function extract_brain_bet(files_list, output_dir, patient_dir)

    % --- First try PDw_WLS1fit_TEzero --- 
    %PDW_WLS1fit_TE_zero is a clean (weighted least squares fit of provided
    %PD images) image created by hMRI toolbox as an ointermediate step. It
    %provides greater anatomical detail than MPMs while having less fitting
    %realted artifacts.
    pd_search_dir = fullfile(patient_dir, 'MPMCalc');
    pd_files= dir(fullfile(pd_search_dir, '*PDw_WLS1fit_TEzero.ni*'));

    if ~isempty(pd_files)
        pd_file = fullfile(pd_files(1).folder, pd_files(1).name);
        fprintf('Found PDw_WLS1fit_TEzero file: %s\n', pd_file);
    else
        % --- Fallback: look in files_list for PD.nii or PD.nii.gz ---
        pd_idx = find(endsWith(files_list, 'PD.nii') | endsWith(files_list, 'PD.nii.gz'), 1);
        if isempty(pd_idx)
            error('No PDw_WLS1fit_TEzero or PD image found.');
        end
        pd_file = string(files_list(pd_idx));
        fprintf('Fallback to PD file from files_list: %s\n', pd_file);
    end

    % --- Prepare mask output path ---
    [~, pd_name, pd_ext] = fileparts(pd_file);
    if strcmp(pd_ext, '.gz') && endsWith(pd_name, '.nii')
        pd_name = extractBefore(pd_name, '.nii');
    end
    pd_output = fullfile(output_dir, [pd_name '_brain.nii.gz']);
    pd_mask   = fullfile(output_dir, [pd_name '_brain_mask.nii.gz']);

    % --- Run BET once to create mask ---
    if ~isfile(pd_mask)
        cmd = sprintf('FSL bet "%s" "%s" -f 0.5 -g 0 -m', pd_file, pd_output);
        fprintf('Running BET: %s\n', cmd);
        [status, result] = system(cmd);
        if status ~= 0
            error('BET failed on reference image:\n%s', result);
        else
            fprintf('✓ Brain mask created: %s\n', pd_mask);
        end
    else
        fprintf('Brain mask already exists: %s. Skipping BET.\n', pd_mask);
    end

    % --- Apply mask to all files in files_list ---
    for k = 1:length(files_list)
        input_file = string(files_list(k));
        [~, name, ext] = fileparts(input_file);
        if strcmp(ext, '.gz') && endsWith(name, '.nii')
            name = extractBefore(name, '.nii');
        end
        output_file = fullfile(output_dir, [char(name) '_brain.nii.gz']);

        if ~isfile(output_file)
            cmd = sprintf('FSL fslmaths "%s" -mas "%s" "%s"', input_file, pd_mask, output_file);
            fprintf('Applying mask to: %s\n', input_file);
            [status, result] = system(cmd);
            if status ~= 0
                fprintf('✗ Masking failed for %s:\n%s\n', input_file, result);
            else
                fprintf('✓ Mask applied: %s\n', output_file);
            end
        else
            fprintf('Masked file already exists: %s. Skipping.\n', output_file);
        end
    end
end



function extract_brain_fsl(data_struct, output_dir)
    threshold = 0.95;
    name = erase(data_struct.R1.name, '_R1.nii');
    mask_file = fullfile(output_dir, [char(name) 'fsl_mask.nii']);
    if ~isfile(mask_file)
        cmd1 =  sprintf('FSL fslmaths "%s" -add "%s" -add "%s" -thr %f -bin "%s" ', data_struct.c1.filepath , data_struct.c2.filepath, data_struct.c3.filepath, threshold, mask_file); 
        cmd2 = sprintf('FSL fslmaths "%s" -kernel 3D -fillh -dilM -dilM -ero -ero -fillh -fillh -dilM "%s" ', mask_file, mask_file); 
        fprintf('Creating mask on: %s\n', data_struct.c1.name);
        [status, result] = system(cmd1);
        if status == 0
            fprintf('  ✓ Mask Created: %s\n', mask_file);
        else
            fprintf('  ✗ Mask creation failed %s:\n%s\n', mask_file, result);
        end
        fprintf('Applying erosion, dilation and hole filling: %s\n', mask_file);
        [status, result] = system(cmd2);
        if status == 0
            fprintf('  ✓ Erosion, dilation completed for : %s\n', mask_file);
        else
            fprintf('  ✗ Erosion, dilation failed for %s:\n%s\n', mask_file, result);
        end
    else
        fprintf("Mask already exists. Skipping mask generation.\n")
    end
    files_list = {data_struct.R1.filepath, data_struct.R2.filepath, data_struct.PD.filepath};
    for k = 1:length(files_list)
        input_file = string(files_list(k));
        [~, name, ext] = fileparts(input_file);
        if strcmp(ext, '.gz') && endsWith(name, '.nii')
            name = extractBefore(name, '.nii');
        end
        output_file = fullfile(output_dir, [char(name) '_brain.nii']); 
        if ~isfile(output_file)
            cmd3 = sprintf('FSL fslmaths "%s" -mas "%s" "%s" ', input_file, mask_file, output_file );
            fprintf('Applying generated mask to: %s\n', input_file);
            [status, result] = system(cmd3);
            if status == 0
                fprintf('  ✓ Applying mask completed for: %s\n', input_file);
            else
                fprintf('  ✗ Applying mask failed for %s:\n%s\n', input_file, result);
            end
        else
            fprintf('Extracted brain already exists. Skipping mask application. \n')
        end
    end
end
