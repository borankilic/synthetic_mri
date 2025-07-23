function register_b1_to_ref(b1_map_path, b1_ref_path, ref_path, output_dir, interp_method)
    % Ensure FSL is in the MATLAB environment
    % setenv('FSLDIR', '/software/fsl/6.0.7.11/debian-bookworm-amd64/bin/fsl');  % Adjust this path
    % setenv('PATH', [getenv('FSLDIR') '/bin:' getenv('PATH')]);
    
    % %Look for the B1 map
    % nii_files_b1 = dir(fullfile(b1_dir, '*B1map*.nii*'));
    % if isempty(nii_files_b1)
    %     error('No B1map file found in directory: %s', b1_dir);
    % end
    % b1_map_path = fullfile(b1_dir, nii_files_b1(1).name);
    % fprintf('Found B1map file: %s\n', b1_map_path);
    % 
    % %Look for the B1 anatomical reference
    % nii_files_b1ref = dir(fullfile(b1_dir, '*B1ref*.nii*'));
    % if isempty(nii_files_b1ref)
    %     error('No B1 reference file found in directory: %s', b1_dir);
    % end
    % b1_ref_path = fullfile(b1_dir, nii_files_b1ref(1).name);
    % fprintf('Found B1 ref file: %s\n', b1_ref_path);
    % 
    % %Look for the reference file
    % nii_files_r1 = dir(fullfile(ref_dir, '*R1_brain.nii*'));
    % if isempty(nii_files_r1)
    %     error('No reference file found in directory: %s', ref_dir);
    % end
    % ref_path = fullfile(ref_dir, nii_files_r1(1).name);
    % fprintf('Found reference file: %s\n', ref_path);
    
    %Set Output file path
    [~, name, ext] = fileparts(b1_map_path);       % Strip extensions for output filename
    if strcmp(ext, '.gz') && endsWith(name, '.nii')
        name = extractBefore(name, '.nii'); % Remove .nii from .nii.gz
    end
    output_path = fullfile(output_dir, [name '_registered.nii.gz']);


    % Temporary transformation matrix file
    matfile = tempname + ".mat";

    % Step 1: Rigid-body registration using FLIRT
    cmd1 = sprintf(['FSL flirt -in "%s" -ref "%s" -out "%s" -omat "%s" ' ...
                    '-interp "%s" -dof 6'], ...
                    b1_ref_path, ref_path, output_path, matfile, interp_method);

    % Step 2: Resample B1 map to match reference exactly
    cmd2 = sprintf(['FSL flirt -in "%s" -applyxfm -init "%s" -ref "%s" ' ...
                    '-out "%s" -interp "%s"'], ...
                    b1_map_path, matfile, ref_path, output_path, interp_method);

    if ~isfile(output_path)
        % Run rigid registration
        [status1, result1] = system(cmd1);
        if status1 ~= 0
            error('FLIRT registration failed:\n%s', result1);
        end
    
        % Apply transformation and resample
        [status2, result2] = system(cmd2);
        if status2 == 0
            fprintf('Registration and resampling complete. Output: %s\n', output_path);
        else
            error('FLIRT resampling failed:\n%s', result2);
        end
    else 
        fprintf('Registered B1 %s already exists. Skipping registeration ', output_path);
    end
end
