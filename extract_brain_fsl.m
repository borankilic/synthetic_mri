function extract_brain_fsl(data_struct, output_dir)
%RUN_BRAIN_EXTRACTION_BATCH Applies FSL's BET to all NIfTI files in a directory
%   run_brain_extraction_batch(input_dir, output_dir)
%   - input_dir: directory containing .nii or .nii.gz files
%   - output_dir: where brain-extracted images will be saved
       
    
    %% Mask Creation and adjustments
     % Set threshold for constructing binary mask from WM, GM, CSF probability maps
    threshold = 0.95;
    name = erase(data_struct.R1.name, "_R1.nii");
    mask_file = fullfile(output_dir, [char(name) 'fsl_mask.nii']);
    
    if ~isfile(mask_file)
        %Add 3 prob maps together and discard voxels whose prob sum is less
        %than the threshold
        cmd1 =  sprintf('FSL fslmaths "%s" -add "%s" -add "%s" -thr %f -bin "%s" ', data_struct.c1.filepath , data_struct.c2.filepath, data_struct.c3.filepath, threshold, mask_file);
        %Apply erosion, dilation, hole filling
        cmd2 = sprintf('FSL fslmaths "%s" -kernel 3D -fillh -dilM -dilM -ero -ero -fillh -fillh -dilM "%s" ', mask_file, mask_file);
        %Creating the mask
        fprintf('Creating mask on: %s\n', data_struct.c1.name);
        [status, result] = system(cmd1);
    
        if status == 0
            fprintf('  ✓ Mask Created: %s\n', mask_file);
        else
            fprintf('  ✗ Mask creation failed %s:\n%s\n', mask_file, result);
        end
    
        % Applying erosion, dilation and hole filling
        fprintf('Applying erosion, dilation and hole filling: %s\n', mask_file);
        [status, result] = system(cmd2);
    
        if status == 0
            fprintf('  ✓ Erosion, dilation completed for : %s\n', mask_file);
        else
            fprintf('  ✗ Erosion, dilation failed for %s:\n%s\n', mask_file, result);
        end
    else
        fprintf("Mask already exists. Skipping mask generation.")
    end
    
    % Loop over each NIfTI file
    files_list = {data_struct.R1.filepath, data_struct.R2.filepath, data_struct.PD.filepath};
    for k = 1:length(files_list)
        
        input_file = string(files_list(k));
        % Strip extensions for output filename
        [~, name, ext] = fileparts(input_file);
        if strcmp(ext, '.gz') && endsWith(name, '.nii')
            name = extractBefore(name, '.nii'); % Remove .nii from .nii.gz
        end

        % Output file path
        
        output_file = fullfile(output_dir, [char(name) '_brain.nii']); 
        cmd3 = sprintf('FSL fslmaths "%s" -mas "%s" "%s" ', input_file, mask_file, output_file );
        
        % Applying the mask to MPMs
        fprintf('Applying generated mask to: %s\n', input_file);
        [status, result] = system(cmd3);
    
        if status == 0
            fprintf('  ✓ Applying mask completed for: %s\n', input_file);
        else
            fprintf('  ✗ Applying mask failed for %s:\n%s\n', input_file, result);
        end


    end
end
