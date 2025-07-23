function extract_brain_bet(files_list, output_dir)
%RUN_BRAIN_EXTRACTION_BATCH Applies FSL's BET to all NIfTI files in a directory
%   run_brain_extraction_batch(input_dir, output_dir)
%   - input_dir: directory containing .nii or .nii.gz files
%   - output_dir: where brain-extracted images will be saved


    % 
    % % Check input directory
    % if ~exist(input_dir, 'dir')
    %     error('Input directory does not exist: %s', input_dir);
    % end
    % 
    % % Create output directory if it doesn't exist
    % if ~exist(output_dir, 'dir')
    %     mkdir(output_dir);
    % end
    % 
    % % Find all NIfTI files (.nii and .nii.gz)
    % nii_files = [dir(fullfile(input_dir, '*.nii')), dir(fullfile(input_dir, '*.nii.gz'))];
    % 
    % if isempty(nii_files)
    %     error('No NIfTI files found in: %s', input_dir);
    % end

    % Loop over each NIfTI file
    for k = 1:length(files_list)
        
        input_file = string(files_list(k));
        % Strip extensions for output filename
        [~, name, ext] = fileparts(input_file);
        if strcmp(ext, '.gz') && endsWith(name, '.nii')
            name = extractBefore(name, '.nii'); % Remove .nii from .nii.gz
        end

        % Output file path
        output_file = fullfile(output_dir, [char(name) '_brain.nii.gz']);
        
        if ~isfile(output_file)
            %% Construct BET command (you can tweak -f and -g)
            cmd = sprintf('bet "%s" "%s" -f 0.5 -g 0 -m', input_file, output_file);
            fprintf('Running BET with: %s\n', cmd);
            
            %% Run the command
            fprintf('Running BET on: %s\n', name);
            [status, result] = system(cmd);

            if status == 0
                fprintf('  ✓ BET completed: %s\n', output_file);
            else
                fprintf('  ✗ BET failed for %s:\n%s\n', nii_files(k).name, result);
            end
        else
            fprintf('Extracted brain %s already exists. Skipping BET ', output_file);
        end
    end
end
