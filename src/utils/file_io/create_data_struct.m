function file_struct = create_data_struct(root_dir, file_struct)
    % FIND_MRI_FILES Searches for MRI files in specified directories and creates a struct
    %
    % Inputs:
    %   root_dir - Root directory path (string or char array)
    %   file_struct - Existing struct to update (optional, creates new if empty/missing)
    %
    % Output:
    %   file_struct - Struct with file information containing fields:
    %                 filepath, name, folder for each found file
    %
    % The function searches for .nii and .nii.gz files in:
    %   - root_dir/Results: R1, R2, PD
    %   - root_dir/B1mapCalc: B1map, B1ref, smap
    %   - root_dir/MPMCalc: c1, c2, c3
    %   - root_dir/extracted_mpms: R1_brain, R1_brain_mask, R2*_brain, R2*_brain_mask, PD_brain, PD_brain_mask, fsl_mask
    %   - root_dir/registered_B1: B1map_registered
    
    % Handle optional input parameter
    if nargin < 2 || isempty(file_struct)
        file_struct = struct();
    end
    
    % Define search configurations: {folder, {patterns}, {field_names}}
    search_configs = {
        {'Results', {'R1', 'R2', 'PD'}, {'R1', 'R2', 'PD'}}, ...
        {'B1mapCalc', {'B1map', 'B1ref', 'smap'}, {'B1map', 'B1ref', 'smap'}}, ...
        {'MPMCalc', {'c1', 'c2', 'c3'}, {'c1', 'c2', 'c3'}}, ...
        {'extracted_mpms', {'R1_brain', 'R1_brain_mask', 'R2.*_brain', 'R2.*_brain_mask', 'PD_brain', 'PD_brain_mask', 'fsl_mask'}, ...
                          {'R1_brain', 'R1_brain_mask', 'R2_brain', 'R2_brain_mask', 'PD_brain', 'PD_brain_mask', 'fsl_mask'}}, ...
        {'registered_B1', {'B1map_registered', 'smap_registered'}, {'B1map_registered','smap_registered'} }, ...
        {'sens_maps', {'smap'}, {'smap'}}
    };
    
    % Process each search configuration
    for config_idx = 1:length(search_configs)
        subfolder   = search_configs{config_idx}{1};
        patterns    = search_configs{config_idx}{2};
        field_names = search_configs{config_idx}{3};
        
        % Construct full path to subfolder
        search_dir = fullfile(root_dir, subfolder);
        
        % Check if directory exists
        if ~exist(search_dir, 'dir')
            fprintf('Skipping - Directory does not exist: %s\n', search_dir);
            continue;
        end
        
        % Search for each pattern in the current subfolder
        for pattern_idx = 1:length(patterns)
            pattern    = patterns{pattern_idx};
            field_name = field_names{pattern_idx};
            
            % Skip if this field is already filled in the struct
            if isfield(file_struct, field_name) && ~isempty(file_struct.(field_name))
                fprintf('Skipping %s - already exists in struct\n', field_name);
                continue;
            end
            
            % Find files matching the pattern in the folder
            found_files = search_files_in_dir(search_dir, pattern);
            
            if isempty(found_files)
                subd_info = dir(fullfile(search_dir, [pattern '*']));
                if ~isempty(subd_info) && subd_info.isdir
                    subd_path = fullfile(subd_info(1).folder, subd_info(1).name);
                    % Read, fix, and rewrite each NIfTI in the subdirectory
                    nii_list = search_files_in_dir(subd_path, '');
                    fixed_files = {};
                    for f = nii_list
                        in_file = fullfile(subd_path, f{1}.name);
                        nii = niftiread(in_file);
                        info = niftiinfo(in_file);
                        info.Datatype = class(nii);
                        info.Description = '';
                        [~, base] = fileparts(f{1}.name);
                        out_file = fullfile(search_dir, [base '_' pattern '.nii']);
                        niftiwrite(nii, out_file, info, 'Compressed', false);
                        % Update found_files to point to the new file
                        file_info.name = [base '_' pattern '.nii'];
                        file_info.folder = search_dir;
                        file_info.filepath = out_file;
                        fixed_files{end+1} = file_info; %#ok<AGROW>
                    end
                    found_files = fixed_files;
                end
            end
            
            % Add found files to struct
            if ~isempty(found_files)
                if isscalar(found_files)
                    % Single file found
                    file_struct.(field_name) = create_file_info(found_files{1});
                    fprintf('Found %s: %s\n', field_name, found_files{1}.name);
                else
                    % Multiple files found - add with numbered suffixes
                    for file_idx = 1:length(found_files)
                        if file_idx == 1
                            current_field = field_name;
                        else
                            current_field = sprintf('%s_%d', field_name, file_idx);
                        end
                        file_struct.(current_field) = create_file_info(found_files{file_idx});
                        fprintf('Found %s: %s\n', current_field, found_files{file_idx}.name);
                    end
                end
            else
                fprintf('No files found for pattern "%s" in %s\n', pattern, search_dir);
            end
        end
    end
end

function found_files = search_files_in_dir(search_dir, pattern)
    % Search for files matching pattern in directory
    
    found_files = {};
    
    % Get all .nii and .nii.gz files
    nii_files    = dir(fullfile(search_dir, '*.nii'));
    nii_gz_files = dir(fullfile(search_dir, '*.nii.gz'));
    all_files    = [nii_files; nii_gz_files];
    
    % Filter files that match the pattern (or all if pattern is empty)
    for i = 1:length(all_files)
        if ~all_files(i).isdir
            filename = all_files(i).name;
            
            if isempty(pattern) || match_pattern(filename, pattern)
                file_info.name     = filename;
                file_info.folder   = search_dir;
                file_info.filepath = fullfile(search_dir, filename);
                found_files{end+1} = file_info; %#ok<AGROW>
            end
        end
    end
end

function matches = match_pattern(filename, pattern)
    % Check if filename matches the given pattern
    
    if contains(pattern, '.*')
        % Regex pattern matching
        matches = ~isempty(regexp(filename, pattern, 'once', 'ignorecase'));
    else
        % Simple substring matching
        matches = contains(filename, pattern, 'IgnoreCase', true);
    end
end

function file_info = create_file_info(file_data)
    % Create standardized file info structure
    
    file_info = struct();
    file_info.filepath = file_data.filepath;
    file_info.name     = file_data.name;
    file_info.folder   = file_data.folder;
end
