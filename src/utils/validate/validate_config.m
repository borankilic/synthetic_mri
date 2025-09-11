function [cfg, flags] = validate_config(cfg)
% Validates configuration parameters and paths
% Returns modified config and validation flags struct

flags = struct();

fprintf('=== Configuration Validation ===\n');

%% Path Validation
fprintf('\nValidating paths...\n');

% Check patient directory
flags.patient_dir_exists = check_path(cfg.paths.patient_dir, 'patient_dir', true);

% Check k-space file
flags.kspace_exists = check_path(cfg.paths.kspace_path, 'kspace_path', false);

% Check magnitude directory
flags.im_mag_dir_exists = check_path(cfg.paths.im_mag_dir, 'im_mag_dir', true);

% Check phase directory  
flags.im_phase_dir_exists = check_path(cfg.paths.im_phase_dir, 'im_phase_dir', true);

% Handle output directory (create copy if files exist)
% cfg.paths.output_dir = handle_output_dir(cfg.paths.output_dir);
% flags.output_dir_ready = true;

% Check ESPIRiT path
flags.espirit_exists = check_path(cfg.paths.espirit_path, 'espirit_path', true);

% Check mapVBVD path
flags.mapvbvd_exists = check_path(cfg.paths.mapvbvd_path, 'mapvbvd_path', false);

% Check ROMEO script
flags.romeo_exists = check_path(cfg.paths.romeo_script_path, 'romeo_script_path', false);

% Will Implement distortion field estimation
flags.imp_B0_corr = flags.im_mag_dir_exists &&  flags.im_phase_dir_exists && flags.romeo_exists;

% Will Implement Synthetic kspace generation
flags.imp_kspace = flags.kspace_exists && flags.mapvbvd_exists && flags.espirit_exists;

%% Parameter Validation
fprintf('\nValidating parameters...\n');

% Validate sequence
valid_sequences = {'GRE', 'MP2RAGE', 'FSE', 'FLAIR'};
[cfg.params.sequence, flags.sequence_valid] = validate_string_param(...
    cfg.params.sequence, valid_sequences, 'GRE', 'sequence');

% Validate extract_brain option
valid_extract = {'bet', 'tissue_maps', 'no'};
[cfg.options.extract_brain, flags.extract_brain_valid] = validate_string_param(...
    cfg.options.extract_brain, valid_extract, 'tissue_maps', 'extract_brain');

% Validate interpolation method
valid_interp = {'trilinear', 'linear', 'nearest', 'cubic'};
[cfg.options.interp_method, flags.interp_method_valid] = validate_string_param(...
    cfg.options.interp_method, valid_interp, 'trilinear', 'interp_method');

% Validate show_figures option
valid_figures = {'none', 'essential', 'detailed'};
[cfg.params.smap_params.show_figures, flags.show_figures_valid] = validate_string_param(...
    cfg.params.smap_params.show_figures, valid_figures, 'none', 'smap_params.show_figures');

%% Threshold Validation
fprintf('\nValidating thresholds...\n');

% Threshold for tissue maps
[cfg.options.threshold_tmap, flags.threshold_tmap_valid] = validate_threshold(...
    cfg.options.threshold_tmap, 'threshold_tmap', 0.95);

% Threshold for BET
[cfg.options.threshold_bet, flags.threshold_bet_valid] = validate_threshold(...
    cfg.options.threshold_bet, 'threshold_bet', 0.5);

% ESPIRiT thresholds
[cfg.params.smap_params.eigThresh1, flags.eigThresh1_valid] = validate_threshold(...
    cfg.params.smap_params.eigThresh1, 'smap_params.eigThresh1', 0.1);

[cfg.params.smap_params.eigThresh2, flags.eigThresh2_valid] = validate_threshold(...
    cfg.params.smap_params.eigThresh2, 'smap_params.eigThresh2', 0.8);

%% Sequence-specific Parameter Validation
fprintf('\nValidating sequence-specific parameters...\n');
flags.sequence_params_valid = validate_sequence_params(cfg);

%% Summary
fprintf('\n=== Validation Summary ===\n');
all_flags = struct2cell(flags);
valid_count = sum([all_flags{:}]);
total_count = length(all_flags);
fprintf('Valid parameters: %d/%d\n', valid_count, total_count);

if valid_count == total_count
    fprintf('✓ All parameters are valid\n');
else
    fprintf('⚠ Some parameters had issues (see above)\n');
end

end

%% Helper Functions

function exists_flag = check_path(path_str, name, is_dir)
    if is_dir
        if isfolder(path_str)
            % Check if directory has files
            files = dir(fullfile(path_str, '**', '*'));
            files = files(~[files.isdir]); % Remove . and ..
            if isempty(files)
                fprintf('⚠ Warning: %s exists but is empty: %s\n', name, path_str);
                exists_flag = false;
            else
                fprintf('✓ %s exists and contains %d files\n', name, length(files));
                exists_flag = true;
            end
        else
            fprintf('✗ Error: %s directory not found: %s\n', name, path_str);
            exists_flag = false;
        end
    else
        if isfile(path_str)
            fprintf('✓ %s file exists\n', name);
            exists_flag = true;
        else
            fprintf('✗ Error: %s file not found: %s\n', name, path_str);
            exists_flag = false;
        end
    end
end

function output_dir = handle_output_dir(original_dir)
    if ~isfolder(original_dir)
        mkdir(original_dir);
        fprintf('✓ Created output directory: %s\n', original_dir);
        output_dir = original_dir;
        return;
    end
    
    % Check if directory has files
    files = dir(fullfile(original_dir, '*'));
    files = files(~[files.isdir]);
    
    if isempty(files)
        fprintf('✓ Output directory exists and is empty: %s\n', original_dir);
        output_dir = original_dir;
    else
        % Create numbered copy
        [path_part, name, ext] = fileparts(original_dir);
        counter = 2;
        while true
            new_dir = sprintf('%s(%d)%s', path_part, counter, '/');
            if ~isfolder(new_dir)
                mkdir(new_dir);
                fprintf('⚠ Output directory had files, created: %s\n', new_dir);
                output_dir = new_dir;
                break;
            end
            counter = counter + 1;
        end
    end
end

function [value, is_valid] = validate_string_param(current_value, valid_options, default_value, param_name)
    if ismember(current_value, valid_options)
        fprintf('✓ %s: %s\n', param_name, current_value);
        value = current_value;
        is_valid = true;
    else
        fprintf('⚠ Warning: Invalid %s "%s", using default "%s"\n', param_name, current_value, default_value);
        value = default_value;
        is_valid = false;
    end
end

function [value, is_valid] = validate_threshold(current_value, param_name, default_value)
    if isnumeric(current_value) && current_value >= 0 && current_value <= 1
        fprintf('✓ %s: %.3f\n', param_name, current_value);
        value = current_value;
        is_valid = true;
    else
        fprintf('⚠ Warning: %s (%.3f) not in [0,1], using default (%.3f)\n', param_name, current_value, default_value);
        value = default_value;
        is_valid = false;
    end
end

function is_valid = validate_sequence_params(cfg)
    is_valid = true;
    
    % Sequence-specific TI validation
    switch cfg.params.sequence
        case 'MP2RAGE'
            if ~iscell(cfg.params.TI_values)
                fprintf('✗ Error: MP2RAGE requires TI_values as cell array of {TI1; TI2} pairs\n');
                is_valid = false;
            else
                % Check each cell contains exactly 2 values
                for i = 1:length(cfg.params.TI_values)
                    if ~iscell(cfg.params.TI_values{i}) || length(cfg.params.TI_values{i}) ~= 2
                        fprintf('✗ Error: MP2RAGE TI_values{%d} must contain exactly 2 values {TI1; TI2}\n', i);
                        is_valid = false;
                    end
                end
                if is_valid
                    fprintf('✓ MP2RAGE TI_values format correct (cell array of pairs)\n');
                end
            end
            
        case {'FLAIR'}
            if iscell(cfg.params.TI_values)
                fprintf('✗ Error: %s requires TI_values as numeric array, not cell array\n', cfg.params.sequence);
                is_valid = false;
            elseif ~isnumeric(cfg.params.TI_values)
                fprintf('✗ Error: %s requires TI_values as numeric array\n', cfg.params.sequence);
                is_valid = false;
            else
                fprintf('✓ %s TI_values format correct\n', cfg.params.sequence);
            end
            
    end
    
    % Check if parameter arrays are non-empty
    if isempty(cfg.params.TE_values)
        fprintf('✗ Error: TE_values cannot be empty\n');
        is_valid = false;
    end
    
    if isempty(cfg.params.TR_values)
        fprintf('✗ Error: TR_values cannot be empty\n');
        is_valid = false;
    end
    
    if isempty(cfg.params.FA_values)
        fprintf('✗ Error: FA_values cannot be empty\n');
        is_valid = false;
    end
end