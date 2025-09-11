function validate_mri_data(data)
%VALIDATE_MRI_DATA Validate that all required MRI data is present
%   validate_mri_data(data) checks that the data struct contains all
%   required fields (PD, R1, R2) and optionally B1 correction data.
%   Throws an error if critical data is missing.

    % Required fields for Ernst equation
    required_fields = {'PD', 'R1', 'R2'};
    
    % Check required fields
    missing_fields = {};
    for i = 1:length(required_fields)
        field = required_fields{i};
        if ~isfield(data, field) || isempty(data.(field))
            missing_fields{end+1} = field; %#ok<AGROW>
        end
    end
    
    if ~isempty(missing_fields)
        error('Missing required MRI data: %s', strjoin(missing_fields, ', '));
    end
    
    % Check data dimensions consistency
    reference_size = size(data.PD.data);
    size_mismatches = {};
    
    for i = 1:length(required_fields)
        field = required_fields{i};
        if ~isequal(size(data.(field)), reference_size)
            size_mismatches{end+1} = sprintf('%s: %s vs reference %s', ...
                field, mat2str(size(data.(field))), mat2str(reference_size)); %#ok<AGROW>
        end
    end
    
    if ~isempty(size_mismatches)
        error('Data size mismatches detected:\n%s', strjoin(size_mismatches, '\n'));
    end
    
    % Check for B1 correction data
    if isfield(data, 'B1map') && ~isempty(data.B1map)
        fprintf('B1 correction data available (size: %s)\n', mat2str(size(data.B1map)));
    else
        warning('No B1 correction data found. Using uniform B1 = 100%%');
    end
    
    % Check for reasonable data ranges
    check_data_ranges(data);
    
    fprintf('Data validation passed for %s data\n', mat2str(reference_size));
end

function check_data_ranges(data)
%CHECK_DATA_RANGES Check if data values are in reasonable ranges
    
    % Define reasonable ranges for each parameter
    ranges = struct(...
        'PD', [0, 10000], ...      % Proton density
        'R1', [0, 10], ...         % R1 relaxation rate (1/s)
        'R2', [0, 100], ...        % R2 relaxation rate (1/s)
        'B1map', [50, 150]);       % B1 map (50-150% of nominal)
    
    field_names = fieldnames(ranges);
    
    for i = 1:length(field_names)
        field = field_names{i};
        if isfield(data, field) && ~isempty(data.(field))
            values = data.(field).data;
            valid_values = values(~isnan(values) & ~isinf(values) & values ~= 0);
            
            if ~isempty(valid_values)
                min_val = min(valid_values(:));
                max_val = max(valid_values(:));
                expected_range = ranges.(field);
                
                if min_val < expected_range(1) || max_val > expected_range(2)
                    warning(['%s values may be outside expected range [%.2f, %.2f]. ' ...
                        'Found range: [%.2f, %.2f]'], ...
                        field, expected_range(1), expected_range(2), min_val, max_val);
                end
            end
        end
    end
end