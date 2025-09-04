function write_json_metadata(output_path, ref_header, TE, TR, FA)
%WRITE_JSON_METADATA Write JSON metadata file with sequence parameters and NIfTI header info
%   write_json_metadata(output_path, ref_header, TE, TR, FA) creates a JSON
%   file containing sequence parameters (TE, TR, FA) and NIfTI header
%   information extracted from ref_header.
%   
%   Inputs:
%   - output_path: Full path for the output JSON file
%   - ref_header: NIfTI header struct from niftiinfo()
%   - TE: Echo time in milliseconds
%   - TR: Repetition time in milliseconds  
%   - FA: Flip angle in degrees
%   
%   The function writes sequence parameters with both short and long key names,
%   plus hMRI toolbox compatible keys:
%   - "TE", "echo_time", and "EchoTime"
%   - "TR", "repetition_time", and "RepetitionTimeExcitation" 
%   - "FA", "flip_angle", and "FlipAngle"

    % Input validation
    if ~ischar(output_path) && ~isstring(output_path)
        error('output_path must be a character array or string');
    end
    if ~isstruct(ref_header)
        error('ref_header must be a struct');
    end
    if ~isnumeric(TE) || TE <= 0
        error('TE must be a positive numeric value');
    end
    if ~isnumeric(TR) || TR <= 0
        error('TR must be a positive numeric value');
    end
    if ~isnumeric(FA) || FA <= 0
        error('FA must be a positive numeric value');
    end
    
    % Ensure output path has .json extension
    [filepath, name, ext] = fileparts(output_path);
    if ~strcmpi(ext, '.json')
        output_path = fullfile(filepath, [name, '.json']);
    end
    
    % Create output directory if it doesn't exist
    if ~isempty(filepath) && ~exist(filepath, 'dir')
        mkdir(filepath);
    end
    
    % Initialize metadata structure with sequence parameters at the top
    metadata = struct();
    
    % Add sequence parameters with multiple key formats for compatibility
    % Original keys
    metadata.TE = TE;
    metadata.echo_time = TE;
    metadata.TR = TR;
    metadata.repetition_time = TR;
    metadata.FA = FA;
    metadata.flip_angle = FA;
    
    % hMRI toolbox compatible keys (convert from ms to seconds for TE and TR)
    metadata.EchoTime = TE / 1000;  % Convert ms to seconds
    metadata.RepetitionTimeExcitation = TR / 1000;  % Convert ms to seconds
    metadata.FlipAngle = FA;  % Already in degrees
    
    % Add units for clarity
    metadata.TE_units = 'ms';
    metadata.echo_time_units = 'ms';
    metadata.TR_units = 'ms';
    metadata.repetition_time_units = 'ms';
    metadata.FA_units = 'degrees';
    metadata.flip_angle_units = 'degrees';
    
    % Add generation timestamp
    metadata.generated_timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    metadata.matlab_version = version;
    
    % Extract and add relevant NIfTI header information
    nifti_info = extract_nifti_metadata(ref_header);
    metadata.nifti_header = nifti_info;
    
    % Convert struct to JSON and write to file
    try
        % Use MATLAB's built-in jsonencode function (available from R2016b)
        json_str = jsonencode(metadata);
        
        % Pretty print the JSON for better readability
        json_str = prettify_json(json_str);
        
        % Write to file
        fid = fopen(output_path, 'w');
        if fid == -1
            error('Could not open file for writing: %s', output_path);
        end
        
        fprintf(fid, '%s', json_str);
        fclose(fid);
        
        fprintf('JSON metadata written to: %s\n', output_path);
        
    catch ME
        error('Failed to write JSON metadata: %s', ME.message);
    end
end

function nifti_info = extract_nifti_metadata(ref_header)
%EXTRACT_NIFTI_METADATA Extract relevant metadata from NIfTI header
    
    nifti_info = struct();
    
    % Image dimensions and voxel information
    if isfield(ref_header, 'ImageSize')
        nifti_info.image_size = ref_header.ImageSize;
    end
    
    if isfield(ref_header, 'PixelDimensions')
        nifti_info.voxel_size = ref_header.PixelDimensions;
    end
    
    if isfield(ref_header, 'SpaceUnits')
        nifti_info.space_units = ref_header.SpaceUnits;
    end
    
    if isfield(ref_header, 'TimeUnits')
        nifti_info.time_units = ref_header.TimeUnits;
    end
    
    % Data type information
    if isfield(ref_header, 'Datatype')
        nifti_info.datatype = ref_header.Datatype;
    end
    
    if isfield(ref_header, 'BitsPerPixel')
        nifti_info.bits_per_pixel = ref_header.BitsPerPixel;
    end
    
    % Spatial transformation matrix
    if isfield(ref_header, 'Transform')
        nifti_info.transform_matrix = ref_header.Transform.T;
        nifti_info.transform_dimensionality = ref_header.Transform.Dimensionality;
    end
    
    % Description and other metadata
    if isfield(ref_header, 'Description')
        nifti_info.description = ref_header.Description;
    end
    
    if isfield(ref_header, 'Version')
        nifti_info.nifti_version = ref_header.Version;
    end
    
    % Raw filename if available
    if isfield(ref_header, 'Filename')
        [~, name, ext] = fileparts(ref_header.Filename);
        nifti_info.source_filename = [name, ext];
    end
    
    % Additional fields that might be present
    optional_fields = {'SliceCode', 'FrequencyDimension', 'PhaseDimension', ...
                      'SliceDimension', 'SliceStart', 'SliceEnd', 'SliceDuration'};
    
    for i = 1:length(optional_fields)
        field = optional_fields{i};
        if isfield(ref_header, field)
            nifti_info.(lower(field)) = ref_header.(field);
        end
    end
end

function pretty_json = prettify_json(json_str)
%PRETTIFY_JSON Add indentation and line breaks to JSON string for readability
    
    % Simple prettification - add newlines and indentation
    pretty_json = '';
    indent_level = 0;
    in_string = false;
    
    for i = 1:length(json_str)
        char = json_str(i);
        
        % Track if we're inside a string (to avoid formatting string content)
        if char == '"' && (i == 1 || json_str(i-1) ~= '\')
            in_string = ~in_string;
        end
        
        if ~in_string
            switch char
                case '{'
                    pretty_json = [pretty_json, char, newline];
                    indent_level = indent_level + 1;
                    pretty_json = [pretty_json, repmat('  ', 1, indent_level)];
                    
                case '}'
                    if ~isempty(pretty_json) && pretty_json(end) ~= newline
                        pretty_json = [pretty_json, newline];
                    end
                    indent_level = indent_level - 1;
                    pretty_json = [pretty_json, repmat('  ', 1, indent_level), char];
                    
                case ','
                    pretty_json = [pretty_json, char, newline];
                    pretty_json = [pretty_json, repmat('  ', 1, indent_level)];
                    
                otherwise
                    pretty_json = [pretty_json, char];
            end
        else
            pretty_json = [pretty_json, char];
        end
    end
    
    % Clean up any trailing whitespace
    pretty_json = strtrim(pretty_json);
end