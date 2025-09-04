function metadata = write_comprehensive_json(input_nii_path, varargin)
%WRITE_COMPREHENSIVE_JSON Write comprehensive JSON metadata with automatic parameter extraction
%   write_comprehensive_json(input_nii_path) creates a JSON file with the same
%   name as the input NIfTI file, extracting header information automatically.
%   
%   write_comprehensive_json(input_nii_path, 'param', value) accepts parameter pairs:
%
%   Parameters:
%       'ref_nii_filepath': Reference NIfTI file to extract missing parameters from
%                          its corresponding JSON file (optional)
%       'TE':              Echo time in milliseconds (optional)
%       'TR':              Repetition time in milliseconds (optional)
%       'FA':              Flip angle in degrees (optional)
%       'fft_details':     FFT processing details struct (optional)
%       'sequence_params': Additional sequence parameters struct (optional)
%       'component_type':  Data component ('magnitude', 'phase', 'real', 'complex') (optional)
%       'data_range':      [min, max] values of the data (optional)
%       'processing_info': Custom processing information struct (optional)
%       'bids_compatible': Generate BIDS-compatible keys (true/false)
%                          Default: true
%       'hmri_compatible': Generate hMRI toolbox compatible keys (true/false)
%                          Default: true
%       'pretty_print':    Format JSON with indentation (true/false)
%                          Default: true
%       'verbose':         Display progress messages (true/false)
%                          Default: true
%
%   Examples:
%       % Basic usage - extracts header info automatically
%       write_comprehensive_json('processed_image.nii');
%       
%       % With sequence parameters
%       write_comprehensive_json('result.nii', 'TE', 30, 'TR', 2000, 'FA', 15);
%       
%       % Using reference file to fill missing parameters
%       write_comprehensive_json('kspace_mag.nii', 'ref_nii_filepath', 'original.nii');
%       
%       % With FFT details and reference
%       write_comprehensive_json('processed.nii', 'fft_details', fft_info, ...
%                               'ref_nii_filepath', 'source.nii');

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'input_nii_path', @(x) ischar(x) || isstring(x));
    addParameter(p, 'ref_nii_filepath', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'TE', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
    addParameter(p, 'TR', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
    addParameter(p, 'FA', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
    addParameter(p, 'fft_details', [], @(x) isempty(x) || isstruct(x));
    addParameter(p, 'sequence_params', struct(), @isstruct);
    addParameter(p, 'component_type', '', @ischar);
    addParameter(p, 'data_range', [], @(x) isempty(x) || (isnumeric(x) && length(x) == 2));
    addParameter(p, 'processing_info', struct(), @isstruct);
    addParameter(p, 'bids_compatible', true, @islogical);
    addParameter(p, 'hmri_compatible', true, @islogical);
    addParameter(p, 'pretty_print', true, @islogical);
    addParameter(p, 'verbose', true, @islogical);
    
    parse(p, input_nii_path, varargin{:});
    
    % Extract parameters
    ref_nii_filepath = p.Results.ref_nii_filepath;
    TE = p.Results.TE;
    TR = p.Results.TR;
    FA = p.Results.FA;
    fft_details = p.Results.fft_details;
    sequence_params = p.Results.sequence_params;
    component_type = p.Results.component_type;
    data_range = p.Results.data_range;
    processing_info = p.Results.processing_info;
    bids_compatible = p.Results.bids_compatible;
    hmri_compatible = p.Results.hmri_compatible;
    pretty_print = p.Results.pretty_print;
    verbose = p.Results.verbose;
    
    % Validate input NIfTI file exists
    if ~exist(input_nii_path, 'file')
        error('Input NIfTI file not found: %s', input_nii_path);
    end
    
    % Generate output JSON path (same name as input NIfTI but with .json extension)
    [filepath, name, ~] = fileparts(input_nii_path);
    output_path = fullfile(filepath, [name '.json']);
    
    if verbose
        fprintf('Processing NIfTI file: %s\n', input_nii_path);
        fprintf('Output JSON will be: %s\n', output_path);
    end
    
    % Load reference parameters from reference JSON file if provided
    ref_params = struct();
    if ~isempty(ref_nii_filepath)
        ref_params = load_reference_parameters(ref_nii_filepath, verbose);
    end
    
    % Extract NIfTI header information
    try
        nii_header = niftiinfo(input_nii_path);
        if verbose
            fprintf('Successfully extracted NIfTI header information\n');
        end
    catch ME
        error('Failed to read NIfTI header from %s: %s', input_nii_path, ME.message);
    end
    
    % Merge parameters: provided parameters override reference parameters
    merged_params = merge_parameters(ref_params, TE, TR, FA, sequence_params, verbose);
    
    % Initialize metadata structure
    metadata = struct();
    
    %% SEQUENCE PARAMETERS SECTION
    seq_section = add_sequence_parameters(merged_params.TE, merged_params.TR, ...
                                         merged_params.FA, merged_params.sequence_params, ...
                                         bids_compatible, hmri_compatible);
    metadata = merge_structs(metadata, seq_section);
    
    %% DATA INFORMATION SECTION  
    data_section = add_data_information(component_type, data_range, input_nii_path);
    metadata = merge_structs(metadata, data_section);
    
    %% FFT/PROCESSING SECTION
    if ~isempty(fft_details)
        fft_section = add_fft_processing_info(fft_details);
        metadata = merge_structs(metadata, fft_section);
    elseif isfield(ref_params, 'fft_details')
        fft_section = add_fft_processing_info(ref_params.fft_details);
        metadata = merge_structs(metadata, fft_section);
    end
    
    if ~isempty(fieldnames(processing_info))
        metadata.processing_details = processing_info;
    elseif isfield(ref_params, 'processing_info')
        metadata.processing_details = ref_params.processing_info;
    end
    
    %% NIFTI HEADER SECTION
    nifti_section = add_nifti_header_info(nii_header, input_nii_path);
    metadata = merge_structs(metadata, nifti_section);
    
    %% PROVENANCE SECTION
    provenance_section = add_provenance_info(input_nii_path, ref_nii_filepath);
    metadata = merge_structs(metadata, provenance_section);
    
    %% WRITE JSON FILE
    write_json_to_file(metadata, output_path, pretty_print, verbose);
    
    if verbose
        fprintf('JSON metadata successfully written to: %s\n', output_path);
    end
end

function ref_params = load_reference_parameters(ref_nii_filepath, verbose)
%LOAD_REFERENCE_PARAMETERS Load parameters from reference JSON file
    
    ref_params = struct();
    
    if ~exist(ref_nii_filepath, 'file')
        if verbose
            fprintf('Warning: Reference NIfTI file not found: %s\n', ref_nii_filepath);
        end
        return;
    end
    
    % Look for corresponding JSON file
    [ref_filepath, ref_name, ~] = fileparts(ref_nii_filepath);
    ref_json_path = fullfile(ref_filepath, [ref_name '.json']);
    
    if ~exist(ref_json_path, 'file')
        if verbose
            fprintf('Warning: Reference JSON file not found: %s\n', ref_json_path);
        end
        return;
    end
    
    try
        % Read and parse reference JSON
        json_text = fileread(ref_json_path);
        ref_data = jsondecode(json_text);
        
        if verbose
            fprintf('Successfully loaded reference parameters from: %s\n', ref_json_path);
        end
        
        % Extract relevant parameters
        ref_params = extract_reference_params(ref_data);
        
    catch ME
        if verbose
            fprintf('Warning: Failed to load reference JSON file: %s\n', ME.message);
        end
    end
end

function params = extract_reference_params(ref_data)
%EXTRACT_REFERENCE_PARAMS Extract relevant parameters from reference JSON data
    
    params = struct();
    
    % Extract TE (try multiple field names)
    te_fields = {'TE', 'echo_time', 'EchoTime'};
    for i = 1:length(te_fields)
        if isfield(ref_data, te_fields{i})
            te_value = ref_data.(te_fields{i});
            % Convert to milliseconds if in seconds
            if te_value < 1  % Assume seconds if < 1
                params.TE = te_value * 1000;
            else
                params.TE = te_value;
            end
            break;
        end
    end
    
    % Extract TR (try multiple field names)
    tr_fields = {'TR', 'repetition_time', 'RepetitionTime', 'RepetitionTimeExcitation'};
    for i = 1:length(tr_fields)
        if isfield(ref_data, tr_fields{i})
            tr_value = ref_data.(tr_fields{i});
            % Convert to milliseconds if in seconds
            if tr_value < 10  % Assume seconds if < 10
                params.TR = tr_value * 1000;
            else
                params.TR = tr_value;
            end
            break;
        end
    end
    
    % Extract FA (try multiple field names)
    fa_fields = {'FA', 'flip_angle', 'FlipAngle'};
    for i = 1:length(fa_fields)
        if isfield(ref_data, fa_fields{i})
            params.FA = ref_data.(fa_fields{i});
            break;
        end
    end
    
    % Extract other relevant fields
    if isfield(ref_data, 'fourier_transform') || isfield(ref_data, 'FourierTransform')
        if isfield(ref_data, 'fourier_transform')
            params.fft_details = ref_data.fourier_transform;
        else
            params.fft_details = ref_data.FourierTransform;
        end
    end
    
    if isfield(ref_data, 'processing_details')
        params.processing_info = ref_data.processing_details;
    end
    
    % Extract additional sequence parameters
    seq_param_fields = {'bandwidth', 'matrix_size', 'field_of_view', 'slice_thickness', ...
                       'number_of_slices', 'acquisition_time'};
    additional_params = struct();
    for i = 1:length(seq_param_fields)
        field = seq_param_fields{i};
        if isfield(ref_data, field)
            additional_params.(field) = ref_data.(field);
        end
    end
    
    if ~isempty(fieldnames(additional_params))
        params.sequence_params = additional_params;
    else
        params.sequence_params = struct();
    end
end

function merged = merge_parameters(ref_params, TE, TR, FA, sequence_params, verbose)
%MERGE_PARAMETERS Merge provided parameters with reference parameters
    
    merged = struct();
    
    % TE: provided parameter takes precedence
    if ~isempty(TE)
        merged.TE = TE;
    elseif isfield(ref_params, 'TE')
        merged.TE = ref_params.TE;
        if verbose
            fprintf('Using TE from reference: %.1f ms\n', merged.TE);
        end
    else
        merged.TE = [];
    end
    
    % TR: provided parameter takes precedence
    if ~isempty(TR)
        merged.TR = TR;
    elseif isfield(ref_params, 'TR')
        merged.TR = ref_params.TR;
        if verbose
            fprintf('Using TR from reference: %.1f ms\n', merged.TR);
        end
    else
        merged.TR = [];
    end
    
    % FA: provided parameter takes precedence
    if ~isempty(FA)
        merged.FA = FA;
    elseif isfield(ref_params, 'FA')
        merged.FA = ref_params.FA;
        if verbose
            fprintf('Using FA from reference: %.1f degrees\n', merged.FA);
        end
    else
        merged.FA = [];
    end
    
    % Sequence parameters: merge both, with provided taking precedence
    merged.sequence_params = struct();
    if isfield(ref_params, 'sequence_params')
        merged.sequence_params = ref_params.sequence_params;
    end
    
    % Override with provided sequence parameters
    provided_fields = fieldnames(sequence_params);
    for i = 1:length(provided_fields)
        field = provided_fields{i};
        merged.sequence_params.(field) = sequence_params.(field);
    end
end

function data_section = add_data_information(component_type, data_range, input_nii_path)
%ADD_DATA_INFORMATION Add data type and range information
    
    data_section = struct();
    
    % Determine component type from filename if not provided
    if isempty(component_type)
        [~, filename, ~] = fileparts(input_nii_path);
        if contains(filename, '_mag', 'IgnoreCase', true)
            component_type = 'magnitude';
        elseif contains(filename, '_phase', 'IgnoreCase', true)
            component_type = 'phase';
        elseif contains(filename, 'kspace', 'IgnoreCase', true) || contains(filename, 'k-space', 'IgnoreCase', true)
            component_type = 'kspace';
        else
            component_type = 'image';
        end
    end
    
    data_section.component_type = component_type;
    data_section.data_domain = get_data_domain(component_type);
    
    if ~isempty(data_range)
        data_section.data_range = data_range;
        data_section.data_min = data_range(1);
        data_section.data_max = data_range(2);
    end
    
    % Add filename information
    [~, filename, ext] = fileparts(input_nii_path);
    data_section.filename = [filename ext];
end

function domain = get_data_domain(component_type)
%GET_DATA_DOMAIN Determine data domain based on component type
    switch lower(component_type)
        case {'magnitude', 'real', 'imag', 'image'}
            domain = 'image_space';
        case 'phase'
            domain = 'image_space';
        case {'kspace', 'k-space', 'k_space', 'magnitude_kspace', 'phase_kspace'}
            domain = 'k_space';
        otherwise
            domain = 'unknown';
    end
end

function seq_section = add_sequence_parameters(TE, TR, FA, sequence_params, bids_compatible, hmri_compatible)
%ADD_SEQUENCE_PARAMETERS Add sequence parameters with multiple naming conventions
    
    seq_section = struct();
    
    % Add TE (Echo Time)
    if ~isempty(TE)
        % Standard keys
        seq_section.TE = TE;
        seq_section.echo_time = TE;
        seq_section.TE_units = 'ms';
        
        % BIDS compatible (seconds)
        if bids_compatible
            seq_section.EchoTime = TE / 1000;
        end
        
        % hMRI compatible (seconds)
        if hmri_compatible
            seq_section.EchoTime_hmri = TE / 1000;
        end
    end
    
    % Add TR (Repetition Time)
    if ~isempty(TR)
        % Standard keys
        seq_section.TR = TR;
        seq_section.repetition_time = TR;
        seq_section.TR_units = 'ms';
        
        % BIDS compatible (seconds)
        if bids_compatible
            seq_section.RepetitionTime = TR / 1000;
        end
        
        % hMRI compatible (seconds)
        if hmri_compatible
            seq_section.RepetitionTimeExcitation = TR / 1000;
        end
    end
    
    % Add FA (Flip Angle)
    if ~isempty(FA)
        % Standard keys
        seq_section.FA = FA;
        seq_section.flip_angle = FA;
        seq_section.FA_units = 'degrees';
        
        % BIDS/hMRI compatible
        if bids_compatible || hmri_compatible
            seq_section.FlipAngle = FA;
        end
    end
    
    % Add additional sequence parameters
    additional_params = fieldnames(sequence_params);
    for i = 1:length(additional_params)
        param_name = additional_params{i};
        seq_section.(param_name) = sequence_params.(param_name);
    end
end

function fft_section = add_fft_processing_info(fft_details)
%ADD_FFT_PROCESSING_INFO Add FFT and processing information
    
    fft_section = struct();
    fft_section.fourier_transform = fft_details;
    
    % Extract key information to top level for easier access
    if isfield(fft_details, 'transformation')
        fft_section.transformation_type = fft_details.transformation;
    end
    
    if isfield(fft_details, 'fft_shift_applied')
        fft_section.fft_shift_applied = fft_details.fft_shift_applied;
    end
    
    if isfield(fft_details, 'kspace_is_complex')
        fft_section.complex_data = fft_details.kspace_is_complex;
    end
end

function nifti_section = add_nifti_header_info(nii_header, input_path)
%ADD_NIFTI_HEADER_INFO Extract and add relevant NIfTI header information
    
    nifti_section = struct();
    
    % Spatial information
    spatial_info = struct();
    if isfield(nii_header, 'ImageSize')
        spatial_info.image_dimensions = nii_header.ImageSize;
        spatial_info.number_of_dimensions = length(nii_header.ImageSize);
    end
    
    if isfield(nii_header, 'PixelDimensions')
        spatial_info.voxel_size = nii_header.PixelDimensions;
        if isfield(nii_header, 'SpaceUnits')
            spatial_info.voxel_units = nii_header.SpaceUnits;
        else
            spatial_info.voxel_units = 'mm';  % Assume mm unless specified
        end
    end
    
    if isfield(nii_header, 'Transform')
        spatial_info.transform_matrix = nii_header.Transform.T;
        spatial_info.transform_dimensionality = nii_header.Transform.Dimensionality;
    end
    
    % Data type information
    data_info = struct();
    if isfield(nii_header, 'Datatype')
        data_info.datatype = nii_header.Datatype;
        data_info.datatype_description = get_datatype_description(nii_header.Datatype);
    end
    
    if isfield(nii_header, 'BitsPerPixel')
        data_info.bits_per_pixel = nii_header.BitsPerPixel;
    end
    
    % Combine sections
    nifti_section.spatial_information = spatial_info;
    nifti_section.data_information = data_info;
    
    % Additional fields
    optional_fields = {'Description', 'Version', 'SliceCode', 'FrequencyDimension', ...
                      'PhaseDimension', 'SliceDimension'};
    
    additional_info = struct();
    for i = 1:length(optional_fields)
        field = optional_fields{i};
        if isfield(nii_header, field)
            additional_info.(lower(field)) = nii_header.(field);
        end
    end
    
    if ~isempty(fieldnames(additional_info))
        nifti_section.additional_header_info = additional_info;
    end
    
    % Source filename
    [~, name, ext] = fileparts(input_path);
    nifti_section.source_filename = [name ext];
end

function description = get_datatype_description(datatype_code)
    % Create your datatype map (keys must be numeric codes)
    datatype_map = containers.Map( ...
        {int32(2), int32(4), int32(8), int32(16), int32(64)}, ...
        {'uint8', 'int16', 'int32', 'float32', 'float64'} ...
    );

    % If datatype_code is a string (from niftiinfo), convert to numeric code
    if ischar(datatype_code) || isstring(datatype_code)
        switch char(datatype_code)
            case 'uint8',  datatype_code = 2;
            case 'int16',  datatype_code = 4;
            case 'int32',  datatype_code = 8;
            case 'single', datatype_code = 16; % float32
            case 'double', datatype_code = 64; % float64
            otherwise
                warning('Unknown datatype string: %s', datatype_code);
                description = 'Unknown';
                return;
        end
    end

    % Now ensure scalar
    if ~isscalar(datatype_code)
        warning('Non-scalar datatype code: using first element');
        datatype_code = datatype_code(1);
    end

    % Lookup
    if datatype_map.isKey(int32(datatype_code))
        description = datatype_map(int32(datatype_code));
    else
        description = 'Unknown';
    end
end
function provenance_section = add_provenance_info(input_path, ref_path)
%ADD_PROVENANCE_INFO Add processing provenance information
    
    provenance_section = struct();
    provenance_section.generated_timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    provenance_section.matlab_version = version;
    
    % Get system information
    if ispc
        provenance_section.operating_system = 'Windows';
    elseif ismac
        provenance_section.operating_system = 'macOS';
    elseif isunix
        provenance_section.operating_system = 'Linux/Unix';
    else
        provenance_section.operating_system = 'Unknown';
    end
    
    % Add function name and version info
    provenance_section.generator_function = 'write_comprehensive_json';
    provenance_section.generator_version = '2.0';
    
    % Add input file information
    [~, input_name, input_ext] = fileparts(input_path);
    provenance_section.input_file = [input_name input_ext];
    
    if ~isempty(ref_path)
        [~, ref_name, ref_ext] = fileparts(ref_path);
        provenance_section.reference_file = [ref_name ref_ext];
    end
end

function merged = merge_structs(struct1, struct2)
%MERGE_STRUCTS Merge two structures, with struct2 fields taking precedence
    merged = struct1;
    fields = fieldnames(struct2);
    for i = 1:length(fields)
        merged.(fields{i}) = struct2.(fields{i});
    end
end

function write_json_to_file(metadata, output_path, pretty_print, verbose)
%WRITE_JSON_TO_FILE Write metadata structure to JSON file
    
    try
        % Use MATLAB's built-in jsonencode
        if pretty_print
            json_str = jsonencode(metadata, 'PrettyPrint', true);
        else
            json_str = jsonencode(metadata);
        end
        
        % Write to file
        fid = fopen(output_path, 'w', 'n', 'UTF-8');
        if fid == -1
            error('Could not open file for writing: %s', output_path);
        end
        
        fprintf(fid, '%s', json_str);
        fclose(fid);
        
    catch ME
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        error('Failed to write JSON metadata: %s', ME.message);
    end
end