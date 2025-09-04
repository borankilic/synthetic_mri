function save_nifti(image_data, output_path, ref_header, varargin)
%SAVE_NIFTI Enhanced NIfTI saving with robust fallback methods
%   save_nifti(image_data, output_path, ref_header) saves image data as NIfTI 
%   file with multiple fallback methods including minimal header and headerless options.
%
%   For complex data, automatically saves magnitude and phase components as
%   separate files with '_mag' and '_phase' suffixes.
%
%   save_nifti(..., 'param', value) accepts parameter pairs:
%
%   Parameters:
%       'data_type':        Output data type ('single', 'double', 'int16', 'uint8')
%                          Default: 'single'
%       'description':      Custom description for NIfTI header
%                          Default: auto-generated from filename
%       'verbose':          Display progress messages (true/false)
%                          Default: true
%       'create_dir':       Create output directory if needed (true/false)
%                          Default: true

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'image_data', @isnumeric);
    addRequired(p, 'output_path', @(x) ischar(x) || isstring(x));
    addRequired(p, 'ref_header', @(x) isstruct(x) || isempty(x));
    addParameter(p, 'data_type', 'single', @ischar);
    addParameter(p, 'description', '', @ischar);
    addParameter(p, 'verbose', true, @islogical);
    addParameter(p, 'create_dir', true, @islogical);
    
    parse(p, image_data, output_path, ref_header, varargin{:});
    
    data_type = p.Results.data_type;
    description = p.Results.description;
    verbose = p.Results.verbose;
    create_dir = p.Results.create_dir;
    
    % Input validation
    if isempty(image_data)
        error('image_data cannot be empty');
    end
    
    % Handle empty ref_header
    if isempty(ref_header)
        ref_header = create_minimal_header(image_data);
    end
    
    % Create output directory if needed
    output_dir = fileparts(output_path);
    if create_dir && ~isempty(output_dir) && ~exist(output_dir, 'dir')
        mkdir(output_dir);
        if verbose
            fprintf('Created directory: %s\n', output_dir);
        end
    end
    
    % Handle complex data by saving magnitude and phase separately
    if ~isreal(image_data)
        if verbose
            fprintf('Complex data detected - saving magnitude and phase components separately\n');
        end
        
        % Get base path without extension
        [filepath, name, ext] = fileparts(output_path);
        
        % Save magnitude component
        mag_data = abs(image_data);
        mag_path = fullfile(filepath, [name '_mag' ext]);
        save_component(mag_data, mag_path, ref_header, data_type, description, verbose, 'magnitude');
        
        % Save phase component
        phase_data = angle(image_data);
        phase_path = fullfile(filepath, [name '_phase' ext]);
        save_component(phase_data, phase_path, ref_header, data_type, description, verbose, 'phase');
        
        if verbose
            fprintf('Successfully saved complex data as:\n');
            fprintf('  Magnitude: %s\n', mag_path);
            fprintf('  Phase: %s\n', phase_path);
        end
        
    else
        % Save real data directly
        save_component(image_data, output_path, ref_header, data_type, description, verbose, 'real');
        
        if verbose
            fprintf('Successfully saved real data: %s\n', output_path);
        end
    end
end

function minimal_header = create_minimal_header(image_data)
%CREATE_MINIMAL_HEADER Create minimal header structure
    minimal_header = struct();
    minimal_header.ImageSize = size(image_data);
    minimal_header.PixelDimensions = ones(1, length(size(image_data)));
    minimal_header.Datatype = 'single';
    minimal_header.BitsPerPixel = 32;
end

function save_component(image_data, output_path, ref_header, data_type, description, verbose, component_type)
%SAVE_COMPONENT Save a single component (magnitude, phase, or real data)
    
    % Convert to specified data type without scaling
    switch lower(data_type)
        case 'single'
            image_data = single(image_data);
            datatype_code = 16;
            bits_per_pixel = 32;
        case 'double'
            image_data = double(image_data);
            datatype_code = 64;
            bits_per_pixel = 64;
        case 'int16'
            image_data = int16(image_data);
            datatype_code = 4;
            bits_per_pixel = 16;
        case 'uint8'
            image_data = uint8(image_data);
            datatype_code = 2;
            bits_per_pixel = 8;
        otherwise
            error('Unsupported data_type: %s', data_type);
    end
    
    if verbose
        data_range = [min(image_data(:)), max(image_data(:))];
        fprintf('Processing %s component - range: [%.3f, %.3f] (%s)\n', ...
                component_type, data_range(1), data_range(2), data_type);
    end
    
    % Prepare output header and save using fallback methods
    save_with_fallback(image_data, output_path, ref_header, datatype_code, ...
                      bits_per_pixel, description, component_type, verbose);
end

function output_header = prepare_matlab_header(ref_header, image_data, datatype_code, bits_per_pixel, description)
%PREPARE_MATLAB_HEADER Create header structure compatible with MATLAB niftiwrite
    
    % Start with a clean header structure for MATLAB
    output_header = struct();
    
    % Copy essential fields from reference header if they exist
    essential_fields = {'PixelDimensions', 'ImageSize', 'Datatype', 'BitsPerPixel', ...
                       'SpaceUnits', 'TimeUnits', 'Transform', 'DisplayIntensityRange'};
    
    for i = 1:length(essential_fields)
        field = essential_fields{i};
        if isfield(ref_header, field)
            output_header.(field) = ref_header.(field);
        end
    end
    
    % Override with new values
    output_header.ImageSize = size(image_data);
    output_header.Datatype = get_matlab_datatype_string(datatype_code);
    output_header.BitsPerPixel = bits_per_pixel;
    
    % Remove any fields that might cause issues with MATLAB's niftiwrite
    problematic_fields = {'Description', 'desc', 'descrip'};
    for i = 1:length(problematic_fields)
        if isfield(output_header, problematic_fields{i})
            output_header = rmfield(output_header, problematic_fields{i});
        end
    end
    
    % Ensure we have basic spatial information
    if ~isfield(output_header, 'PixelDimensions')
        output_header.PixelDimensions = ones(1, length(size(image_data)));
    end
    
    % Create basic transform if missing
    if ~isfield(output_header, 'Transform')
        ndims_data = length(size(image_data));
        if ndims_data >= 3
            % Create basic affine transform
            transform_matrix = eye(4);
            if isfield(output_header, 'PixelDimensions')
                voxel_size = output_header.PixelDimensions;
                for dim = 1:min(3, length(voxel_size))
                    transform_matrix(dim, dim) = voxel_size(dim);
                end
            end
            output_header.Transform = struct('T', transform_matrix, 'Dimensionality', 3);
        end
    end
    
    % Set data range
    data_min = min(image_data(:));
    data_max = max(image_data(:));
    output_header.DisplayIntensityRange = [data_min, data_max];
end

function output_header = prepare_minimal_matlab_header(image_data, datatype_code)
%PREPARE_MINIMAL_MATLAB_HEADER Create minimal header for MATLAB niftiwrite
    output_header = struct();
    output_header.ImageSize = size(image_data);
    output_header.Datatype = get_matlab_datatype_string(datatype_code);
    output_header.PixelDimensions = ones(1, length(size(image_data)));
    
    % Only add transform for 3D+ data
    if length(size(image_data)) >= 3
        output_header.Transform = struct('T', eye(4), 'Dimensionality', 3);
    end
end

function datatype_str = get_matlab_datatype_string(datatype_code)
%GET_MATLAB_DATATYPE_STRING Convert numeric datatype code to MATLAB string
    switch datatype_code
        case 2
            datatype_str = 'uint8';
        case 4
            datatype_str = 'int16';
        case 8
            datatype_str = 'int32';
        case 16
            datatype_str = 'single';
        case 64
            datatype_str = 'double';
        case 256
            datatype_str = 'int8';
        case 512
            datatype_str = 'uint16';
        case 768
            datatype_str = 'uint32';
        otherwise
            datatype_str = 'single';  % Default fallback
    end
end

function save_with_fallback(image_data, output_path, ref_header, datatype_code, bits_per_pixel, description, component_type, verbose)
%SAVE_WITH_FALLBACK Try multiple NIfTI saving methods with robust fallbacks
    
    success = false;
    methods_tried = {};
    last_error = '';
    
    % Method 1: Try MATLAB built-in with full header
    if ~success && exist('niftiwrite', 'file') == 2
        try
            matlab_header = prepare_matlab_header(ref_header, image_data, datatype_code, bits_per_pixel, description);
            niftiwrite(image_data, output_path, matlab_header);
            success = true;
            methods_tried{end+1} = 'MATLAB built-in (full header)';
            if verbose
                fprintf('    Saved using MATLAB built-in niftiwrite (full header)\n');
            end
        catch ME
            last_error = ME.message;
            methods_tried{end+1} = sprintf('MATLAB built-in full (failed: %s)', ME.message);
        end
    end
    
    % Method 2: Try MATLAB built-in with minimal header
    if ~success && exist('niftiwrite', 'file') == 2
        try
            minimal_header = prepare_minimal_matlab_header(image_data, datatype_code);
            niftiwrite(image_data, output_path, minimal_header);
            success = true;
            methods_tried{end+1} = 'MATLAB built-in (minimal header)';
            if verbose
                fprintf('    Saved using MATLAB built-in niftiwrite (minimal header)\n');
            end
        catch ME
            last_error = ME.message;
            methods_tried{end+1} = sprintf('MATLAB built-in minimal (failed: %s)', ME.message);
        end
    end
    
    % Method 3: Try MATLAB built-in with no extra header info
    if ~success && exist('niftiwrite', 'file') == 2
        try
            niftiwrite(image_data, output_path);
            success = true;
            methods_tried{end+1} = 'MATLAB built-in (no header)';
            if verbose
                fprintf('    Saved using MATLAB built-in niftiwrite (no header)\n');
            end
        catch ME
            last_error = ME.message;
            methods_tried{end+1} = sprintf('MATLAB built-in no header (failed: %s)', ME.message);
        end
    end
    
    % Method 4: Try NIfTI toolbox
    if ~success && exist('save_nii', 'file') == 2
        try
            nii_struct = create_nii_struct(image_data, ref_header, datatype_code, bits_per_pixel, description, component_type);
            save_nii(nii_struct, output_path);
            success = true;
            methods_tried{end+1} = 'NIfTI toolbox';
            if verbose
                fprintf('    Saved using NIfTI toolbox\n');
            end
        catch ME
            last_error = ME.message;
            methods_tried{end+1} = sprintf('NIfTI toolbox (failed: %s)', ME.message);
        end
    end
    
    % Method 5: Try SPM
    if ~success && exist('spm_write_vol', 'file') == 2
        try
            V = create_spm_vol_struct(image_data, ref_header, datatype_code, output_path, description);
            spm_write_vol(V, image_data);
            success = true;
            methods_tried{end+1} = 'SPM';
            if verbose
                fprintf('    Saved using SPM\n');
            end
        catch ME
            last_error = ME.message;
            methods_tried{end+1} = sprintf('SPM (failed: %s)', ME.message);
        end
    end
    
    % Method 6: Try manual NIfTI creation (basic header)
    if ~success
        try
            write_basic_nifti(image_data, output_path, datatype_code);
            success = true;
            methods_tried{end+1} = 'Manual NIfTI (basic)';
            if verbose
                fprintf('    Saved using manual NIfTI writer (basic header)\n');
            end
        catch ME
            last_error = ME.message;
            methods_tried{end+1} = sprintf('Manual NIfTI basic (failed: %s)', ME.message);
        end
    end
    
    % Method 7: Try raw binary save as last resort
    if ~success
        try
            [filepath, name, ~] = fileparts(output_path);
            raw_path = fullfile(filepath, [name '.raw']);
            
            % Save as raw binary
            fid = fopen(raw_path, 'w');
            if fid == -1
                error('Cannot open file for writing');
            end
            fwrite(fid, image_data, class(image_data));
            fclose(fid);
            
            % Create a simple info file
            info_path = fullfile(filepath, [name '_info.txt']);
            fid = fopen(info_path, 'w');
            if fid ~= -1
                fprintf(fid, 'Raw binary data file: %s\n', [name '.raw']);
                fprintf(fid, 'Data type: %s\n', class(image_data));
                fprintf(fid, 'Dimensions: %s\n', mat2str(size(image_data)));
                fprintf(fid, 'Data range: [%.6f, %.6f]\n', min(image_data(:)), max(image_data(:)));
                fclose(fid);
            end
            
            success = true;
            methods_tried{end+1} = 'Raw binary fallback';
            if verbose
                fprintf('    Saved as raw binary: %s (with info file: %s)\n', raw_path, info_path);
                warning('Could not save as NIfTI format. Saved as raw binary with info file instead.');
            end
        catch ME
            last_error = ME.message;
            methods_tried{end+1} = sprintf('Raw binary (failed: %s)', ME.message);
        end
    end
    
    if ~success
        error(['Failed to save file using any available method.\n' ...
               'Methods tried: %s\n' ...
               'Last error: %s\n' ...
               'Please check file permissions and available disk space.'], ...
               strjoin(methods_tried, ', '), last_error);
    end
end

function write_basic_nifti(image_data, output_path, datatype_code)
%WRITE_BASIC_NIFTI Write a basic NIfTI file manually
    
    % Create minimal NIfTI header (348 bytes)
    header = zeros(348, 1, 'uint8');
    
    % Header size
    header(1:4) = typecast(int32(348), 'uint8');
    
    % Magic string "n+1" at end
    header(345:348) = uint8(['n+1' char(0)]);
    
    % Data dimensions
    img_dims = size(image_data);
    ndims = length(img_dims);
    
    % dim[0] = number of dimensions
    header(41:42) = typecast(int16(ndims), 'uint8');
    
    % dim[1:8] = dimensions
    for i = 1:min(8, ndims)
        start_idx = 43 + (i-1)*2;
        header(start_idx:start_idx+1) = typecast(int16(img_dims(i)), 'uint8');
    end
    
    % Datatype and bits per pixel
    header(71:72) = typecast(int16(datatype_code), 'uint8');
    
    switch datatype_code
        case 2  % uint8
            bitpix = 8;
        case 4  % int16
            bitpix = 16;
        case 8  % int32
            bitpix = 32;
        case 16 % single
            bitpix = 32;
        case 64 % double
            bitpix = 64;
        otherwise
            bitpix = 32;
    end
    header(73:74) = typecast(int16(bitpix), 'uint8');
    
    % Voxel offset (where data starts)
    header(109:112) = typecast(single(352), 'uint8');
    
    % Scale slope and intercept
    header(113:116) = typecast(single(1), 'uint8'); % slope
    header(117:120) = typecast(single(0), 'uint8'); % intercept
    
    % Write header and data
    fid = fopen(output_path, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', output_path);
    end
    
    % Write header
    fwrite(fid, header, 'uint8');
    
    % Write padding to reach offset 352
    padding = zeros(4, 1, 'uint8');
    fwrite(fid, padding, 'uint8');
    
    % Write data
    fwrite(fid, image_data, class(image_data));
    
    fclose(fid);
end

function nii_struct = create_nii_struct(image_data, ref_header, datatype_code, bits_per_pixel, description, component_type)
%CREATE_NII_STRUCT Create structure compatible with NIfTI toolbox save_nii
    
    nii_struct = struct();
    nii_struct.img = image_data;
    
    % Create header structure compatible with NIfTI toolbox
    nii_struct.hdr = struct();
    nii_struct.hdr.hk = struct();  % Header key
    nii_struct.hdr.dime = struct(); % Data dimension
    nii_struct.hdr.hist = struct(); % History
    
    % Header key section
    nii_struct.hdr.hk.sizeof_hdr = 348;
    nii_struct.hdr.hk.data_type = '';
    nii_struct.hdr.hk.db_name = '';
    nii_struct.hdr.hk.extents = 0;
    nii_struct.hdr.hk.session_error = 0;
    nii_struct.hdr.hk.regular = 'r';
    nii_struct.hdr.hk.dim_info = 0;
    
    % Data dimensions section
    img_size = size(image_data);
    nii_struct.hdr.dime.dim = [length(img_size), img_size, ones(1, 8-length(img_size))];
    nii_struct.hdr.dime.intent_p1 = 0;
    nii_struct.hdr.dime.intent_p2 = 0;
    nii_struct.hdr.dime.intent_p3 = 0;
    nii_struct.hdr.dime.intent_code = 0;
    nii_struct.hdr.dime.datatype = datatype_code;
    nii_struct.hdr.dime.bitpix = bits_per_pixel;
    nii_struct.hdr.dime.slice_start = 0;
    
    % Pixel dimensions - fix the voxel_size variable scope issue
    if isfield(ref_header, 'PixelDimensions')
        voxel_size = ref_header.PixelDimensions;
        pixdim = [1, voxel_size, ones(1, 8-length(voxel_size))];
    else
        voxel_size = ones(1, length(img_size));
        pixdim = [1, voxel_size, ones(1, 8-length(img_size))];
    end
    nii_struct.hdr.dime.pixdim = pixdim;
    
    nii_struct.hdr.dime.vox_offset = 352;
    nii_struct.hdr.dime.scl_slope = 1;
    nii_struct.hdr.dime.scl_inter = 0;
    nii_struct.hdr.dime.slice_end = 0;
    nii_struct.hdr.dime.slice_code = 0;
    nii_struct.hdr.dime.xyzt_units = 10;  % mm and seconds
    nii_struct.hdr.dime.cal_max = max(image_data(:));
    nii_struct.hdr.dime.cal_min = min(image_data(:));
    nii_struct.hdr.dime.slice_duration = 0;
    nii_struct.hdr.dime.toffset = 0;
    nii_struct.hdr.dime.glmax = max(image_data(:));
    nii_struct.hdr.dime.glmin = min(image_data(:));
    
    % History section
    if ~isempty(description)
        nii_struct.hdr.hist.descrip = description;
    else
        nii_struct.hdr.hist.descrip = sprintf('Processed %s component', component_type);
    end
    
    nii_struct.hdr.hist.aux_file = '';
    nii_struct.hdr.hist.qform_code = 1;
    nii_struct.hdr.hist.sform_code = 1;
    nii_struct.hdr.hist.quatern_b = 0;
    nii_struct.hdr.hist.quatern_c = 0;
    nii_struct.hdr.hist.quatern_d = 0;
    nii_struct.hdr.hist.qoffset_x = 0;
    nii_struct.hdr.hist.qoffset_y = 0;
    nii_struct.hdr.hist.qoffset_z = 0;
    
    % Use the voxel_size variable that's now in scope
    if length(voxel_size) >= 3
        nii_struct.hdr.hist.srow_x = [voxel_size(1), 0, 0, 0];
        nii_struct.hdr.hist.srow_y = [0, voxel_size(2), 0, 0];
        nii_struct.hdr.hist.srow_z = [0, 0, voxel_size(3), 0];
    else
        nii_struct.hdr.hist.srow_x = [1, 0, 0, 0];
        nii_struct.hdr.hist.srow_y = [0, 1, 0, 0];
        nii_struct.hdr.hist.srow_z = [0, 0, 1, 0];
    end
    
    nii_struct.hdr.hist.intent_name = '';
    nii_struct.hdr.hist.magic = 'n+1';
    
    % Copy spatial transformation if available
    if isfield(ref_header, 'Transform') && isfield(ref_header.Transform, 'T')
        transform = ref_header.Transform.T;
        if size(transform, 1) >= 4 && size(transform, 2) >= 4
            nii_struct.hdr.hist.srow_x = transform(1, :);
            nii_struct.hdr.hist.srow_y = transform(2, :);
            nii_struct.hdr.hist.srow_z = transform(3, :);
        end
    end
end

function V = create_spm_vol_struct(image_data, ref_header, datatype_code, output_path, description)
%CREATE_SPM_VOL_STRUCT Create volume structure compatible with SPM
    
    V = struct();
    V.fname = output_path;
    V.dim = size(image_data);
    V.dt = [datatype_code, 0];
    
    % Spatial transformation matrix
    if isfield(ref_header, 'Transform') && isfield(ref_header.Transform, 'T')
        V.mat = ref_header.Transform.T;
    elseif isfield(ref_header, 'PixelDimensions')
        % Create basic transformation matrix from voxel size
        voxel_size = ref_header.PixelDimensions;
        if length(voxel_size) >= 3
            V.mat = diag([voxel_size(1:3), 1]);
        else
            V.mat = eye(4);
        end
    else
        V.mat = eye(4);  % Identity matrix as fallback
    end
    
    % Description
    if ~isempty(description)
        V.descrip = description;
    else
        V.descrip = 'Processed image';
    end
    
    % Additional SPM fields
    V.n = [1, 1];
    V.pinfo = [1; 0; 0];
end