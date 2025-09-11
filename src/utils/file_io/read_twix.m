function data = read_twix(twix_filename, output_path, varargin)
%READ_TWIX  Convert Siemens twix file to NIfTI format, with caching.
%
%   data = READ_TWIX(twix_filename, output_path, ...)
%
%   If the NIfTI file at output_path exists, it loads the data from it and returns.
%   Otherwise, reads the Siemens twix file (.dat) using mapVBVD, processes it,
%   and saves it as a NIfTI file at output_path, including header info.


    % Parse inputs
    p = inputParser;
    addRequired(p, 'twix_filename', @ischar);
    addRequired(p, 'output_path', @ischar);
    addParameter(p,'modality', 'image', @ischar)
    addParameter(p, 'data_struct', @isstruct);
    parse(p, twix_filename, output_path, varargin{:});

    data_struct = p.Results.data_struct;
    modality = p.Results.modality;
    
    % Path to the folder containing mapVBVD
    mapvbvd_folder = '../../scripts/'; % <-- CHANGE THIS
    addpath(mapvbvd_folder);

    if ~exist('mapVBVD', 'file')
        error('mapVBVD.m not found. Please check the mapvbvd_folder path.');
    end
    

    if strcmp(modality,'image')
        % If output file exists, load and return
        magnitude_path = replace(output_path, '.nii', '_magnitude.nii');
        phase_path = replace(output_path, '.nii', '_phase.nii');
        if exist(magnitude_path, 'file') && exist(phase_path, 'file')
            fprintf('Loading existing NIfTI from: %s\n', output_path);
            mag = niftiread(magnitude_path);
            phase = niftiread(phase_path);
            data = mag .* exp(1i * phase);
            return;
        end
    
    
    
        % Check if twix file exists
        if ~exist(twix_filename, 'file')
            error('Twix file not found: %s', twix_filename);
        end
    
        fprintf('Reading twix file: %s\n', twix_filename);
    
        % Read twix file
        try
            twix_data = mapVBVD(twix_filename);
        catch ME
            error('Failed to read twix file with mapVBVD: %s', ME.message);
        end
    
        % If multiple measurements, take the last one
        if iscell(twix_data)
            if length(twix_data) > 1
                fprintf('Multiple measurements found. Using the last one.\n');
            end
            twix_data = twix_data{end};
        end
    
        % Extract image data
        if isfield(twix_data, 'image')
            data = twix_data.image;
            data.flagRemoveOS = true;
            fprintf('Using image data from twix file.\n');
        elseif isfield(twix_data, 'data')
            data = twix_data.data();
            fprintf('Using raw data from twix file.\n');
        else
            error('No image or data field found in twix structure');
        end
    
        % Display original dims
        original_size = data.dataSize;
        dim_names = {'Col', 'Cha', 'Lin', 'Par', 'Sli', 'Ave', 'Phs', 'Eco', ...
                     'Rep', 'Set', 'Seg', 'Ida', 'Idb', 'Idc', 'Idd', 'Ide'};
        fprintf('Original data dimensions:\n');
        for i = 1:length(original_size)
            if i <= length(dim_names)
                fprintf('  %s: %d\n', dim_names{i}, original_size(i));
            else
                fprintf('  Dim%d: %d\n', i, original_size(i));
            end
        end
    
        data = squeeze(data);
        data = permute(squeeze(data(:,:,:,:,1)), [1,3,4,2]);
    
        final_size = size(data);
        fprintf('Final data dimensions after squeezing: %s\n', mat2str(final_size));
    
        % % ---- Extract NIfTI header info from twix_data ----
        % hdr = niftiinfoFromTwix(twix_data, final_size);
    
        % Handle complex data - save magnitude and phase separately
        if ~isreal(data)
            % Extract magnitude and phase
            magnitude_data = abs(data);
            phase_data = angle(data);
            
            % Extract headder info from the MPMs
            ref_header = niftiinfo(data_struct.R1.filepath);
            hdr = extractHeaderFromData(ref_header, magnitude_data);
            
            % Write magnitude NIfTI file with minimal header (let MATLAB create defaults)
            niftiwrite(single(magnitude_data), magnitude_path,hdr);
            
            % Write phase NIfTI file with minimal header (let MATLAB create defaults)
            niftiwrite(single(phase_data), phase_path, hdr);
            
            fprintf('Successfully converted CFL to NIfTI:\n');
            fprintf('  Magnitude: %s\n', magnitude_path);
            fprintf('  Phase: %s\n', phase_path);
            
        else
            % Data is real, save as single file
            output_path = fullfile(output_dir, [filename '.nii']);
    
            % Extract headder info from the twix
            ref_header = niftiinfo(data_struct.R1.filepath);
            hdr = extractHeaderFromData(ref_header, data);
    
            
            % Write NIfTI file with minimal header (let MATLAB create defaults)
            niftiwrite(single(data),output_path, hdr);
            
            fprintf('Successfully converted CFL to NIfTI: %s\n', output_path);
        end
        
    elseif strcmp(modality,'noise')
        if exist(output_path, 'file') 
            fprintf('Loading existing noise scan from: %s\n', output_path);
            data = load(output_path);
            return;
        end
        
         % Check if twix file exists
        if ~exist(twix_filename, 'file')
            error('Twix file not found: %s', twix_filename);
        end
    
        fprintf('Reading twix file: %s\n', twix_filename);
    
        % Read twix file
        try
            twix_data = mapVBVD(twix_filename);
        catch ME
            error('Failed to read twix file with mapVBVD: %s', ME.message);
        end
    
        % If multiple measurements, take the last one
        if iscell(twix_data)
            if length(twix_data) > 1
                fprintf('Multiple measurements found. Using the last one.\n');
            end
            twix_data = twix_data{end};
        end
    
        % Extract noise data
        isfield(twix_data, 'image')
        data = twix_data.noise;
        data.flagRemoveOS = true;
        fprintf('Using noise data from twix file.\n');
        data = data(:,:);
        save(output_path, 'data');

   end



end
function hdr_out = extractHeaderFromData(hdr_in, data)
% extractHeaderFromData  Header updater for 4D k-space data (multi-coil MRI)
%
% hdr_out = extractHeaderFromData(hdr_in, data)
%
% Key modifications for 4D k-space:
% 1. Forces 4D ImageSize while preserving spatial geometry (qform/sform)
% 2. Sets 4th PixelDimension to 1 (unitary for coil dimension)
% 3. Handles complex datatype for k-space
% 4. Adds coil count to description
% 5. Maintains NIfTI-1 compatibility

    if nargin<2, error('Provide hdr_in and data'); end
    if ~isstruct(hdr_in), error('hdr_in must be a struct (niftiinfo-like)'); end

    hdr_out = hdr_in; % shallow copy
    
    % --- Get data dimensions ---
    data_sz = size(data);
    ndims_data = numel(data_sz);
    data_sz = double(data_sz(:)'); % row vector
    
    % --- Force ImageSize to match data ---
    hdr_out.ImageSize = data_sz;
    
    % --- Handle PixelDimensions ---
    if isfield(hdr_out, 'PixelDimensions')
        % Preserve spatial resolution (first 3 dims)
        spatial_pd = hdr_out.PixelDimensions(1:min(3,end));
        
        % Initialize with default 1.0 for higher dimensions
        new_pd = ones(1, ndims_data);
        new_pd(1:length(spatial_pd)) = spatial_pd;
        
        hdr_out.PixelDimensions = new_pd;
    else
        hdr_out.PixelDimensions = ones(1, ndims_data);
    end
    
    % --- Update NIfTI raw header fields ---
    if isfield(hdr_out, 'raw')
        % dim[0] = number of dimensions (4 for 4D)
        hdr_out.raw.dim = [ndims_data, data_sz, ones(1, 8-ndims_data)];
        
        % pixdim[0] unused, then resolutions (matches PixelDimensions)
        hdr_out.raw.pixdim = [0, hdr_out.PixelDimensions];
        
        % Handle complex k-space data
        if ~isreal(data)
            hdr_out.raw.datatype = 32;  % NIfTI complex float32
            hdr_out.raw.bitpix = 64;    % 32 real + 32 imag
        end
    end
    
    % --- Preserve spatial transforms (CRITICAL) ---
    % Ensure qform/sform matrices remain unchanged
    if isfield(hdr_out, 'Transform')
        % Do not modify - keeps spatial alignment
    elseif isfield(hdr_out, 'raw')
        if isfield(hdr_out.raw, 'qform_code') && hdr_out.raw.qform_code > 0
            % qform exists - preserve it
        elseif isfield(hdr_out.raw, 'sform_code') && hdr_out.raw.sform_code > 0
            % sform exists - preserve it
        end
    end
    
    % --- Update description with coil info ---
    coil_count = 1;
    if ndims_data >= 4
        coil_count = data_sz(4);
        desc_suffix = sprintf(' | %d-channel k-space', coil_count);
        
        % Update all possible description fields (safely)
        desc_fields = {'Description', 'descrip', 'Raw.descrip', 'raw.descrip'};
        for f = 1:length(desc_fields)
            try
                parts = strsplit(desc_fields{f}, '.');
                if numel(parts) == 1
                    % Top-level field
                    old_desc = hdr_out.(parts{1});
                    new_desc = [old_desc desc_suffix];
                    hdr_out.(parts{1}) = new_desc(1:min(79,end)); % Truncate
                else
                    % Nested field (e.g., Raw.descrip)
                    old_desc = hdr_out.(parts{1}).(parts{2});
                    new_desc = [old_desc desc_suffix];
                    hdr_out.(parts{1}).(parts{2}) = new_desc(1:min(79,end));
                end
                break; % Stop after first successful update
            catch
                continue;
            end
        end
    end
    
    % --- Ensure valid datatype ---
    if isfield(hdr_out, 'Datatype')
        if ~isreal(data)
            hdr_out.Datatype = 'single'; % Complex stored as single
        else
            hdr_out.Datatype = class(data);
        end
    end
end