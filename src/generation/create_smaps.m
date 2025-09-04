function  create_smaps(twix_filename, sens_maps_out_file, varargin)
% CREATE_SMAPS  Compute ESPIRiT sensitivity maps (simplified, preserving original logic)
%
%   smaps = create_smaps(twix_filename, sens_maps_out_file, ...)
%
% Required:
%   twix_filename         - path to Siemens TWIX (.dat)
%   sens_maps_out_file    - output filename for sensitivity maps (NIfTI or .mat)
%
% Optional name/value pairs (defaults match your original script):
%   'espirit_path'  : '/data/u_kilicb_software/ESPIRiT/'
%   'ncalib'        : 40
%   'ksize'         : [6 6]
%   'eigThresh1'    : 0.02
%   'eigThresh2'    : 0.95
%   'show_figures'  : 'none'    % options: 'none','essential','detailed'
%   'verbose'       : true
%
% Notes:
% - Assumes an existing read_twix function to read data into DATA.
% - DATA is expected 4D: [n_slice, nx, ny, n_channel]
% - Final smaps returned and saved with shape: [nx, n_slice, ny, n_channel]

% -------------------- parse inputs --------------------
p = inputParser;
addRequired(p,'twix_filename', @(x)ischar(x)||isstring(x));
addRequired(p,'sens_maps_out_file', @(x)ischar(x)||isstring(x));
addParameter(p,'espirit_path','/data/u_kilicb_software/ESPIRiT/', @(x)ischar(x)||isstring(x));
addParameter(p,'ncalib',20,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'ksize',[5 5],@(x)isnumeric(x)&&numel(x)==2);
addParameter(p,'eigThresh1',0.1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'eigThresh2',0.98,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'show_figures','none', @(x)any(validatestring(x,{'none','essential','detailed'})));
addParameter(p,'verbose',true,@islogical);
addParameter(p, 'data_struct', @(x)isstruct(x));
addParameter(p,'smoothing_radius',25,@(x)isnumeric(x)&&isscalar(x));
parse(p, twix_filename, sens_maps_out_file, varargin{:});

espirit_path = char(p.Results.espirit_path);
ncalib       = p.Results.ncalib;
ksize        = p.Results.ksize;
eigThresh1   = p.Results.eigThresh1;
eigThresh2   = p.Results.eigThresh2;
show_figures = validatestring(p.Results.show_figures,{'none','essential','detailed'});
verbose      = p.Results.verbose;
data_struct =  p.Results.data_struct;
smoothing_radius = p.Results.smoothing_radius;

% -------------------- simple checks --------------------
if ~exist(twix_filename,'file')
    error('TWIX file not found: %s', twix_filename);
end

% If output exists, load and return quickly
if exist(sens_maps_out_file,'file')
    if verbose, fprintf('Output exists. Loading %s and returning.\n', sens_maps_out_file); end
    try
        smaps = niftiread(sens_maps_out_file);
        return;
    catch
        % fallback: if .mat file
        try
            tmp = load(sens_maps_out_file);
            if isfield(tmp,'smaps'), smaps = tmp.smaps; return; end
        catch
            warning('Could not load existing file (%s). Will recompute.', sens_maps_out_file);
        end
    end
end

% -------------------- add ESPIRiT path and run setPath.m (if available) ----
if exist(espirit_path,'dir')
    addpath(genpath(espirit_path));
    sp = fullfile(espirit_path,'setPath.m');
    if exist(sp,'file')
        try
            run(sp);
            if verbose, fprintf('Ran setPath.m in %s\n', espirit_path); end
        catch
            warning('setPath.m found but could not be run. Continuing anyway.');
        end
    else
        if verbose, warning('setPath.m not found in ESPIRiT path.'); end
    end
else
    if verbose, warning('ESPIRiT path not found: %s. Ensure dat2Kernel/kernelEig are on path.', espirit_path); end
end

% -------------------- check helper functions --------------------
if exist('read_twix','file')~=2
    error('read_twix function not found on path. Please provide it.');
end
if exist('dat2Kernel','file')~=2
    error('dat2Kernel not found on path. Please ensure ESPIRiT utilities are available.');
end
if exist('kernelEig','file')~=2
    error('kernelEig not found on path. Please ensure ESPIRiT utilities are available.');
end
% ESPIRiT constructor is optional; we'll warn if missing.
have_ESPIRiT = (exist('ESPIRiT','file')==2);
if ~(exist('ESPIRiT','file')==2) && verbose
    warning('ESPIRiT constructor not found. Using fallback projection if needed.');
end

% -------------------- read data using existing read_twix --------------------

twix_file_nii = replace(sens_maps_out_file, 'smap', 'kspace_raw');
magnitude_path = replace(twix_file_nii, '.nii.gz', '_magnitude.nii.gz');
phase_path = replace(twix_file_nii, '.nii.gz', '_phase.nii.gz');
if exist(magnitude_path, 'file') && exist(phase_path, 'file')
    fprintf('Reading ')
    data_mag = niftiread(magnitude_path); 
    data_phase = niftiread(phase_path);
    DATA = data_mag .* exp(1i* data_phase);
else
    try
        % assume read_twix returns numeric DATA or a struct with .data/.image; preserve original naming
        if verbose, fprintf('Reading TWIX data using read_twix: %s\n', twix_filename); end
        DATA = read_twix(twix_filename, twix_file_nii, 'data_struct', data_struct);
    catch ME
        error('read_twix failed: %s', ME.message);
    end
end

% % -------------------- optional: extract header snippet via mapVBVD ----
twix_hdr_snip = '';
% if exist('mapVBVD','file')==2
%     try
%         t = mapVBVD(twix_filename);
%         if iscell(t), t = t{end}; end
%         if isfield(t,'hdr')
%             % stringify small part
%             try
%                 twix_hdr_snip = jsonencode(t.hdr);
%             catch
%                 twix_hdr_snip = 'TWIX header present (stringify failed)';
%             end
%         end
%     catch
%         % ignore, keep empty snippet
%     end
% end

% -------------------- data shape expectation & small validation -------
% Expect DATA to be 4D: [n_slice, nx, ny, n_channel]
if ndims(DATA) ~= 4
    error('Expected DATA to be 4D [n_slice, nx, ny, n_channel], but got shape: %s', mat2str(size(DATA)));
end
[n_slice, nx, ny, n_channel] = size(DATA);
if verbose
    fprintf('DATA shape: n_slice=%d, nx=%d, ny=%d, n_channel=%d\n', n_slice, nx, ny, n_channel);
end

% -------------------- ifft in readout direction (original code did ifftc(DATA,1)) ----
if exist('ifftc','file')==2
    kspace_ifft_readout = ifftc(DATA, 1);  % keep original convention
else
    % fallback: use fftshifted ifft along dim1
    if verbose, warning('ifftc not found. Using ifft along dim 1 with fftshift fallback.'); end
    kspace_ifft_readout = ifft(ifftshift(DATA,1), [], 1); % simple fallback (not identical to ifftc)
end

% -------------------- allocate smaps in original orientation ----------------
% Will store as smaps(slice, x, y, coil) as in original then permute at end.
smaps = zeros(n_slice, nx, ny, n_channel, 'like', DATA);




% Progress tracking setup
progress_interval = max(1, floor(n_slice / 20)); % Print progress every ~5%
processed_slices = 0;
start_time = tic;

fprintf('Processing %d slices...\n', n_slice);
fprintf('Progress: [');
for i = 1:20, fprintf(' '); end
fprintf('] 0%%\n');




% -------------------- per-slice processing (close to original code) ----------
for slice_idx = 1:n_slice
    try
        kspace_slice = squeeze(kspace_ifft_readout(slice_idx,:,:,:));  % [nx, ny, nc]
        [sx, sy, Nc] = size(kspace_slice); 
        
        % calib crop (center)
        try
            calib = crop(kspace_slice, [ncalib, ncalib, Nc]);
        catch
            % very simple fallback central crop if crop missing
            cx = round(sx/2); cy = round(sy/2);
            xr = max(1,cx-floor(ncalib/2)) : min(sx, cx+ceil(ncalib/2)-1);
            yr = max(1,cy-floor(ncalib/2)) : min(sy, cy+ceil(ncalib/2)-1);
            calib = kspace_slice(xr, yr, :);
            if verbose, warning('crop failed for slice %d; used central crop.', slice_idx); end
        end
        
        % display coil magnitude if essential/detailed
        
        if ~strcmp(show_figures,'none') && (strcmp(show_figures,'essential') || strcmp(show_figures,'detailed'))
            try
                im = ifft2c(kspace_slice);  % coil images
                figure('Name',sprintf('Slice %d - coil magnitude',slice_idx));
                imshow3(abs(im),[],[4,8]);
                title(sprintf('Slice %d - coil magnitude', slice_idx));
            catch
                if verbose, warning('Could not display coil magnitude for slice %d', slice_idx); end
            end
        end
        
        % dat2Kernel and first SVD etc.
        [k, S] = dat2Kernel(calib, ksize);
        idx = find(S >= S(1)*eigThresh1, 1, 'last');
        if isempty(idx), idx = size(k,4); end
        
        % optional detailed plots: singular values & kernel display
        if strcmp(show_figures,'detailed')
            try
                kdisp = reshape(k, [ksize(1)*ksize(2)*Nc, ksize(1)*ksize(2)*Nc]);
                figure('Name',sprintf('Slice %d - singular values',slice_idx));
                subplot(2,1,1); plot(S,'LineWidth',2); hold on;
                plot([1,numel(S)],[S(1)*eigThresh1 S(1)*eigThresh1],'r-','LineWidth',2);
                title('Singular values');
                subplot(2,1,2); imagesc(abs(kdisp)); colormap(gray(256)); title('Kernel vectors');
            catch
                if verbose, warning('Could not display kernel diagnostics for slice %d', slice_idx); end
            end
        end
        
        % kernelEig -> M (eigenvectors), W (eigenvalues)
        [M, W] = kernelEig(k(:,:,:,1:idx), [sx, sy]);
        
        % detailed: show eigenvalue image and eigenvector magnitude
        if strcmp(show_figures,'detailed')
            try
                figure('Name',sprintf('Slice %d - eigenvalues',slice_idx));
                imshow3(abs(W),[],[4,8]); title('Eigenvalues (magnitude)');
                figure('Name',sprintf('Slice %d - eigenvector magnitudes',slice_idx));
                imshow3(abs(M),[],[32,32]); title('Eigenvector magnitudes');
            catch
                if verbose, warning('Could not show eigen maps for slice %d', slice_idx); end
            end
        end
        
        % % projection (ESPIRiT if available else fallback)
        % if (exist('ESPIRiT','file')==2)
        %     try
        %         ESP = ESPIRiT(M);
        %         P = ESP' * im;
        %     catch
        %         P = sum(repmat(im,[1,1,1,size(M,4)]) .* conj(M), 3);
        %     end
        % else
        %     P = sum(repmat(im,[1,1,1,size(M,4)]) .* conj(M), 3);
        % end
        
        % select final maps: last eigenvector * mask of eigenvalue>eigThresh2
        mask = W(:,:,end) > eigThresh2;
        maps = M(:,:,:,end) .* repmat(mask, [1,1,Nc]);
        
        % show final map magnitudes (essential/detailed)
        if ~strcmp(show_figures,'none') && (strcmp(show_figures,'essential') || strcmp(show_figures,'detailed'))
            try
                figure('Name',sprintf('Slice %d - final smaps',slice_idx));
                imshow3(abs(maps),[],[4,8]); title(sprintf('Slice %d - sensitivity magnitude', slice_idx));
            catch
                if verbose, warning('Could not display final smaps for slice %d', slice_idx); end
            end
        end
        
        % store in original slice-first ordering
        smaps(slice_idx,:,:,:) = maps;  % maps is [sx, sy, Nc]
        
        % Thread-safe progress tracking
        processed_slices = processed_slices + 1;
        if mod(slice_idx, progress_interval) == 0 || slice_idx == n_slice
            elapsed_time = toc(start_time);
            percent_complete = (slice_idx / n_slice) * 100;
            estimated_total = elapsed_time * n_slice / slice_idx;
            remaining_time = estimated_total - elapsed_time;
            
            % Create progress bar
            progress_chars = floor(percent_complete / 5);
            progress_bar = ['[' repmat('=', 1, progress_chars) repmat(' ', 1, 20-progress_chars) ']'];
            
            fprintf('\rProgress: %s %.1f%% (%d/%d slices) | Elapsed: %.1fs | ETA: %.1fs', ...
                    progress_bar, percent_complete, slice_idx, n_slice, elapsed_time, remaining_time);
        end


    catch MEslice
        warning('Error processing slice %d: %s. Skipping slice.', slice_idx, MEslice.message);
        smaps(slice_idx,:,:,:) = 0;
    end
end

total_time = toc(start_time);
fprintf('\nSlice processing completed in %.2f seconds (%.3f sec/slice).\n', total_time, total_time/n_slice);

% -------------------- permute to final requested shape --------------------
% original: smaps is [n_slice, sx, sy, nc]
% desired: [nx, n_slice, ny, n_channel] == [sx, n_slice, sy, nc]

try
    smaps = permute(smaps, [2,1,3,4]);
catch
    warning('Failed to permute smaps to [nx,n_slice,ny,n_channel]. Returning original ordering.');
end
% backup save
save(sens_maps_out_file, 'smaps', '-v7.3');
niftiwrite(abs(smaps), replace(sens_maps_out_file, 'smap', 'smap_mag')); % Save only the magnitude into to NIFTI
niftiwrite(angle(smaps), replace(sens_maps_out_file, 'smap', 'smap_phase')); % Save only the magnitude into to NIFTI


%--------------------- Smooth the sensitivity map --------------------------
% Build the fslmaths command
cmd_smooth = sprintf('FSL fslmaths "%s" -s %d "%s"', sens_maps_out_file, smoothing_radius, sens_maps_out_file);

% Run the command from MATLAB
[status, cmdout] = system(cmd_smooth);

% Check status
if status == 0
    fprintf('fslmaths command executed successfully.\n');
else
    fprintf('Error: fslmaths command failed.\n', cmdout);


end
% -------------------- save NIfTI with short description --------------------
desc = sprintf('ESPIRiT smaps mag from %s on %s. ksize=[%d %d], ncalib=%d, eigThr1=%.3g, eigThr2=%.3g', ...
    twix_filename, datestr(now), ksize(1), ksize(2), ncalib, eigThresh1, eigThresh2);
if ~isempty(twix_hdr_snip)
    desc = [desc, sprintf('\nTWIX header snippet: %s', twix_hdr_snip(1:min(1000,end)))];
end

try
    hdr_R1 = niftiinfo(data_struct.R1.filepath);
    info = extractHeaderFromData( hdr_R1, smaps);
    niftiwrite(abs(smaps), replace(sens_maps_out_file, 'smap', 'smap_mag'), info); % Can try to add info later
    niftiwrite(angle(smaps), replace(sens_maps_out_file, 'smap', 'smap_phase'), info);
    if verbose, fprintf('Saved smaps to %s\n', sens_maps_out_file); end
catch MEw
    warning('niftiwrite failed (%s). Saving as MAT file instead.', ME.message);
    try
        save(sens_maps_out_file, 'smaps', '-v7.3');
        if verbose, fprintf('Saved smaps as MAT: %s\n', sens_maps_out_file); end
    catch MEsave
        error('Failed to save smaps to disk: %s', MEsave.message);
    end
end

end








function data = read_twix(twix_filename, output_path, varargin)
%READ_TWIX  Convert Siemens twix file to NIfTI format, with caching.
%
%   data = READ_TWIX(twix_filename, output_path, ...)
%
%   If the NIfTI file at output_path exists, it loads the data from it and returns.
%   Otherwise, reads the Siemens twix file (.dat) using mapVBVD, processes it,
%   and saves it as a NIfTI file at output_path, including header info.

    % Path to the folder containing mapVBVD
    mapvbvd_folder = '/data/u_kilicb_software/mapVBVD-main/'; % <-- CHANGE THIS
    addpath(mapvbvd_folder);

    if ~exist('mapVBVD', 'file')
        error('mapVBVD.m not found. Please check the mapvbvd_folder path.');
    end

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

    % Parse inputs
    p = inputParser;
    addRequired(p, 'twix_filename', @ischar);
    addRequired(p, 'output_path', @ischar);
    addParameter(p, 'subset', [], @(x) ischar(x) || iscell(x));
    addParameter(p, 'data_struct', @isstruct);
    parse(p, twix_filename, output_path, varargin{:});
    subset_params = p.Results.subset;
    data_struct = p.Results.data_struct;

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
        
    % % Save as NIfTI with header
    % fprintf('Saving processed data to NIfTI: %s\n', output_path);
    % niftiwrite(data, output_path);



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