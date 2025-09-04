function kspace_to_nifti(input_dir, varargin)
%KSPACE_TO_NIFTI Convert k-space data back to image space using forward FFT
%   Reconstructs images from both magnitude and phase components of k-space data
%
%   KSPACE_TO_NIFTI(input_dir) processes k-space files in the input
%   directory, converts them back to image space using forward FFT, and saves
%   them with modified naming convention.
%
%   KSPACE_TO_NIFTI(input_dir, 'param', value) allows additional parameters:
%
%   Parameters:
%       input_dir:       Input directory containing k-space NIfTI files
%       'output_dir':    Custom output directory (optional)
%       'fft_shift':     Apply fftshift (default: true)
%       'normalize':     Normalize reconstructed image (default: false)
%       'verbose':       Display progress messages (default: true)
%       'output_suffix': Suffix for output files (default: 'fft_image')
%
%   Expected input files:
%       - *_synthkspace_magnitude_*.nii (magnitude component)
%       - *_synthkspace_phase_*.nii (phase component)
%   OR
%       - *_synthkspace_*.nii (magnitude only or real data)
%
%   Output files:
%       - *_synthkspace_*_fft_image.nii (reconstructed images)

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'input_dir', @ischar);
    addParameter(p, 'output_dir', '', @ischar);
    addParameter(p, 'fft_shift', true, @islogical);
    addParameter(p, 'normalize', false, @islogical);
    addParameter(p, 'verbose', true, @islogical);
    addParameter(p, 'output_suffix', 'recon_image', @ischar);
    parse(p, input_dir, varargin{:});
    
    fft_shift = p.Results.fft_shift;
    normalize = p.Results.normalize;
    verbose = p.Results.verbose;
    output_suffix = p.Results.output_suffix;
    
    % Ensure input directory exists
    if ~exist(input_dir, 'dir')
        error('Input directory does not exist: %s', input_dir);
    end
    
    % Create output directory path
    if isempty(p.Results.output_dir)
        output_dir = strrep(input_dir, 'kspace', 'reconstructed');
    else
        output_dir = p.Results.output_dir;
    end
    
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
        if verbose
            fprintf('Created output directory: %s\n', output_dir);
        end
    end
    
    % Find k-space files
    [mag_files, phase_files, single_files] = find_kspace_files(input_dir);
    
    total_pairs = length(mag_files);
    total_single = length(single_files);
    
    if total_pairs == 0 && total_single == 0
        warning('No k-space files found in directory: %s', input_dir);
        return;
    end
    
    if verbose
        fprintf('Found %d magnitude/phase pairs\n', total_pairs);
        fprintf('Found %d single k-space files\n', total_single);
        fprintf('Input directory:  %s\n', input_dir);
        fprintf('Output directory: %s\n', output_dir);
        fprintf('FFT shift: %s\n', mat2str(fft_shift));
        fprintf('Normalize: %s\n', mat2str(normalize));
        fprintf('Processing files...\n\n');
    end
    
    success_count = 0;
    error_count = 0;
    
    % Process magnitude/phase pairs
    for i = 1:total_pairs
        try
            if verbose
                fprintf('[%d/%d] Processing pair: %s + %s\n', i, total_pairs, ...
                    mag_files(i).name, phase_files(i).name);
            end
            
            % Load magnitude and phase data
            mag_filepath = fullfile(input_dir, mag_files(i).name);
            phase_filepath = fullfile(input_dir, phase_files(i).name);
            
            mag_nii = load_nii(mag_filepath);
            phase_nii = load_nii(phase_filepath);
            
            magnitude = double(mag_nii.img);
            phase = double(phase_nii.img);
            
            if verbose
                fprintf('         Magnitude range: [%.2e, %.2e]\n', min(magnitude(:)), max(magnitude(:)));
                fprintf('         Phase range: [%.2f, %.2f] radians\n', min(phase(:)), max(phase(:)));
            end
            
            % Reconstruct complex k-space data
            kspace_complex = magnitude .* exp(1i * phase);
            
            % Convert back to image space
            if verbose
                fprintf('         Computing forward FFT...\n');
            end
            tic;
            [reconstructed_image, fft_details] = convert_kspace_to_image(kspace_complex, fft_shift, normalize);
            fft_time = toc;
            
            if verbose
                fprintf('         FFT completed in %.2f seconds\n', fft_time);
                fprintf('         Reconstructed image range: [%.2f, %.2f]\n', ...
                    min(reconstructed_image(:)), max(reconstructed_image(:)));
            end
            
            % Generate output filename
            base_name = extract_base_name(mag_files(i).name);
            output_filename = sprintf('%s_%s.nii', base_name, output_suffix);
            output_filepath = fullfile(output_dir, output_filename);
            
            % Save reconstructed image
            save_reconstructed_image(mag_nii, reconstructed_image, output_filepath, verbose);
            
            % Process JSON metadata
            process_reconstruction_json(mag_filepath, output_dir, output_filename, fft_details, verbose);
            
            success_count = success_count + 1;
            if verbose
                fprintf('         ✓ Saved reconstructed image: %s\n\n', output_filename);
            end
            
        catch ME
            error_count = error_count + 1;
            if verbose
                fprintf('         ✗ Error processing pair: %s\n', ME.message);
                fprintf('\n');
            end
        end
    end
    
    % Process single k-space files (magnitude only or real data)
    for i = 1:total_single
        try
            if verbose
                fprintf('[%d/%d] Processing single file: %s\n', i, total_single, single_files(i).name);
            end
            
            % Load k-space data
            kspace_filepath = fullfile(input_dir, single_files(i).name);
            kspace_nii = load_nii(kspace_filepath);
            kspace_data = double(kspace_nii.img);
            
            if verbose
                fprintf('         K-space range: [%.2e, %.2e]\n', min(kspace_data(:)), max(kspace_data(:)));
            end
            
            % Convert back to image space (treating as magnitude only for complex reconstruction)
            if verbose
                fprintf('         Computing forward FFT (magnitude only)...\n');
            end
            tic;
            [reconstructed_image, fft_details] = convert_kspace_to_image(kspace_data, fft_shift, normalize);
            fft_time = toc;
            
            if verbose
                fprintf('         FFT completed in %.2f seconds\n', fft_time);
                fprintf('         Reconstructed image range: [%.2f, %.2f]\n', ...
                    min(reconstructed_image(:)), max(reconstructed_image(:)));
            end
            
            % Generate output filename
            base_name = extract_base_name(single_files(i).name);
            output_filename = sprintf('%s_%s.nii', base_name, output_suffix);
            output_filepath = fullfile(output_dir, output_filename);
            
            % Save reconstructed image
            save_reconstructed_image(kspace_nii, reconstructed_image, output_filepath, verbose);
            
            % Process JSON metadata
            process_reconstruction_json(kspace_filepath, output_dir, output_filename, fft_details, verbose);
            
            success_count = success_count + 1;
            if verbose
                fprintf('         ✓ Saved reconstructed image: %s\n\n', output_filename);
            end
            
        catch ME
            error_count = error_count + 1;
            if verbose
                fprintf('         ✗ Error processing file: %s\n', ME.message);
                fprintf('\n');
            end
        end
    end
    
    % Summary
    if verbose
        fprintf('Reconstruction complete!\n');
        fprintf('Successfully processed: %d files\n', success_count);
        fprintf('Errors encountered: %d files\n', error_count);
        fprintf('Output directory: %s\n', output_dir);
    end
end

function [mag_files, phase_files, single_files] = find_kspace_files(input_dir)
%FIND_KSPACE_FILES Find and pair magnitude/phase files, identify single files
    
    % Find all k-space files
    patterns = {"*synthkspace*.nii", "*kspace*.nii", "*synthkspc*.nii", "*synthetic_kspace*.nii"};
    all_files = [];
    for p = 1:numel(patterns)
        files = dir(fullfile(input_dir, patterns{p}));
        all_files = [all_files; files];
    end

    
    % Separate magnitude, phase, and single files
    mag_files = [];
    phase_files = [];
    single_files = [];
    
    for i = 1:length(all_files)
        filename = all_files(i).name;
        if contains(filename, 'mag')
            mag_files = [mag_files; all_files(i)];
        elseif contains(filename, 'phase')
            phase_files = [phase_files; all_files(i)];
        else
            single_files = [single_files; all_files(i)];
        end
    end
    
    % Sort files to ensure proper pairing
    if ~isempty(mag_files)
        [~, mag_idx] = sort({mag_files.name});
        mag_files = mag_files(mag_idx);
    end
    
    if ~isempty(phase_files)
        [~, phase_idx] = sort({phase_files.name});
        phase_files = phase_files(phase_idx);
    end
    
    % Verify pairing
    if length(mag_files) ~= length(phase_files)
        warning('Mismatch between magnitude (%d) and phase (%d) files', ...
            length(mag_files), length(phase_files));
        % Keep only matching pairs
        min_len = min(length(mag_files), length(phase_files));
        mag_files = mag_files(1:min_len);
        phase_files = phase_files(1:min_len);
    end
end

function base_name = extract_base_name(filename)
%EXTRACT_BASE_NAME Extract base name for output file generation
    % Remove file extension
    [~, name, ~] = fileparts(filename);
    
    % Remove magnitude/phase suffixes
    base_name = strrep(name, '_magnitude', '');
    base_name = strrep(base_name, '_phase', '');
    
    % Ensure we have synthkspace in the name
    if ~contains(base_name, 'synthkspace')
        base_name = strrep(name, 'synthim', 'synthkspace');
    end
end

function [reconstructed_image, fft_details] = convert_kspace_to_image(kspace_data, apply_fft_shift, normalize)
%CONVERT_KSPACE_TO_IMAGE Convert k-space data back to image space using forward FFT
    
    % Store FFT processing details
    fft_details = struct();
    fft_details.transformation = 'forward_fft';
    fft_details.domain_from = 'k_space';
    fft_details.domain_to = 'image_space';
    fft_details.fft_shift_applied = apply_fft_shift;
    fft_details.normalization_applied = normalize;
    fft_details.processing_timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fft_details.input_was_complex = ~isreal(kspace_data);
    
    % Handle different dimensionalities
    ndims_data = ndims(kspace_data);
    data_size = size(kspace_data);
    fft_details.kspace_dimensions = data_size;
    fft_details.number_of_dimensions = ndims_data;
    
    if ndims_data == 2
        % 2D k-space
        fft_details.fft_type = '2D_forward_fft';
        if apply_fft_shift
            reconstructed_image = fftshift(ifft2(ifftshift(kspace_data)));
        else
            reconstructed_image = ifft2(kspace_data);
        end
        
    elseif ndims_data == 3
        % 3D k-space
        fft_details.fft_type = '3D_forward_fft';
        if apply_fft_shift
            reconstructed_image = fftshift(ifftn(ifftshift(kspace_data)));
        else
            reconstructed_image = ifftn(kspace_data);
        end
        
    elseif ndims_data == 4
        % 4D k-space data - apply 3D forward FFT to spatial dimensions only
        fft_details.fft_type = '3D_forward_fft_per_volume';
        fft_details.temporal_dimension = size(kspace_data, 4);
        reconstructed_image = zeros(size(kspace_data));
        for t = 1:size(kspace_data, 4)
            vol = kspace_data(:,:,:,t);
            if apply_fft_shift
                reconstructed_image(:,:,:,t) = fftshift(ifftn(ifftshift(vol)));
            else
                reconstructed_image(:,:,:,t) = ifftn(vol);
            end
        end
        
    else
        error('Unsupported k-space dimensionality: %d', ndims_data);
    end
    
    % Take real part if result is complex (common in reconstruction)
    if ~isreal(reconstructed_image)
        fft_details.output_was_complex = true;
        fft_details.complex_magnitude_range = [min(abs(reconstructed_image(:))), max(abs(reconstructed_image(:)))];
        fft_details.complex_phase_range = [min(angle(reconstructed_image(:))), max(angle(reconstructed_image(:)))];
        
        % For medical images, typically take the magnitude
        reconstructed_image = abs(reconstructed_image);
        fft_details.final_output = 'magnitude_only';
    else
        fft_details.output_was_complex = false;
        fft_details.final_output = 'real_data';
    end
    
    % Apply normalization if requested
    if normalize
        reconstructed_image = reconstructed_image / max(reconstructed_image(:));
        fft_details.normalization_factor = 1 / max(reconstructed_image(:));
    end
    
    % Add reconstruction statistics
    fft_details.reconstructed_image_range = [min(reconstructed_image(:)), max(reconstructed_image(:))];
    fft_details.reconstruction_quality = struct();
    fft_details.reconstruction_quality.mean_intensity = mean(reconstructed_image(:));
    fft_details.reconstruction_quality.std_intensity = std(reconstructed_image(:));
end

function save_reconstructed_image(original_nii, image_data, output_filepath, verbose)
%SAVE_RECONSTRUCTED_IMAGE Save reconstructed image as NIfTI
    
    % Prepare output NIfTI structure
    output_nii = original_nii;
    output_nii.img = single(image_data); % Use single precision
    
    % Update header information
    output_nii.hdr.dime.datatype = 16; % Single precision float
    output_nii.hdr.dime.bitpix = 32;
    
    % Compute statistics
    output_nii.hdr.dime.cal_max = max(image_data(:));
    output_nii.hdr.dime.cal_min = min(image_data(:));
    
    % Save NIfTI file
    if verbose
        fprintf('         Writing reconstructed image...\n');
    end
    tic;
    save_nii(image_data, output_filepath);
    save_time = toc;
    
    if verbose
        fprintf('         Image saved in %.2f seconds\n', save_time);
    end
end

function process_reconstruction_json(input_filepath, output_dir, output_filename, fft_details, verbose)
%PROCESS_RECONSTRUCTION_JSON Create JSON metadata for reconstructed images
    
    [input_dir, base_name, ~] = fileparts(input_filepath);
    json_input_filepath = fullfile(input_dir, [base_name '.json']);
    
    % Try to read original JSON if it exists
    original_data = struct();
    if exist(json_input_filepath, 'file')
        try
            json_text = fileread(json_input_filepath);
            original_data = jsondecode(json_text);
        catch
            if verbose
                fprintf('         Warning: Could not read original JSON file\n');
            end
        end
    end
    
    % Create enhanced JSON for reconstructed image
    enhanced_data = original_data;
    enhanced_data.DataType = 'reconstructed_image';
    enhanced_data.ReconstructionMethod = 'FFT_from_kspace';
    enhanced_data.FourierTransform = fft_details;
    
    % Add reconstruction history
    if ~isfield(enhanced_data, 'ProcessingHistory')
        enhanced_data.ProcessingHistory = [];
    end
    
    reconstruction_step = struct();
    reconstruction_step.Step = 'KSpace_to_Image_Reconstruction';
    reconstruction_step.ProcessingDate = fft_details.processing_timestamp;
    reconstruction_step.InputDomain = fft_details.domain_from;
    reconstruction_step.OutputDomain = fft_details.domain_to;
    reconstruction_step.ReconstructionQuality = fft_details.reconstruction_quality;
    
    if isempty(enhanced_data.ProcessingHistory)
        enhanced_data.ProcessingHistory = reconstruction_step;
    else
        enhanced_data.ProcessingHistory(end+1) = reconstruction_step;
    end
    
    % Write JSON file
    [~, json_base, ~] = fileparts(output_filename);
    json_output_filepath = fullfile(output_dir, [json_base '.json']);
    
    try
        json_text = jsonencode(enhanced_data, 'PrettyPrint', true);
    catch
        json_text = jsonencode(enhanced_data);
    end
    
    fid = fopen(json_output_filepath, 'w');
    if fid ~= -1
        fprintf(fid, '%s', json_text);
        fclose(fid);
        if verbose
            fprintf('         ✓ JSON metadata created\n');
        end
    end
end

% Include the same helper functions for NIfTI I/O
function nii = load_nii(filename)
%LOAD_NII Load NIfTI file using available toolbox
    try
        nii = load_nii(filename);
    catch
        try
            V = spm_vol(filename);
            img = spm_read_vols(V);
            nii.hdr = V;
            nii.img = img;
            nii.hdr.dime.datatype = V.dt(1);
            nii.hdr.dime.bitpix = spm_type(V.dt(1), 'bits');
            nii.hdr.dime.cal_max = max(img(:));
            nii.hdr.dime.cal_min = min(img(:));
        catch
            try
                img = niftiread(filename);
                info = niftiinfo(filename);
                nii.img = img;
                nii.hdr = info;
            catch
                error('Unable to read NIfTI file. Please ensure you have NIfTI toolbox, SPM, or MATLAB R2017b+ installed.');
            end
        end
    end
end

function save_nii(nii, filename)
%SAVE_NII Save NIfTI file using available toolbox
    try
        save_nii(nii, filename);
    catch
        try
            V = nii.hdr;
            V.fname = filename;
            V.dt = [nii.hdr.dime.datatype 0];
            spm_write_vol(V, nii.img);
        catch
            try
                niftiwrite(nii.img, filename, nii.hdr);
            catch
                error('Unable to save NIfTI file. Please ensure you have NIfTI toolbox, SPM, or MATLAB R2017b+ installed.');
            end
        end
    end
end