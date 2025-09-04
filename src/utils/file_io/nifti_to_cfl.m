function nifti_to_cfl(input_nifti_path, output_dir)
%NIFTI_TO_CFL Convert NIfTI file to CFL format
%   nifti_to_cfl(input_nifti_path, output_dir) converts a NIfTI file to CFL
%   format and saves it in the specified output directory with the same
%   base filename.
%
%   Inputs:
%   - input_nifti_path: Full path to the input NIfTI file (.nii or .nii.gz)
%   - output_dir: Directory where the output CFL files will be saved
%
%   Example:
%   nifti_to_cfl('/path/to/data/image.nii', '/output/folder/')
%   This will create '/output/folder/image.cfl' and '/output/folder/image.hdr'

    % Input validation
    if ~ischar(input_nifti_path) && ~isstring(input_nifti_path)
        error('input_nifti_path must be a character array or string');
    end
    if ~ischar(output_dir) && ~isstring(output_dir)
        error('output_dir must be a character array or string');
    end
    
    % Convert to char if string
    input_nifti_path = char(input_nifti_path);
    output_dir = char(output_dir);
    
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Check if NIfTI file exists
    if ~exist(input_nifti_path, 'file')
        error('NIfTI file not found: %s', input_nifti_path);
    end
    
    % Extract filename without extension (handle .nii.gz)
    [~, filename, ext] = fileparts(input_nifti_path);
    if strcmpi(ext, '.gz')
        [~, filename, ~] = fileparts(filename);
    end
    
    try
        % Read NIfTI data and header
        data = niftiread(input_nifti_path);
        header_info = niftiinfo(input_nifti_path);
        
        % Create output CFL filename (without extension)
        output_path = fullfile(output_dir, filename);
        
        % Write CFL data using writecfl function
        % Note: This assumes you have the writecfl function available
        % If not available, you'll need to implement CFL writing or use BART toolbox
        writecfl(output_path, data);
        
        fprintf('Successfully converted NIfTI to CFL: %s.cfl/.hdr\n', output_path);
        
    catch ME
        error('Failed to convert NIfTI to CFL: %s', ME.message);
    end
end
function writecfl(filename, data)
%WRITECFL Write CFL file format (simplified version)
%   This is a basic implementation. For full BART compatibility,
%   use the official BART toolbox writecfl function.

    % Get dimensions
    dims = size(data);
    
    % Write header file
    hdr_file = [filename '.hdr'];
    fid = fopen(hdr_file, 'w');
    if fid == -1
        error('Cannot create header file: %s', hdr_file);
    end
    
    fprintf(fid, '# Dimensions\n');
    fprintf(fid, '%d ', dims);
    fprintf(fid, '\n');
    fclose(fid);
    
    % Write data file
    cfl_file = [filename '.cfl'];
    fid = fopen(cfl_file, 'w');
    if fid == -1
        error('Cannot create CFL file: %s', cfl_file);
    end
    
    % Convert complex data to interleaved real/imaginary
    data_flat = data(:);
    raw_data = zeros(2 * length(data_flat), 1);
    raw_data(1:2:end) = real(data_flat);
    raw_data(2:2:end) = imag(data_flat);
    
    % Write as float32
    fwrite(fid, raw_data, 'float32');
    fclose(fid);
end