function cfl_to_nifti(input_cfl_path, output_dir)
%CFL_TO_NIFTI Convert CFL file to NIfTI format
%   cfl_to_nifti(input_cfl_path, output_dir) converts a CFL file to NIfTI
%   format and saves it in the specified output directory with the same
%   base filename.
%
%   Inputs:
%   - input_cfl_path: Full path to the input CFL file (without extension)
%   - output_dir: Directory where the output NIfTI file will be saved
%
%   Example:
%   cfl_to_nifti('/path/to/data/image', '/output/folder/')
%   This will create '/output/folder/image_magnitude.nii' and '/output/folder/image_phase.nii'

    % Input validation
    if ~ischar(input_cfl_path) && ~isstring(input_cfl_path)
        error('input_cfl_path must be a character array or string');
    end
    if ~ischar(output_dir) && ~isstring(output_dir)
        error('output_dir must be a character array or string');
    end
    
    % Convert to char if string
    input_cfl_path = char(input_cfl_path);
    output_dir = char(output_dir);
    
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Extract filename without extension
    [~, filename, ~] = fileparts(input_cfl_path);
    
    % Check if CFL files exist
    cfl_file = [input_cfl_path '.cfl'];
    hdr_file = [input_cfl_path '.hdr'];
    
    if ~exist(cfl_file, 'file')
        error('CFL file not found: %s', cfl_file);
    end
    if ~exist(hdr_file, 'file')
        error('HDR file not found: %s', hdr_file);
    end
    
    try
        % Read CFL data using readcfl function
        % Note: This assumes you have the readcfl function available
        % If not available, you'll need to implement CFL reading or use BART toolbox
        data = readcfl(input_cfl_path);
        
        % Handle complex data - save magnitude and phase separately
        if ~isreal(data)
            % Extract magnitude and phase
            magnitude_data = abs(data);
            phase_data = angle(data);
            
            % Create output paths
            magnitude_path = fullfile(output_dir, [filename '_magnitude.nii']);
            phase_path = fullfile(output_dir, [filename '_phase.nii']);
            
            % Write magnitude NIfTI file with minimal header (let MATLAB create defaults)
            niftiwrite(single(magnitude_data), magnitude_path);
            
            % Write phase NIfTI file with minimal header (let MATLAB create defaults)
            niftiwrite(single(phase_data), phase_path);
            
            fprintf('Successfully converted CFL to NIfTI:\n');
            fprintf('  Magnitude: %s\n', magnitude_path);
            fprintf('  Phase: %s\n', phase_path);
            
        else
            % Data is real, save as single file
            output_path = fullfile(output_dir, [filename '.nii']);
            
            % Write NIfTI file with minimal header (let MATLAB create defaults)
            niftiwrite(single(data), output_path);
            
            fprintf('Successfully converted CFL to NIfTI: %s\n', output_path);
        end
        
    catch ME
        error('Failed to convert CFL to NIfTI: %s', ME.message);
    end
end


% Helper functions for CFL I/O (basic implementations)
% Note: For full functionality, consider using BART toolbox functions

function data = readcfl(filenameBase)
%READCFL Read complex data from file.
%   READCFL(filenameBase) read in reconstruction data stored in filenameBase.cfl 
%   (complex float) based on dimensions stored in filenameBase.hdr.
%
%   Parameters:
%       filenameBase:   path and filename of cfl file (without extension)
%
%   Written to edit data with the Berkeley Advanced Reconstruction Toolbox (BART).
%
% Copyright 2016. CBClab, Maastricht University.
% 2016 Tim Loderhose (t.loderhose@student.maastrichtuniversity.nl)

    dims = readReconHeader(filenameBase);

    filename = strcat(filenameBase,'.cfl');
    fid = fopen(filename);

    data_r_i = fread(fid, prod([2 dims]), '*float32');
    data_r_i = reshape(data_r_i, [2 dims]);
    data = complex(zeros(dims,'single'),0);
    data(:) = complex(data_r_i(1,:),data_r_i(2,:));

    fclose(fid);
end

function dims = readReconHeader(filenameBase)
    filename = strcat(filenameBase,'.hdr');
    fid = fopen(filename);
    
    line = getNextLine(fid);
    dims = str2num(line);
    
    fclose(fid);
end

function line = getNextLine(fid)
    line = fgetl(fid);
    while(line(1) == '#')
        line = fgetl(fid);
    end
end