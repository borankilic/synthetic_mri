twix_filename = '/data/tu_kilicb/KP2BK/meas_MID00251_t1w_kp_mtflash3d_v1l3x_0p5_20210924-131453.dat';
output_filename = '/data/tu_kilicb/derivatives/meas_MID00251_t1w_kp_mtflash3d_v1l3x_0p5_20210924-131453.dat_kspace'; % Use chars (not strings) for filenames

twix_to_cfl(twix_filename, output_filename)

function twix_to_cfl(twix_filename, output_filename, varargin)
%TWIX_TO_CFL Convert Siemens twix file to CFL format
%   TWIX_TO_CFL(twix_filename, output_filename) reads a Siemens twix file
%   (.dat) using mapVBVD and writes it to CFL format using writecfl.
%
%   TWIX_TO_CFL(twix_filename, output_filename, 'subset', subset_params)
%   allows extracting a subset of the data before writing to CFL.
%
%   Parameters:
%       twix_filename:   path to input .dat file
%       output_filename: path and base filename for output (without extension)
%       subset_params:   optional parameter for data subsetting
%                       Can be a string specifying dimension names or
%                       a cell array with dimension specifications
%
%   Data dimensions: [Col Cha Lin Par Sli Ave Phs Eco Rep Set Seg Ida Idb Idc Idd Ide]
%   Original size:   [768 32  335 288 1   1   1   6   1   1   1   1   1   1   1   1  ]
%
%   Examples:
%       % Convert entire dataset
%       twix_to_cfl('data.dat', 'output')
%       
%       % Extract specific echo (assuming Eco is dimension 8)
%       twix_to_cfl('data.dat', 'output', 'subset', 'Eco=3')
%       
%       % Extract multiple subsets
%       twix_to_cfl('data.dat', 'output', 'subset', {'Eco=1:3', 'Cha=1:16'})

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'twix_filename', @ischar);
    addRequired(p, 'output_filename', @ischar);
    addParameter(p, 'subset', [], @(x) ischar(x) || iscell(x));
    parse(p, twix_filename, output_filename, varargin{:});
    
    subset_params = p.Results.subset;
    
    % Check if twix file exists
    if ~exist(twix_filename, 'file')
        error('Twix file not found: %s', twix_filename);
    end
    
    fprintf('Reading twix file: %s\n', twix_filename);
    
    % Read twix file using mapVBVD
    mapVBVD_path = '/data/u_kilicb_software/mapVBVD-main/';
    addpath(mapVBVD_path)
    try
        twix_data = mapVBVD(twix_filename);
    catch ME
        error('Failed to read twix file with mapVBVD: %s', ME.message);
    end
    
    % Handle case where mapVBVD returns a cell array (multiple measurements)
    if iscell(twix_data)
        if length(twix_data) > 1
            fprintf('Multiple measurements found. Using the last measurement.\n');
        end
        twix_data = twix_data{end};
    end
    
    % Extract image data
    if isfield(twix_data, 'image')
        data = twix_data.image();
        fprintf('Using image data from twix file.\n');
    elseif isfield(twix_data, 'data')
        data = twix_data.data();
        fprintf('Using raw data from twix file.\n');
    else
        error('No image or data field found in twix structure');
    end
    
    % Display original data dimensions
    original_size = size(data);
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
    
    % Apply subsetting if specified
    if ~isempty(subset_params)
        data = apply_subset(data, subset_params, dim_names);
        fprintf('Data after subsetting: %s\n', mat2str(size(data)));
    end
    
    % Squeeze singleton dimensions
    data = squeeze(data);
    final_size = size(data);
    
    fprintf('Final data dimensions after squeezing: %s\n', mat2str(final_size));
    fprintf('Data type: %s\n', class(data));
    fprintf('Data range: [%.6f + %.6fj, %.6f + %.6fj]\n', ...
            real(min(data(:))), imag(min(data(:))), ...
            real(max(data(:))), imag(max(data(:))));
    
    % Write to CFL format
    fprintf('Writing to CFL files: %s.cfl and %s.hdr\n', output_filename, output_filename);
    try
        writecfl(output_filename, data);
        fprintf('Successfully wrote CFL files.\n');
    catch ME
        error('Failed to write CFL files: %s', ME.message);
    end
end

function data_subset = apply_subset(data, subset_params, dim_names)
%APPLY_SUBSET Apply subsetting parameters to data
    
    % Convert single string to cell array
    if ischar(subset_params)
        subset_params = {subset_params};
    end
    
    % Create indexing cell array for all dimensions
    indices = cell(1, ndims(data));
    for i = 1:ndims(data)
        indices{i} = ':';
    end
    
    % Parse each subset parameter
    for i = 1:length(subset_params)
        param = subset_params{i};
        
        % Parse dimension=range format (e.g., 'Eco=1:3' or 'Cha=5')
        parts = strsplit(param, '=');
        if length(parts) ~= 2
            error('Invalid subset parameter format. Use "DimName=range" (e.g., "Eco=1:3")');
        end
        
        dim_name = strtrim(parts{1});
        range_str = strtrim(parts{2});
        
        % Find dimension index
        dim_idx = find(strcmp(dim_names, dim_name), 1);
        if isempty(dim_idx)
            error('Unknown dimension name: %s. Valid names: %s', ...
                  dim_name, strjoin(dim_names, ', '));
        end
        
        if dim_idx > ndims(data)
            error('Dimension %s (index %d) exceeds data dimensions (%d)', ...
                  dim_name, dim_idx, ndims(data));
        end
        
        % Parse range
        try
            range_indices = eval(range_str);
        catch
            error('Invalid range specification: %s', range_str);
        end
        
        % Validate range
        max_idx = size(data, dim_idx);
        if any(range_indices < 1) || any(range_indices > max_idx)
            error('Range [%s] out of bounds for dimension %s (size: %d)', ...
                  range_str, dim_name, max_idx);
        end
        
        indices{dim_idx} = range_indices;
        fprintf('Subsetting %s (dim %d): %s\n', dim_name, dim_idx, range_str);
    end
    
    % Apply subsetting
    data_subset = data(indices{:});
end

function writecfl(filenameBase,data)
%WRITECFL  Write complex data to file.
%   WRITECFL(filenameBase, data) writes reconstruction data to 
%   filenameBase.cfl (complex float) and its dimensions to filenameBase.hdr.
%
%   Written to edit data with the Berkeley Advanced Reconstruction Toolbox (BART).
%
%   Parameters:
%       filenameBase:   path and filename of cfl file (without extension)
%       data:           array/matrix to be written
%
% Copyright 2013. Joseph Y Cheng.
% Copyright 2016. CBClab, Maastricht University.
% 2012 Joseph Y Cheng (jycheng@mrsrl.stanford.edu).
% 2016 Tim Loderhose (t.loderhose@student.maastrichtuniversity.nl).
    dims = size(data);
    writeReconHeader(filenameBase,dims);
    filename = strcat(filenameBase,'.cfl');
    fid = fopen(filename,'w');
    
    data = data(:);
    
    fwrite(fid,[real(data)'; imag(data)'],'float32');
    fclose(fid);
end

function writeReconHeader(filenameBase,dims)
    filename = strcat(filenameBase,'.hdr');
    fid = fopen(filename,'w');
    fprintf(fid,'# Dimensions\n');
    for N=1:length(dims)
        fprintf(fid,'%d ',dims(N));
    end
    if length(dims) < 5
        for N=1:(5-length(dims))
            fprintf(fid,'1 ');
        end
    end
    fprintf(fid,'\n');
    
    fclose(fid);
end