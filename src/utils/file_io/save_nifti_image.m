function save_nifti_image(image_data, output_path, ref_header)
%SAVE_NIFTI_IMAGE Save image data as NIfTI file
%   save_nifti_image(image_data, output_path, ref_header) saves the image
%   data to a NIfTI file using the reference header for spatial information.

    % Input validation
    if ~isnumeric(image_data)
        error('image_data must be numeric');
    end
    if ~ischar(output_path) && ~isstring(output_path)
        error('output_path must be a character array or string');
    end
    if ~isstruct(ref_header)
        error('ref_header must be a struct');
    end
    
    % Ensure data is in single precision
    image_data = single(image_data);
    
    % Create output directory if it doesn't exist
    output_dir = fileparts(output_path);
    if ~isempty(output_dir) && ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Update header information for the new image
    output_header = ref_header;
    output_header.Datatype = 'single';
    output_header.BitsPerPixel = 32;
    output_header.ImageSize = size(image_data);
    
    % Add description
    [~, filename, ~] = fileparts(output_path);
    params = split(filename, '_');
    params = params(end-2:end);
    params = strjoin(params, ', ');
    output_header.Description = sprintf('Synthetic image: %s', params);
    
    try
        % Save the NIfTI file
        niftiwrite_numbered(image_data, output_path, output_header);
        
    catch ME
        error('Failed to save NIfTI file "%s": %s', output_path, ME.message);
    end
end

function niftiwrite_numbered(img, filename, info, varargin)
    % niftiwrite_numbered: Save a NIfTI file with automatic numbering
    % If filename.nii exists, it will save as filename(1).nii, filename(2).nii, etc.
    
    [pathstr, name, ext] = fileparts(filename);
    if isempty(ext)
        ext = '.nii';
    end
    
    fullFile = fullfile(pathstr, [name ext]);
    count = 1;
    
    % Increment number until a free filename is found
    while exist(fullFile, 'file')
        fullFile = fullfile(pathstr, sprintf('%s(%d)%s', name, count, ext));
        count = count + 1;
    end
    
    % Save using niftiwrite
    niftiwrite(img, fullFile, info, varargin{:});
    
    fprintf('Saved NIfTI file as: %s\n', fullFile);
end
