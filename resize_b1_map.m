

function B1_upscaled = resize_b1_map(B1map_path, ref_array_path, interp_method)
%RESIZE_B1_MAP Resize B1 map to match reference array using spatial referencing
%
%   B1_upscaled = resize_b1_map(B1map_path, ref_array_path)
%   B1_upscaled = resize_b1_map(..., interp_method)
%
%   - Reads both NIfTIs, builds imref3d for source & target in real mm space
%   - Uses imwarp to resample the B1 map into the target grid
%   - Applies an affine transformation from NIfTI S-form/Q-form to align
%     the images in physical space (if available).
%   - Falls back to pixel dimensions if affine transform is not present/valid.
%   - Clamps values to [10,200] (%)
%
    if nargin < 3, interp_method = 'linear'; end
    
      % --- Input Validation and Default Handling ---
    valid_interp_methods = {'nearest', 'linear', 'cubic', 'spline'}; % Added 'spline'
    if ~any(strcmpi(interp_method, valid_interp_methods))
        if strcmpi(interp_method, 'nearestneighbour')
            interp_method = 'nearest';
        elseif  strcmpi(interp_method, 'trilinear')
            interp_method = 'linear';
        else    
        warning('Invalid interpolation method "%s", using linear.', interp_method);
        end
    end

    % --- Read NIfTI Data and Headers ---
    try
        B1vol = niftiread(B1map_path);
        B1info = niftiinfo(B1map_path); % Use niftiinfo to get full header
    catch ME
        error('Failed to read B1 map NIfTI file "%s": %s', B1map_path, ME.message);
    end

    try
        refvol = niftiread(ref_array_path);
        refinfo = niftiinfo(ref_array_path); % Use niftiinfo to get full header
    catch ME
        error('Failed to read reference NIfTI file "%s": %s', ref_array_path, ME.message);
    end
    
    % --- Handle Empty B1 Map ---
    if isempty(B1vol) || all(size(B1vol) == 0) % Check for truly empty volume
        warning('Empty B1 map: returning 100%% uniform map with reference volume dimensions.');
        B1_upscaled = ones(size(refvol), 'like', refvol) * 100;
        return;
    end
    
    % --- Build Spatial References with NIfTI Affine Transformations ---
    % Use the improved niftiCustomSpatialReference that handles missing transforms
    srcR = niftiCustomSpatialReference(B1info, B1map_path); % Pass path for better error messages
    tgtR = niftiCustomSpatialReference(refinfo, ref_array_path);
               
    % --- Run the warp/resample ---
    try
        B1_upscaled = imwarp(B1vol, srcR, 'OutputView', tgtR, ...
                             'InterpolationMethod', interp_method);
    catch ME
        % Fallback to imresize3 is less ideal for medical images as it
        % doesn't respect physical space, but keep as a last resort if imwarp fails.
        warning('imwarp failed (%s). This usually indicates an issue with spatial referencing. Falling back to imresize3. The output might be misaligned.', E.message);
        B1_upscaled = imresize3(B1vol, size(refvol), 'Method', interp_method);
    end
    
    % --- Clamp to reasonable B1 percentages ---
    % Ensure data type consistency if needed, though imwarp usually handles this
    % based on the input B1vol's type.
    B1_upscaled = cast(B1_upscaled, 'like', B1vol); % Maintain original data type
    B1_upscaled = max(B1_upscaled, 10);
    B1_upscaled = min(B1_upscaled, 200);
end

% Helper function to create imref3d from niftiinfo, with fallback
function R = niftiCustomSpatialReference(nifti_info, filepath)
    % filepath is used for more informative warning messages

    imageSize = nifti_info.ImageSize;
    
    % Attempt to use the full affine transform first
    if isfield(nifti_info, 'Transform') && ~isempty(nifti_info.Transform) && isa(nifti_info.Transform, 'affine3d')
        try
            % This is the preferred method as it handles origin, orientation, and spacing
            R = imref3d(imageSize, nifti_info.Transform);
            return; % Successfully created R using transform, exit
        catch ME
            warning('Failed to create imref3d from NIfTI Transform for "%s" (%s). Falling back to PixelDimensions.', filepath, ME.message);
            % Fall through to the next section if transform fails
        end
    else
        warning('NIfTI Transform is missing or invalid for "%s". Falling back to PixelDimensions for spatial referencing. Output may be misaligned if orientation/origin differs.', filepath);
    end
    
    % Fallback: Create imref3d using ImageSize and PixelDimensions
    % This assumes the image is axis-aligned and its origin is at (0,0,0) or
    % implicitly handled by imref3d's default behavior, which is often not
    % robust for medical images with complex coordinate systems.
    
    if isfield(nifti_info, 'PixelDimensions') && numel(nifti_info.PixelDimensions) >= 3
        pixdim_x = nifti_info.PixelDimensions(1);
        pixdim_y = nifti_info.PixelDimensions(2);
        pixdim_z = nifti_info.PixelDimensions(3);
        
        % The order for imref3d(size, deltaY, deltaX, deltaZ) or (size, deltaX, deltaY, deltaZ)
        % can be tricky. NIfTI PixelDimensions are typically (X, Y, Z).
        % For 3D images, imref3d(size, pixelExtentInWorldX, pixelExtentInWorldY, pixelExtentInWorldZ) is robust.
        R = imref3d(imageSize, pixdim_x, pixdim_y, pixdim_z);
    else
        % Last resort: no spatial info, create a default imref3d
        warning('Neither NIfTI Transform nor PixelDimensions found for "%s". Creating an identity imref3d. This will likely result in misalignment.', filepath);
        R = imref3d(imageSize); % This assumes 1mm isotropic voxels and (0,0,0) origin
    end
end