function view_brain_slices(volume)
% view_brain_slices - Visualize orthogonal slices from a 3D brain volume
%
% Syntax: view_brain_slices(volume)
%
% Input:
%   volume - 3D matrix (e.g., from a NIfTI brain scan), size: [X Y Z]
%
% Example usage:
%   nifti = niftiread('brain.nii');
%   view_brain_slices(nifti);

    if ndims(volume) ~= 3
        error('Input must be a 3D volume.');
    end

    % Get dimensions
    [x, y, z] = size(volume);

    % Find middle slices
    mid_x = round(x / 2);
    mid_y = round(y / 2);
    mid_z = round(z / 2);

    % Create a figure with 3 subplots
    figure;
    
    % Axial view (top-down)
    subplot(1, 3, 1);
    imagesc(squeeze(volume(:, :, mid_z)));
    axis image off;
    colormap gray;
    title('Axial (Z)');

    % Coronal view (front-back)
    subplot(1, 3, 2);
    imagesc(squeeze(volume(:, mid_y, :))');
    axis image off;
    colormap gray;
    title('Coronal (Y)');

    % Sagittal view (side-side)
    subplot(1, 3, 3);
    imagesc(squeeze(volume(mid_x, :, :))');
    axis image off;
    colormap gray;
    title('Sagittal (X)');

end
