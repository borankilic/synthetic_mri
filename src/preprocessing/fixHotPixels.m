function fixedImage = fixHotPixels(inputImage, threshold_mad, threshold_mehod)
%FIXHOTPIXELS3D_MAD_Capped Identifies and caps "hot pixels" in a 3D image using MAD.
%   fixedImage = FIXHOTPIXELS3D_MAD_Capped(inputImage, threshold_mad)
%
%   This function takes a 3D input image and identifies voxels whose values
%   are significantly higher than the image median, based on the Median
%   Absolute Deviation (MAD). These "hot pixels" are then directly set to
%   the calculated hot pixel threshold.
%
%   Inputs:
%       inputImage: A 3D numerical array representing the image volume.
%                   Can be single or double precision.
%       threshold_mad: A scalar value indicating how many MADs away from
%                      the median a pixel must be considered a hot pixel.
%                      A common robust threshold is 3, but for very extreme
%                      outliers, 5 or more might be appropriate.
%
%   Output:
%       fixedImage: The 3D image with hot pixels capped at the threshold.
%
%   Example:
%       % Create a sample 3D image with a hot pixel group
%       img3D = rand(20, 20, 20) * 0.5; % Base values around 0.5
%       img3D(9:11, 9:11, 9:11) = 1e10; % A group of hot pixels
%       
%       % Fix the hot pixel using a MAD threshold (e.g., 5 MADs)
%       fixed_img3D = fixHotPixels3D_MAD_Capped(img3D, 5); 
%       
%       % You can visualize a slice (optional)
%       figure;
%       subplot(1,2,1); imagesc(img3D(:,:,10)); title('Original Slice 10'); colorbar;
%       subplot(1,2,2); imagesc(fixed_img3D(:,:,10)); title('Fixed Slice 10'); colorbar;

if nargin < 2
    threshold_mad = 5; % Default MAD threshold if not provided
end

if ndims(inputImage) ~= 3 || ~isnumeric(inputImage)
    error('Input image must be a 3D numerical array.');
end

fixedImage = double(inputImage); % Work with double precision

% --- Calculate Median and MAD (Robust Statistics) ---
% Flatten the image for median and MAD calculation
imageDataFlat = fixedImage(:);

% Remove NaN values for robust calculation
imageDataFlat(isnan(imageDataFlat)) = []; 

if strcmp(threshold_mehod, 'median')
    img_median = median(imageDataFlat);
    
    % Calculate Median Absolute Deviation (MAD)
    abs_deviations = abs(imageDataFlat - img_median);
    img_mad = median(abs_deviations);
    
    % Define the hot pixel threshold
    hot_pixel_threshold = img_median + threshold_mad * img_mad;
    fprintf('Image Median: %.4g\n', img_median);
else
    img_mean = mean(imageDataFlat);
    
    % Calculate Median Absolute Deviation (MAD)
    abs_deviations = abs(imageDataFlat - img_mean);
    img_mad = median(abs_deviations);
    
    % Define the hot pixel threshold
    hot_pixel_threshold = img_mean + threshold_mad * img_mad;
    fprintf('Image Median: %.4g\n', img_mean);
end



fprintf('Image MAD: %.4g\n', img_mad);
fprintf('Hot Pixel Detection Threshold (Median + %.1f * MAD): %.4g\n', threshold_mad, hot_pixel_threshold);

% Find hot pixels: values greater than the calculated threshold
hotPixelMask = fixedImage > hot_pixel_threshold;

% --- Core change: Directly set hot pixels to the threshold ---
numHotPixels = sum(hotPixelMask(:)); % Count identified hot pixels
fprintf('Found %d potential hot pixels in the 3D volume.\n', numHotPixels);

fixedImage(hotPixelMask) = 0;

fprintf('All identified hot pixels have been capped at the threshold value (%.4g).\n', hot_pixel_threshold);

end