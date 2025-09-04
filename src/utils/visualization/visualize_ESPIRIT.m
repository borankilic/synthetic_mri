    % Input directory - modify this path as needed
input_dir = '/data/tu_kilicb/derivatives/slice_pngs'; % Change this to your actual directory

% Initialize variables
num_slices = 32;
images = cell(num_slices, 1);
first_image_loaded = false;

% Load all images from slice 00 to 31
for i = 0:31
    % Format filename with zero-padding
    filename = sprintf('MID00130_FID06609_t1w_kp_mtflash3d_v1s_0p65_smaps_kund_soft_slice190-%04d.png', i);
    filepath = fullfile(input_dir, filename);

    whole_image_filepath = fullfile(input_dir, 'calib_mag.nii');
    
    % Check if file exists and load
    if exist(filepath, 'file')
        img = imread(filepath);
        
        % Convert to double for calculations
        if size(img, 3) == 3
            img = rgb2gray(img); % Convert RGB to grayscale if needed
        end
        img = double(img);
        
        images{i+1} = img;
        
        % Initialize RMS array with first image dimensions
        if ~first_image_loaded
            [rows, cols] = size(img);
            rms_sum = zeros(rows, cols);
            first_image_loaded = true;
        end
        
        fprintf('Loaded: %s\n', filename);
    else
        fprintf('Warning: File not found: %s\n', filename);
        images{i+1} = [];
    end
end

% Calculate Root Mean Square (RMS) of all images
valid_images = ~cellfun(@isempty, images);
num_valid = sum(valid_images);

if num_valid > 0
    % Sum of squares
    sum_of_squares = zeros(rows, cols);
    for i = 1:num_slices
        if ~isempty(images{i})
            sum_of_squares = sum_of_squares + images{i}.^2;
        end
    end
    
    % Calculate RMS
    rms_image = sqrt(sum_of_squares / num_valid);
    
    fprintf('RMS calculated from %d valid images\n', num_valid);
else
    error('No valid images found!');
end

% Display all slices in 4x8 grid
figure('Name', 'All Channels (4x8 Grid)', 'Position', [100, 100, 1200, 600]);
for i = 1:num_slices
    subplot(4, 8, i);
    
    if ~isempty(images{i})
        imagesc(images{i});
        colormap gray;
        axis image off;
        title(sprintf('Channel %02d', i-1), 'FontSize', 8);
    else
        % Show empty plot for missing images
        plot(0, 0);
        axis off;
        title(sprintf('Channel %02d\n(Missing)', i-1), 'FontSize', 8, 'Color', 'red');
    end
end
sgtitle('All Calibration Slices', 'FontSize', 14, 'FontWeight', 'bold');

% Display RMS image in separate figure
figure('Name', 'RMS Image', 'Position', [200, 200, 600, 500]);
imagesc(rms_image);
colormap gray;
colorbar;
axis image;
title(sprintf('Root Mean Square Image\n(From %d channel)', num_valid), ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('Pixel Column');
ylabel('Pixel Row');

% Display statistics
fprintf('\n--- Image Statistics ---\n');
fprintf('Number of valid images: %d/%d\n', num_valid, num_slices);
fprintf('Image dimensions: %d x %d\n', rows, cols);
fprintf('RMS image - Min: %.2f, Max: %.2f, Mean: %.2f\n', ...
        min(rms_image(:)), max(rms_image(:)), mean(rms_image(:)));