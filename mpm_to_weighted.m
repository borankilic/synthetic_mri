clear all;

function B1_upscaled = resize_B1map(B1map, refArray, interpMethod)
%RESIZE_B1MAP Resize low-res B1 map to match the size of a reference array
%   B1_up = resize_B1map(B1map, refArray) uses linear interpolation by default.
%   B1_up = resize_B1map(B1map, refArray, interpMethod) specifies 'nearest',
%   'linear', 'cubic', etc.

    % Default interpolation
    if ~exist('interpMethod','var') || isempty(interpMethod)
        interpMethod = 'linear';
    end

    % Desired output size
    outSize = size(refArray);

    % Use 3D resize
    B1_upscaled = imresize3(B1map, outSize, interpMethod);
end




%---- Script to create weighted images from MT PD, R1, R2 maps ----

%%Define Hyperparameters of the system


% Specify input directories
patient_dir = '/data/dt_transfer/KP2BK/20241028_27993.cb/mpm/Results/';

% Specify output directory


% Define the range of TE, TR and flip angle values


%% Load nifti data and convert them into matrices

%LOAD_NII_MATRICES Load .nii files containing MT, PD, R1, R2 into matrices
%   data = LOAD_NII_MATRICES(patient_dir) scans the given directory for
%   NIfTI (.nii) files whose filenames contain 'MT', 'PD', 'R1', or 'R2',
%   reads them, and stores the resulting 3D matrices in the output struct
%   data with fields MT, PD, R1, R2.

% Validate input directory
if ~exist('patient_dir', 'var') || isempty(patient_dir) || ~isfolder(patient_dir)
    error('Please provide a valid patient directory.');
end

% Get list of .nii files in directory
nii_files = dir(fullfile(patient_dir, '*.nii'));
if isempty(nii_files)
    error('No .nii files found in the directory: %s', patient_dir);
end

% Initialize output struct
data = struct('MT', [], 'PD', [], 'R1', [], 'R2', [], 'B1map',[],  'B1ref', []);

% Loop over files and load matching ones
for k = 1:numel(nii_files)
    filename = nii_files(k).name;
    filepath = fullfile(patient_dir, filename);

    % Determine type based on filename
    if contains(filename, 'MT', 'IgnoreCase', true)
        data.MT = niftiread(filepath);
    elseif contains(filename, 'PD', 'IgnoreCase', true)
        data.PD = niftiread(filepath);
         ref_header = niftiinfo(filepath);
    elseif contains(filename, 'R1', 'IgnoreCase', true)
        data.R1 = niftiread(filepath);
    elseif contains(filename, 'R2', 'IgnoreCase', true)
        data.R2 = niftiread(filepath);
    elseif contains(filename, 'B1map', 'IgnoreCase', true)
        data.B1map = niftiread(filepath);    
    elseif contains(filename, 'B1ref', 'IgnoreCase', true)
        data.B1ref = niftiread(filepath);  
    else
        warning('Skipping file with unrecognized type: %s', filename);
    end
end

% Check that all fields have been loaded
fields = {'MT','PD','R1','R2'};
for i = 1:numel(fields)
    if isempty(data.(fields{i}))
        warning('No file found for %s. Data.%s is empty.', fields{i}, fields{i});
    end
end







%% Implement the Ernst Equation 
% Define TE, TR, alpha values
TE_vals = [10];    % in ms
TR_vals = [500]; % in ms
alpha_vals = [5];  % in degrees
signal_constant =2000;

% Convert degrees to radians for computation
deg2rad = @(x) x * pi / 180;

% Ernst equation: S = K * [PD] * sin(a*(B1)) * (1 - exp(-TR*R1)) / (1 - cos(a*(B1)) * exp(-TR*R1)) * exp(-TE*R2)
for TE = TE_vals
    for TR = TR_vals
        for alpha = alpha_vals
            B1_scaled = resize_B1map(data.B1map, data.PD, 'linear');
            % Adjust flip angle by B1 inhomogeneity
            eff_alpha = deg2rad(alpha) .* double(B1_scaled) ./ 100 ;

            % Compute signal using element-wise operations
            E1 = exp(-TR * data.R1 / 1000); % TR in ms -> s
            E2 = exp(-TE * data.R2 / 1000); % TE in ms -> s

            numerator = sin(eff_alpha) .* (1 - E1);
            denominator = 1 - cos(eff_alpha) .* E1;

            signal = signal_constant .* data.PD .* (numerator ./ denominator) .* E2;

            % Convert signal matrix to NIfTI and save
            output_name = sprintf('%s_TE=%d_fa=%d_TR=%d.nii','20241028_27993' , TE, alpha, TR);
            niftiwrite(single(signal), fullfile(patient_dir, output_name), ref_header);
        end
    end
end

fprintf('Finished generating synthetic images using Ernst equation.\n');



