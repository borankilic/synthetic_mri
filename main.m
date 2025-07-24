clear; clc;

%% Set up paths
patient_dir = '/data/dt_transfer/KP2BK/20241028_27993.cb/mpm/';
output_dir = '/data/tu_kilicb/synthetic_imgs_generated_fsl_mask/';

% %% Set up TR, TE, FA values
% TE_values = [5, 10, 15, 20];           % Echo times (ms)
% TR_values = [300, 500, 800, 1200];     % Repetition times (ms)
% alpha_values = [5, 10, 15, 20, 30];    % Flip angles (degrees)

%% Set up TR, TE, FA values
TE_values = [5, 10, 15, 20];           % Echo times (ms)
TR_values = [300, 500, 800, 1200];     % Repetition times (ms)
alpha_values = [5, 10, 15, 20, 30];    % Flip angles (degrees)

%% Basic usage with default parameters
try
    % Generate images with custom parameters
    generate_synthetic_images(patient_dir, output_dir, ...
        'TE_vals', TE_values, ...
        'TR_vals', TR_values, ...
        'alpha_vals', alpha_values, ...
        'signal_constant', 1500, ...  
        'interp_method', 'nearestneighbour', ... ;  % {trilinear,nearestneighbour,sinc,spline}
        'extract_brain', 'fsl')  % Brain extraction using 'fsl' (c1+c2+c3>0.95 and erode&dilate) , 'bet' or 'no' for no extraction
catch ME
    fprintf('Error: %s\n', ME.message);
    return;
end