%% B1 Map Registration to Quantitative Maps using hMRI DARTEL
% This script registers B1 maps to match the resolution of R1, R2, PD maps
% using diffeomorphic deformation via hMRI toolbox's DARTEL implementation
%
% Requirements:
% - hMRI toolbox installed and in MATLAB path
% - SPM12 with DARTEL toolbox
% - B1 map and anatomical B1 reference
% - R1, R2, PD quantitative maps (target resolution)
% - Tissue probability maps (c1, c2, c3 files)

%% Initialize
clear; clc;


%% Define file paths - MODIFY THESE PATHS FOR YOUR DATA
% Subject directories
subj_dir = "/data/dt_transfer/KP2BK/20241028_27993.cb/mpm/";  % Main subject directory
b1_dir = fullfile(subj_dir, 'B1mapCalc/');                      % B1 map directory
quant_dir = fullfile(subj_dir, 'Results/');         % Quantitative maps directory
seg_dir = fullfile(subj_dir ,'MPMCalc/');           % Segmentation directory
output_dir = "/data/tu_kilicb/DARTEL_registered_imgs/";          % Output directory

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Input files
b1_map = fullfile(b1_dir, 'MF2024-10-28_15-29-000140-00001-00002-0_B1map.nii');                 % Your B1 map
b1_ref = fullfile(b1_dir, 'MF2024-10-28_15-29-000140-00001-00002-0_B1map.nii');     % Anatomical B1 reference
r1_map = fullfile(quant_dir, 'MF2024-10-28_15-29-000140-00001-00001-0_R1.nii');             % R1 map (target resolution)
r2_map = fullfile(quant_dir, 'MF2024-10-28_15-29-000140-00001-00001-0_R2s_WLS1.nii');             % R2 map (target resolution)
pd_map = fullfile(quant_dir, 'MF2024-10-28_15-29-000140-00001-00001-0_R1.nii');             % PD map (target resolution)

% Tissue probability maps (segmentation results)
c1_file = fullfile(seg_dir, 'c1MF2024-10-28_15-29-000140-00001-00001-0_R1_outer_suppressed.nii');           % Grey matter probability
c2_file = fullfile(seg_dir, 'c2MF2024-10-28_15-29-000140-00001-00001-0_R1_outer_suppressed.nii');           % White matter probability
c3_file = fullfile(seg_dir, 'c3MF2024-10-28_15-29-000140-00001-00001-0_R1_outer_suppressed.nii');           % CSF probability


%% Input validation
fprintf('Validating input files...\n');

% Check if all required files exist
required_files = {b1_map, b1_ref, r1_map, r2_map, pd_map, c1_file, c2_file, c3_file};
file_names = {'B1 map', 'B1 reference', 'R1 map', 'R2 map', 'PD map', 'c1 tissue', 'c2 tissue', 'c3 tissue'};

for i = 1:length(required_files)
    if ~exist(required_files{i}, 'file')
        error('Required file not found: %s (%s)', file_names{i}, required_files{i});
    else
        fprintf('✓ Found %s: %s\n', file_names{i}, required_files{i});
    end
end

fprintf('All input files validated successfully.\n\n');
% Check if we have the required tissue class images
fprintf('Checking tissue class images...\n');

% Check if c1, c2 files exist
if ~exist(c1_file, 'file') || ~exist(c2_file, 'file')
    error('Tissue probability maps c1 and c2 not found. Please run SPM segmentation first.');
end

% DARTEL-imported tissue class images (these will be created by DARTEL)
rc1_file = fullfile(seg_dir, 'rc1_tissue.nii');
rc2_file = fullfile(seg_dir, 'rc2_tissue.nii');

fprintf('Found tissue class images: %s and %s\n', c1_file, c2_file);

%% Step 2: Create DARTEL Template and get deformation fields
% This step creates population-specific templates and generates deformation fields
template_dir = fullfile(output_dir, 'templates');
if ~exist(template_dir, 'dir')
    mkdir(template_dir);
end

fprintf('Creating DARTEL templates and deformation fields...\n');



spm('defaults', 'FMRI');
spm_jobman('initcfg');
% Job for creating DARTEL templates
% Note: This will create both templates and the rc1, rc2 files automatically
matlabbatch{1}.spm.tools.dartel.warp.images{1} = {c1_file};  % GM tissue class
matlabbatch{1}.spm.tools.dartel.warp.images{2} = {c2_file};  % WM tissue class
matlabbatch{1}.spm.tools.dartel.warp.settings.template = 'Template';
matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;

% Run DARTEL template creation
fprintf('Running DARTEL template creation...\n');

% Debug: Check if batch structure is valid
try
    spm_jobman('run', matlabbatch);
    fprintf('DARTEL template creation completed successfully.\n');
catch ME
    fprintf('Error in DARTEL template creation:\n');
    fprintf('Error message: %s\n', ME.message);
    fprintf('Error identifier: %s\n', ME.identifier);
    
    % Display batch structure for debugging
    fprintf('\nBatch structure:\n');
    disp(matlabbatch{1});
    
    rethrow(ME);
end
clear matlabbatch;

%% Step 3: Get the deformation field
% After DARTEL template creation, a deformation field is generated
% Look for the u_rc1_*.nii file (deformation field)
fprintf('Looking for deformation field...\n');

% Wait a moment for file system to update
pause(2);

% Look for deformation field with different possible naming patterns
deformation_field = spm_select('FPList', seg_dir, '^u_rc1.*\.nii');

%% Step 4: Normalize B1 map to MNI space using DARTEL deformation
fprintf('Normalizing B1 map using DARTEL deformation...\n');

% Get reference image (one of your quantitative maps) for output space
ref_img = r1_map;  % Use R1 map as reference for output space

% Job for normalizing B1 map to MNI space
matlabbatch{1}.spm.tools.dartel.mni_norm.template = {fullfile(spm('dir'), 'toolbox', 'dartel', 'icbm152_t1_tal_nlin_sym_09a.nii')};
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = {deformation_field};
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.images = {b1_map; b1_ref};
matlabbatch{1}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN];  % Keep original voxel size
matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];  % Keep original bounding box
matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;  % Don't preserve amounts (for non-modulated normalization)
matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];  % No smoothing

% Run normalization
spm_jobman('run', matlabbatch);
clear matlabbatch;

%% Step 5: Reslice normalized B1 map to match quantitative maps resolution
fprintf('Reslicing B1 map to match quantitative maps resolution...\n');

% Find the normalized B1 map (should have 'w' prefix)
normalized_b1 = spm_select('FPList', b1_dir, '^w.*B1_map.*\.nii$');
if isempty(normalized_b1)
    % Try alternative naming convention
    normalized_b1 = spm_select('FPList', b1_dir, '^w.*\.nii$');
end

if isempty(normalized_b1)
    error('Normalized B1 map not found. Check DARTEL normalization output.');
end

% Job for reslicing to match quantitative maps
matlabbatch{1}.spm.spatial.coreg.write.ref = {r1_map};  % Reference image (target resolution)
matlabbatch{1}.spm.spatial.coreg.write.source = {normalized_b1};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;  % 4th degree B-spline interpolation
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

% Run reslicing
spm_jobman('run', matlabbatch);
clear matlabbatch;

%% Step 6: Copy final registered B1 map to output directory
fprintf('Copying final registered B1 map to output directory...\n');

% Find the final registered B1 map (should have 'rw' prefix)
final_b1 = spm_select('FPList', b1_dir, '^rw.*B1_map.*\.nii');

%% Step 7: Verify registration quality
fprintf('Verifying registration quality...\n');

% Load images for comparison
V_b1 = spm_vol(output_b1);
V_r1 = spm_vol(r1_map);

% Check if dimensions match
if ~isequal(V_b1.dim, V_r1.dim)
    warning('Dimensions do not match perfectly:');
    fprintf('B1 map dimensions: %d x %d x %d\n', V_b1.dim);
    fprintf('R1 map dimensions: %d x %d x %d\n', V_r1.dim);
end

% Check if voxel sizes match (within tolerance)
vox_tol = 0.01;  % 0.01 mm tolerance
if max(abs(V_b1.mat(:) - V_r1.mat(:))) > vox_tol
    warning('Voxel sizes/orientations do not match perfectly');
    fprintf('B1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_b1.mat(1:3,1:3).^2)));
    fprintf('R1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_r1.mat(1:3,1:3).^2)));
else
    fprintf('Registration successful! B1 map now matches quantitative maps resolution.\n');
end

%% Step 8: Optional - Create a simple overlay for visual inspection
fprintf('Creating overlay for visual inspection...\n');

% Read image data
Y_b1 = spm_read_vols(V_b1);
Y_r1 = spm_read_vols(V_r1);

% Create a simple overlay (middle slice)
mid_slice = round(size(Y_b1, 3) / 2);
figure;
subplot(1,3,1);
imagesc(Y_b1(:,:,mid_slice)); colormap(gray); title('Registered B1 Map');
axis image; colorbar;

subplot(1,3,2);
imagesc(Y_r1(:,:,mid_slice)); colormap(gray); title('R1 Map (Reference)');
axis image; colorbar;

subplot(1,3,3);
imagesc(Y_b1(:,:,mid_slice) + Y_r1(:,:,mid_slice)); colormap(gray); title('Overlay');
axis image; colorbar;

sgtitle('Registration Quality Check');

fprintf('\nRegistration complete!\n');
fprintf('Output file: %s\n', output_b1);
fprintf('Please visually inspect the overlay to verify registration quality.\n');

if isempty(final_b1)
    final_b1 = spm_select('FPList', b1_dir, '^rw.*\.nii');

%% Step 7: Verify registration quality
fprintf('Verifying registration quality...\n');

% Load images for comparison
V_b1 = spm_vol(output_b1);
V_r1 = spm_vol(r1_map);

% Check if dimensions match
if ~isequal(V_b1.dim, V_r1.dim)
    warning('Dimensions do not match perfectly:');
    fprintf('B1 map dimensions: %d x %d x %d\n', V_b1.dim);
    fprintf('R1 map dimensions: %d x %d x %d\n', V_r1.dim);
end

% Check if voxel sizes match (within tolerance)
vox_tol = 0.01;  % 0.01 mm tolerance
if max(abs(V_b1.mat(:) - V_r1.mat(:))) > vox_tol
    warning('Voxel sizes/orientations do not match perfectly');
    fprintf('B1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_b1.mat(1:3,1:3).^2)));
    fprintf('R1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_r1.mat(1:3,1:3).^2)));
else
    fprintf('Registration successful! B1 map now matches quantitative maps resolution.\n');
end

%% Step 8: Optional - Create a simple overlay for visual inspection
fprintf('Creating overlay for visual inspection...\n');

% Read image data
Y_b1 = spm_read_vols(V_b1);
Y_r1 = spm_read_vols(V_r1);

% Create a simple overlay (middle slice)
mid_slice = round(size(Y_b1, 3) / 2);
figure;
subplot(1,3,1);
imagesc(Y_b1(:,:,mid_slice)); colormap(gray); title('Registered B1 Map');
axis image; colorbar;

subplot(1,3,2);
imagesc(Y_r1(:,:,mid_slice)); colormap(gray); title('R1 Map (Reference)');
axis image; colorbar;

subplot(1,3,3);
imagesc(Y_b1(:,:,mid_slice) + Y_r1(:,:,mid_slice)); colormap(gray); title('Overlay');
axis image; colorbar;

sgtitle('Registration Quality Check');

fprintf('\nRegistration complete!\n');
fprintf('Output file: %s\n', output_b1);
fprintf('Please visually inspect the overlay to verify registration quality.\n');

end

if ~isempty(final_b1)
    output_b1 = fullfile(output_dir, 'B1_map_registered.nii');
    copyfile(final_b1, output_b1);
    fprintf('Registered B1 map saved as: %s\n', output_b1);
else
    error('Final registered B1 map not found.');
end

%% Step 7: Verify registration quality
fprintf('Verifying registration quality...\n');

% Load images for comparison
V_b1 = spm_vol(output_b1);
V_r1 = spm_vol(r1_map);

% Check if dimensions match
if ~isequal(V_b1.dim, V_r1.dim)
    warning('Dimensions do not match perfectly:');
    fprintf('B1 map dimensions: %d x %d x %d\n', V_b1.dim);
    fprintf('R1 map dimensions: %d x %d x %d\n', V_r1.dim);
end

% Check if voxel sizes match (within tolerance)
vox_tol = 0.01;  % 0.01 mm tolerance
if max(abs(V_b1.mat(:) - V_r1.mat(:))) > vox_tol
    warning('Voxel sizes/orientations do not match perfectly');
    fprintf('B1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_b1.mat(1:3,1:3).^2)));
    fprintf('R1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_r1.mat(1:3,1:3).^2)));
else
    fprintf('Registration successful! B1 map now matches quantitative maps resolution.\n');
end

%% Step 8: Optional - Create a simple overlay for visual inspection
fprintf('Creating overlay for visual inspection...\n');

% Read image data
Y_b1 = spm_read_vols(V_b1);
Y_r1 = spm_read_vols(V_r1);

% Create a simple overlay (middle slice)
mid_slice = round(size(Y_b1, 3) / 2);
figure;
subplot(1,3,1);
imagesc(Y_b1(:,:,mid_slice)); colormap(gray); title('Registered B1 Map');
axis image; colorbar;

subplot(1,3,2);
imagesc(Y_r1(:,:,mid_slice)); colormap(gray); title('R1 Map (Reference)');
axis image; colorbar;

subplot(1,3,3);
imagesc(Y_b1(:,:,mid_slice) + Y_r1(:,:,mid_slice)); colormap(gray); title('Overlay');
axis image; colorbar;

sgtitle('Registration Quality Check');

fprintf('\nRegistration complete!\n');
fprintf('Output file: %s\n', output_b1);
fprintf('Please visually inspect the overlay to verify registration quality.\n');

if isempty(deformation_field)
    % Try alternative pattern
    deformation_field = spm_select('FPList', seg_dir, '^u_.*\.nii');

%% Step 4: Normalize B1 map to MNI space using DARTEL deformation
fprintf('Normalizing B1 map using DARTEL deformation...\n');

% Get reference image (one of your quantitative maps) for output space
ref_img = r1_map;  % Use R1 map as reference for output space

% Job for normalizing B1 map to MNI space
matlabbatch{1}.spm.tools.dartel.mni_norm.template = {fullfile(spm('dir'), 'toolbox', 'dartel', 'icbm152_t1_tal_nlin_sym_09a.nii')};
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = {deformation_field};
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.images = {b1_map; b1_ref};
matlabbatch{1}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN];  % Keep original voxel size
matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];  % Keep original bounding box
matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;  % Don't preserve amounts (for non-modulated normalization)
matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];  % No smoothing

% Run normalization
spm_jobman('run', matlabbatch);
clear matlabbatch;

%% Step 5: Reslice normalized B1 map to match quantitative maps resolution
fprintf('Reslicing B1 map to match quantitative maps resolution...\n');

% Find the normalized B1 map (should have 'w' prefix)
normalized_b1 = spm_select('FPList', b1_dir, '^w.*B1_map.*\.nii$');
if isempty(normalized_b1)
    % Try alternative naming convention
    normalized_b1 = spm_select('FPList', b1_dir, '^w.*\.nii$');
end

if isempty(normalized_b1)
    error('Normalized B1 map not found. Check DARTEL normalization output.');
end

% Job for reslicing to match quantitative maps
matlabbatch{1}.spm.spatial.coreg.write.ref = {r1_map};  % Reference image (target resolution)
matlabbatch{1}.spm.spatial.coreg.write.source = {normalized_b1};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;  % 4th degree B-spline interpolation
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

% Run reslicing
spm_jobman('run', matlabbatch);
clear matlabbatch;

%% Step 6: Copy final registered B1 map to output directory
fprintf('Copying final registered B1 map to output directory...\n');

% Find the final registered B1 map (should have 'rw' prefix)
final_b1 = spm_select('FPList', b1_dir, '^rw.*B1_map.*\.nii');

%% Step 7: Verify registration quality
fprintf('Verifying registration quality...\n');

% Load images for comparison
V_b1 = spm_vol(output_b1);
V_r1 = spm_vol(r1_map);

% Check if dimensions match
if ~isequal(V_b1.dim, V_r1.dim)
    warning('Dimensions do not match perfectly:');
    fprintf('B1 map dimensions: %d x %d x %d\n', V_b1.dim);
    fprintf('R1 map dimensions: %d x %d x %d\n', V_r1.dim);
end

% Check if voxel sizes match (within tolerance)
vox_tol = 0.01;  % 0.01 mm tolerance
if max(abs(V_b1.mat(:) - V_r1.mat(:))) > vox_tol
    warning('Voxel sizes/orientations do not match perfectly');
    fprintf('B1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_b1.mat(1:3,1:3).^2)));
    fprintf('R1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_r1.mat(1:3,1:3).^2)));
else
    fprintf('Registration successful! B1 map now matches quantitative maps resolution.\n');
end

%% Step 8: Optional - Create a simple overlay for visual inspection
fprintf('Creating overlay for visual inspection...\n');

% Read image data
Y_b1 = spm_read_vols(V_b1);
Y_r1 = spm_read_vols(V_r1);

% Create a simple overlay (middle slice)
mid_slice = round(size(Y_b1, 3) / 2);
figure;
subplot(1,3,1);
imagesc(Y_b1(:,:,mid_slice)); colormap(gray); title('Registered B1 Map');
axis image; colorbar;

subplot(1,3,2);
imagesc(Y_r1(:,:,mid_slice)); colormap(gray); title('R1 Map (Reference)');
axis image; colorbar;

subplot(1,3,3);
imagesc(Y_b1(:,:,mid_slice) + Y_r1(:,:,mid_slice)); colormap(gray); title('Overlay');
axis image; colorbar;

sgtitle('Registration Quality Check');

fprintf('\nRegistration complete!\n');
fprintf('Output file: %s\n', output_b1);
fprintf('Please visually inspect the overlay to verify registration quality.\n');

if isempty(final_b1)
    final_b1 = spm_select('FPList', b1_dir, '^rw.*\.nii');

%% Step 7: Verify registration quality
fprintf('Verifying registration quality...\n');

% Load images for comparison
V_b1 = spm_vol(output_b1);
V_r1 = spm_vol(r1_map);

% Check if dimensions match
if ~isequal(V_b1.dim, V_r1.dim)
    warning('Dimensions do not match perfectly:');
    fprintf('B1 map dimensions: %d x %d x %d\n', V_b1.dim);
    fprintf('R1 map dimensions: %d x %d x %d\n', V_r1.dim);
end

% Check if voxel sizes match (within tolerance)
vox_tol = 0.01;  % 0.01 mm tolerance
if max(abs(V_b1.mat(:) - V_r1.mat(:))) > vox_tol
    warning('Voxel sizes/orientations do not match perfectly');
    fprintf('B1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_b1.mat(1:3,1:3).^2)));
    fprintf('R1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_r1.mat(1:3,1:3).^2)));
else
    fprintf('Registration successful! B1 map now matches quantitative maps resolution.\n');
end

%% Step 8: Optional - Create a simple overlay for visual inspection
fprintf('Creating overlay for visual inspection...\n');

% Read image data
Y_b1 = spm_read_vols(V_b1);
Y_r1 = spm_read_vols(V_r1);

% Create a simple overlay (middle slice)
mid_slice = round(size(Y_b1, 3) / 2);
figure;
subplot(1,3,1);
imagesc(Y_b1(:,:,mid_slice)); colormap(gray); title('Registered B1 Map');
axis image; colorbar;

subplot(1,3,2);
imagesc(Y_r1(:,:,mid_slice)); colormap(gray); title('R1 Map (Reference)');
axis image; colorbar;

subplot(1,3,3);
imagesc(Y_b1(:,:,mid_slice) + Y_r1(:,:,mid_slice)); colormap(gray); title('Overlay');
axis image; colorbar;

sgtitle('Registration Quality Check');

fprintf('\nRegistration complete!\n');
fprintf('Output file: %s\n', output_b1);
fprintf('Please visually inspect the overlay to verify registration quality.\n');

end

if ~isempty(final_b1)
    output_b1 = fullfile(output_dir, 'B1_map_registered.nii');
    copyfile(final_b1, output_b1);
    fprintf('Registered B1 map saved as: %s\n', output_b1);
else
    error('Final registered B1 map not found.');
end

%% Step 7: Verify registration quality
fprintf('Verifying registration quality...\n');

% Load images for comparison
V_b1 = spm_vol(output_b1);
V_r1 = spm_vol(r1_map);

% Check if dimensions match
if ~isequal(V_b1.dim, V_r1.dim)
    warning('Dimensions do not match perfectly:');
    fprintf('B1 map dimensions: %d x %d x %d\n', V_b1.dim);
    fprintf('R1 map dimensions: %d x %d x %d\n', V_r1.dim);
end

% Check if voxel sizes match (within tolerance)
vox_tol = 0.01;  % 0.01 mm tolerance
if max(abs(V_b1.mat(:) - V_r1.mat(:))) > vox_tol
    warning('Voxel sizes/orientations do not match perfectly');
    fprintf('B1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_b1.mat(1:3,1:3).^2)));
    fprintf('R1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_r1.mat(1:3,1:3).^2)));
else
    fprintf('Registration successful! B1 map now matches quantitative maps resolution.\n');
end

%% Step 8: Optional - Create a simple overlay for visual inspection
fprintf('Creating overlay for visual inspection...\n');

% Read image data
Y_b1 = spm_read_vols(V_b1);
Y_r1 = spm_read_vols(V_r1);

% Create a simple overlay (middle slice)
mid_slice = round(size(Y_b1, 3) / 2);
figure;
subplot(1,3,1);
imagesc(Y_b1(:,:,mid_slice)); colormap(gray); title('Registered B1 Map');
axis image; colorbar;

subplot(1,3,2);
imagesc(Y_r1(:,:,mid_slice)); colormap(gray); title('R1 Map (Reference)');
axis image; colorbar;

subplot(1,3,3);
imagesc(Y_b1(:,:,mid_slice) + Y_r1(:,:,mid_slice)); colormap(gray); title('Overlay');
axis image; colorbar;

sgtitle('Registration Quality Check');

fprintf('\nRegistration complete!\n');
fprintf('Output file: %s\n', output_b1);
fprintf('Please visually inspect the overlay to verify registration quality.\n');

end

if isempty(deformation_field)
    % List all files in seg_dir to help debug
    fprintf('Files in segmentation directory:\n');
    dir_contents = dir(seg_dir);
    for i = 1:length(dir_contents)
        if ~dir_contents(i).isdir
            fprintf('  %s\n', dir_contents(i).name);
        end
    end
    error('Deformation field not found. Make sure DARTEL template creation was successful.');
end

fprintf('Found deformation field: %s\n', deformation_field);

%% Step 4: Normalize B1 map to MNI space using DARTEL deformation
fprintf('Normalizing B1 map using DARTEL deformation...\n');

% Get reference image (one of your quantitative maps) for output space
ref_img = r1_map;  % Use R1 map as reference for output space

% Job for normalizing B1 map to MNI space
matlabbatch{1}.spm.tools.dartel.mni_norm.template = {fullfile(spm('dir'), 'toolbox', 'dartel', 'icbm152_t1_tal_nlin_sym_09a.nii')};
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = {deformation_field};
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.images = {b1_map; b1_ref};
matlabbatch{1}.spm.tools.dartel.mni_norm.vox = [NaN NaN NaN];  % Keep original voxel size
matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [NaN NaN NaN; NaN NaN NaN];  % Keep original bounding box
matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;  % Don't preserve amounts (for non-modulated normalization)
matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = [0 0 0];  % No smoothing

% Run normalization
spm_jobman('run', matlabbatch);
clear matlabbatch;

%% Step 5: Reslice normalized B1 map to match quantitative maps resolution
fprintf('Reslicing B1 map to match quantitative maps resolution...\n');

% Find the normalized B1 map (should have 'w' prefix)
normalized_b1 = spm_select('FPList', b1_dir, '^w.*B1_map.*\.nii$');
if isempty(normalized_b1)
    % Try alternative naming convention
    normalized_b1 = spm_select('FPList', b1_dir, '^w.*\.nii$');
end

if isempty(normalized_b1)
    error('Normalized B1 map not found. Check DARTEL normalization output.');
end

% Job for reslicing to match quantitative maps
matlabbatch{1}.spm.spatial.coreg.write.ref = {r1_map};  % Reference image (target resolution)
matlabbatch{1}.spm.spatial.coreg.write.source = {normalized_b1};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;  % 4th degree B-spline interpolation
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

% Run reslicing
spm_jobman('run', matlabbatch);
clear matlabbatch;

%% Step 6: Copy final registered B1 map to output directory
fprintf('Copying final registered B1 map to output directory...\n');

% Find the final registered B1 map (should have 'rw' prefix)
final_b1 = spm_select('FPList', b1_dir, '^rw.*B1_map.*\.nii');

%% Step 7: Verify registration quality
fprintf('Verifying registration quality...\n');

% Load images for comparison
V_b1 = spm_vol(output_b1);
V_r1 = spm_vol(r1_map);

% Check if dimensions match
if ~isequal(V_b1.dim, V_r1.dim)
    warning('Dimensions do not match perfectly:');
    fprintf('B1 map dimensions: %d x %d x %d\n', V_b1.dim);
    fprintf('R1 map dimensions: %d x %d x %d\n', V_r1.dim);
end

% Check if voxel sizes match (within tolerance)
vox_tol = 0.01;  % 0.01 mm tolerance
if max(abs(V_b1.mat(:) - V_r1.mat(:))) > vox_tol
    warning('Voxel sizes/orientations do not match perfectly');
    fprintf('B1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_b1.mat(1:3,1:3).^2)));
    fprintf('R1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_r1.mat(1:3,1:3).^2)));
else
    fprintf('Registration successful! B1 map now matches quantitative maps resolution.\n');
end

%% Step 8: Optional - Create a simple overlay for visual inspection
fprintf('Creating overlay for visual inspection...\n');

% Read image data
Y_b1 = spm_read_vols(V_b1);
Y_r1 = spm_read_vols(V_r1);

% Create a simple overlay (middle slice)
mid_slice = round(size(Y_b1, 3) / 2);
figure;
subplot(1,3,1);
imagesc(Y_b1(:,:,mid_slice)); colormap(gray); title('Registered B1 Map');
axis image; colorbar;

subplot(1,3,2);
imagesc(Y_r1(:,:,mid_slice)); colormap(gray); title('R1 Map (Reference)');
axis image; colorbar;

subplot(1,3,3);
imagesc(Y_b1(:,:,mid_slice) + Y_r1(:,:,mid_slice)); colormap(gray); title('Overlay');
axis image; colorbar;

sgtitle('Registration Quality Check');

fprintf('\nRegistration complete!\n');
fprintf('Output file: %s\n', output_b1);
fprintf('Please visually inspect the overlay to verify registration quality.\n');

if isempty(final_b1)
    final_b1 = spm_select('FPList', b1_dir, '^rw.*\.nii');

%% Step 7: Verify registration quality
fprintf('Verifying registration quality...\n');

% Load images for comparison
V_b1 = spm_vol(output_b1);
V_r1 = spm_vol(r1_map);

% Check if dimensions match
if ~isequal(V_b1.dim, V_r1.dim)
    warning('Dimensions do not match perfectly:');
    fprintf('B1 map dimensions: %d x %d x %d\n', V_b1.dim);
    fprintf('R1 map dimensions: %d x %d x %d\n', V_r1.dim);
end

% Check if voxel sizes match (within tolerance)
vox_tol = 0.01;  % 0.01 mm tolerance
if max(abs(V_b1.mat(:) - V_r1.mat(:))) > vox_tol
    warning('Voxel sizes/orientations do not match perfectly');
    fprintf('B1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_b1.mat(1:3,1:3).^2)));
    fprintf('R1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_r1.mat(1:3,1:3).^2)));
else
    fprintf('Registration successful! B1 map now matches quantitative maps resolution.\n');
end

%% Step 8: Optional - Create a simple overlay for visual inspection
fprintf('Creating overlay for visual inspection...\n');

% Read image data
Y_b1 = spm_read_vols(V_b1);
Y_r1 = spm_read_vols(V_r1);

% Create a simple overlay (middle slice)
mid_slice = round(size(Y_b1, 3) / 2);
figure;
subplot(1,3,1);
imagesc(Y_b1(:,:,mid_slice)); colormap(gray); title('Registered B1 Map');
axis image; colorbar;

subplot(1,3,2);
imagesc(Y_r1(:,:,mid_slice)); colormap(gray); title('R1 Map (Reference)');
axis image; colorbar;

subplot(1,3,3);
imagesc(Y_b1(:,:,mid_slice) + Y_r1(:,:,mid_slice)); colormap(gray); title('Overlay');
axis image; colorbar;

sgtitle('Registration Quality Check');

fprintf('\nRegistration complete!\n');
fprintf('Output file: %s\n', output_b1);
fprintf('Please visually inspect the overlay to verify registration quality.\n');

end

if ~isempty(final_b1)
    output_b1 = fullfile(output_dir, 'B1_map_registered.nii');
    copyfile(final_b1, output_b1);
    fprintf('Registered B1 map saved as: %s\n', output_b1);
else
    error('Final registered B1 map not found.');
end

%% Step 7: Verify registration quality
fprintf('Verifying registration quality...\n');

% Load images for comparison
V_b1 = spm_vol(output_b1);
V_r1 = spm_vol(r1_map);

% Check if dimensions match
if ~isequal(V_b1.dim, V_r1.dim)
    warning('Dimensions do not match perfectly:');
    fprintf('B1 map dimensions: %d x %d x %d\n', V_b1.dim);
    fprintf('R1 map dimensions: %d x %d x %d\n', V_r1.dim);
end

% Check if voxel sizes match (within tolerance)
vox_tol = 0.01;  % 0.01 mm tolerance
if max(abs(V_b1.mat(:) - V_r1.mat(:))) > vox_tol
    warning('Voxel sizes/orientations do not match perfectly');
    fprintf('B1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_b1.mat(1:3,1:3).^2)));
    fprintf('R1 voxel size: [%.3f %.3f %.3f]\n', sqrt(sum(V_r1.mat(1:3,1:3).^2)));
else
    fprintf('Registration successful! B1 map now matches quantitative maps resolution.\n');
end

%% Step 8: Optional - Create a simple overlay for visual inspection
fprintf('Creating overlay for visual inspection...\n');

% Read image data
Y_b1 = spm_read_vols(V_b1);
Y_r1 = spm_read_vols(V_r1);

% Create a simple overlay (middle slice)
mid_slice = round(size(Y_b1, 3) / 2);
figure;
subplot(1,3,1);
imagesc(Y_b1(:,:,mid_slice)); colormap(gray); title('Registered B1 Map');
axis image; colorbar;

subplot(1,3,2);
imagesc(Y_r1(:,:,mid_slice)); colormap(gray); title('R1 Map (Reference)');
axis image; colorbar;

subplot(1,3,3);
imagesc(Y_b1(:,:,mid_slice) + Y_r1(:,:,mid_slice)); colormap(gray); title('Overlay');
axis image; colorbar;

sgtitle('Registration Quality Check');

fprintf('\nRegistration complete!\n');
fprintf('Output file: %s\n', output_b1);
fprintf('Please visually inspect the overlay to verify registration quality.\n');
