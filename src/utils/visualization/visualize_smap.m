if ~exist('smoothed_smaps', 'var') || ~exist('synth_image', 'var')
    smoothed_smaps = niftiread('/data/tu_kilicb/synthetic_imgs_generated_fsl_mask/sens_maps/MF2024-10-28_15-29-000140-00001-00001-0_smap_mag_gaussian_smoothed_s25.nii.gz');
    synth_image = niftiread('/data/tu_kilicb/synthetic_imgs_generated_fsl_mask/MF2024-10-28_15-29-000140-00001-00001-0_synthim_TE4_FA10_TR25.nii');
end

fprintf('smoothed_smaps size: %s\n', mat2str(size(smoothed_smaps)));
fprintf('synth_image size: %s\n', mat2str(size(synth_image)));

% Match synth_image to smoothed_smaps dimensions
synth_image = synth_image(1:size(smoothed_smaps,1), ...
                          1:size(smoothed_smaps,2), ...
                          1:size(smoothed_smaps,3));

% Combine
sens_incorp = smoothed_smaps(:,210,:,:) .* synth_image(:,210,:,:);

% Remove singleton dimension → size: 335 x 288 x 32
sens_incorp = squeeze(sens_incorp)/100;
sens_incorp = permute(flip(sens_incorp,1), [2,1,3]);
% Visualize all 32 channels
figure;
imshow3(sens_incorp,[], [4,8])

