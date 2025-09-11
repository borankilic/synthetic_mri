%gtPath = '/data/tu_kilicb/hMRI-toolbox-demodata-main/raw-data/t1w_mfc_3dflash_v1i_R4_0015/anon_s2018-02-28_18-26-190921-00001-00224-1.nii';
%synthPath = '/data/tu_kilicb/synthetic_demodata_tissue_map_mask/synthetic_images/anon_s2018-02-28_18-26-185345-00001-00224-1_GRE_synthim_TE2.200000e+00_FA20_TR1.870000e+01.nii';
function metrics = synthImageEval(gtPath, synthPath)
% Compare ground truth and synthetic NIfTI images using various similarity metrics
% Returns struct with metrics: linear_fit, nrmse, pearson_r, ssim, mutual_info

    % Load images
    [I_gt, I_synth] = load_nifti_images(gtPath, synthPath);
    
    % Coregister synthetic to ground truth using FSL
    I_synth_reg = coregister_fsl(I_synth, I_gt);
    
    % Create mask and extract values
    [gt_vals, synth_vals] = extract_masked_values(I_gt, I_synth_reg);
    
    % Compute metrics
    metrics = compute_metrics(gt_vals, synth_vals, I_gt, I_synth_reg);
    
    % Display results
    display_results(I_gt, I_synth_reg, metrics);
end

function [I_gt, I_synth] = load_nifti_images(gtPath, synthPath)
    I_gt = double(niftiread(gtPath));
    I_synth = double(niftiread(synthPath));
    
    % Handle 4D images
    if ndims(I_gt) == 4, I_gt = squeeze(I_gt(:,:,:,1)); end
    if ndims(I_synth) == 4, I_synth = squeeze(I_synth(:,:,:,1)); end
end

function I_synth_reg = coregister_fsl(I_synth, I_gt)
    % Save temporary files for FSL
    temp_gt = 'temp_gt.nii.gz';
    temp_synth = 'temp_synth.nii.gz';
    temp_reg = 'temp_reg.nii.gz';
    temp_mat = 'temp_transform.mat';
    
    niftiwrite(I_gt, temp_gt);
    niftiwrite(I_synth, temp_synth);
    
    % FSL rigid registration
    cmd = sprintf('FSL flirt -in %s -ref %s -out %s -omat %s -dof 6', ...
        temp_synth, temp_gt, temp_reg, temp_mat);
    
    if system(cmd) == 0
        I_synth_reg = double(niftiread(temp_reg));
    else
        warning('FSL registration failed, using original image');
        I_synth_reg = I_synth;
    end
    
    % Cleanup
    delete(temp_gt); delete(temp_synth); delete(temp_reg); delete(temp_mat);
end

function [gt_vals, synth_vals] = extract_masked_values(I_gt, I_synth_reg)
    % Create initial mask from nonzero synthetic voxels
    mask = (I_synth_reg ~= 0) & isfinite(I_synth_reg);
    
    % Apply hole filling using FSL
    mask = holefill_fsl(mask);
    
    gt_vals = I_gt(mask);
    synth_vals = I_synth_reg(mask);
    
    % Remove invalid pairs
    valid = isfinite(gt_vals) & isfinite(synth_vals);
    gt_vals = gt_vals(valid);
    synth_vals = synth_vals(valid);
end

function mask_filled = holefill_fsl(mask)
    % Use FSL to fill holes in binary mask
    temp_mask = 'temp_mask.nii.gz';
    temp_filled = 'temp_mask_filled.nii.gz';
    
    % Save mask as NIfTI
    niftiwrite(double(mask), temp_mask);
    
    % FSL hole filling command
    cmd = sprintf('FSL fslmaths %s -fillh %s', temp_mask, temp_filled);
    
    if system(cmd) == 0
        mask_filled = logical(niftiread(temp_filled));
    else
        warning('FSL hole filling failed, using original mask');
        mask_filled = mask;
    end
    
    % Cleanup
    delete(temp_mask); 
    delete(temp_filled);
end

function metrics = compute_metrics(gt_vals, synth_vals, I_gt, I_synth_reg)
    % Linear fit
    coef = [gt_vals(:), ones(length(gt_vals),1)] \ synth_vals(:);
    metrics.linear_fit.a = coef(1);
    metrics.linear_fit.b = coef(2);
    
    % NRMSE
    residuals = synth_vals - (coef(1) * gt_vals + coef(2));
    rmse = sqrt(mean(residuals.^2));
    metrics.nrmse = rmse / (max(gt_vals) - min(gt_vals));
    
    % Pearson correlation
    R = corrcoef(gt_vals, synth_vals);
    metrics.pearson_r = R(1,2);
    
    % SSIM (slice-wise average)
    metrics.ssim = compute_ssim_3d(I_gt, I_synth_reg);
    
    % Mutual information
    metrics.mutual_info = mutual_info(gt_vals, synth_vals, 64);
    
    metrics.n_voxels = length(gt_vals);
end

function mean_ssim = compute_ssim_3d(I_gt, I_synth_reg)
    mask = (I_synth_reg ~= 0);
    ssim_vals = [];
    
    for z = 1:size(I_gt,3)
        mask_slice = mask(:,:,z);
        if nnz(mask_slice) < 100, continue; end
        
        gt_slice = I_gt(:,:,z);
        synth_slice = I_synth_reg(:,:,z);
        
        % Normalize to [0,1]
        gt_n = normalize_slice(gt_slice, mask_slice);
        synth_n = normalize_slice(synth_slice, mask_slice);
        
        gt_n(~mask_slice) = 0;
        synth_n(~mask_slice) = 0;
        
        try
            ssim_vals(end+1) = ssim(gt_n, synth_n);
        catch
            continue;
        end
    end
    
    mean_ssim = mean(ssim_vals);
end

function img_n = normalize_slice(img, mask)
    vals = img(mask);
    p_lo = prctile(vals, 1);
    p_hi = prctile(vals, 99);
    img_n = (img - p_lo) / max(1e-12, p_hi - p_lo);
    img_n = max(0, min(1, img_n));
end

function mi = mutual_info(x, y, nbins)
    xed = linspace(min(x), max(x), nbins+1);
    yed = linspace(min(y), max(y), nbins+1);
    
    H = histcounts2(x, y, xed, yed) + eps;
    Pxy = H / sum(H(:));
    Px = sum(Pxy, 2);
    Py = sum(Pxy, 1);
    
    Hx = -sum(Px .* log(Px));
    Hy = -sum(Py .* log(Py));
    Hxy = -sum(Pxy(:) .* log(Pxy(:)));
    
    mi = Hx + Hy - Hxy;
end

function display_results(I_gt, I_synth_reg, metrics)
    fprintf('Linear fit: a=%.3g, b=%.3g\n', metrics.linear_fit.a, metrics.linear_fit.b);
    fprintf('NRMSE: %.3g | Pearson r: %.3g\n', metrics.nrmse, metrics.pearson_r);
    fprintf('SSIM: %.3g | MI: %.3g\n', metrics.ssim, metrics.mutual_info);
    
    % Display center slice (sagittal view)
    cz = round(size(I_gt,2)/2+50);
    mask_slice = (squeeze(I_synth_reg(:,cz,:)) ~= 0);
    
    slice_gt = squeeze(I_gt(:,cz,:));
    slice_synth = squeeze(I_synth_reg(:,cz,:));
    slice_gt(~mask_slice) = 0;
    
    % Scale synthetic image using linear regression parameters
    slice_synth_scaled = metrics.linear_fit.a * slice_synth + metrics.linear_fit.b;
    
    figure('Color','w');
    subplot(1,2,1); 
    imagesc(mat2gray(slice_gt)); 
    axis image off; 
    title('GT (masked)');
    colormap(gca, 'gray');
    
    subplot(1,2,2); 
    imagesc(mat2gray(slice_synth_scaled)); 
    axis image off; 
    title('Synthetic (scaled)');
    colormap(gca, 'gray');
end