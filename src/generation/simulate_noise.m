function kspace_noisy = simulate_noise(kspace, twix_file_path, cfg)
%ADD_AND_CORRELATE_NOISE Adds noise and correlation to k-space data.
%
%   Args:
%       kspace: K-space data of size [nx, n_readout, ny, n_coil].
%       twix_file_path: Path to the TWIX file containing noise scans.
%       SNR_dB: Desired SNR in decibels.
%
%   Returns:
%       kspace_noisy: K-space data with added Gaussian and correlated noise.

    output_dir = cfg.paths.output_dir;
    SNR_dB = cfg.params.noise_params.SNR_dB;
    if SNR_dB == -1 % This corresponds to no noise
        kspace_noisy = kspace;
        return;
    end
    
    psi_sqrt_path = fullfile(output_dir, 'synthetic_kspace/corrolated_noise_matrix.mat');
    if exist(psi_sqrt_path, "file")
        psi_sqrt = load(psi_sqrt_path);
    else
    % 1. Get Noise Scans from TWIX file to calculate the correlation matrix
    noise_output_path = fullfile(output_dir, 'synthetic_kspace/noise_scan.mat');
    try
        noise_scans = read_twix(twix_file_path, noise_output_path, 'modality', 'noise'); % Shape: (n_readout x n_coils)
    catch
        error('Could not read noise scans from TWIX file. Cannot calculate noise correlation.');
    end
    
    
    % 1. Calculate Noise Covariance Matrix (Psi) and its square root
    noise_covariance_matrix = cov(noise_scans);
    
    % For numerical stability, add a small regularization term
    noise_covariance_matrix = noise_covariance_matrix + 1e-6 * eye(size(noise_scans,2));
    
    % Calculate the square root of the covariance matrix (Psi^1/2)
    [V, D] = eig(noise_covariance_matrix);
    psi_sqrt = V * sqrt(D) * V';    
    save(psi_sqrt_path, 'psi_sqrt');
    end
    
    % Convert SNR from dB to linear scale
    SNR_linear = 10^(SNR_dB / 10);
    
    % Compute signal power
    signal_power = mean(abs(kspace(:)).^2);
    
    % Compute noise power needed
    noise_power = signal_power / SNR_linear;
    
    % Standard deviation of the complex Gaussian noise
    sigma = sqrt(noise_power / 2);
    
    % Generate i.i.d. complex Gaussian noise
    noise_uncorrelated = (sigma * randn(size(kspace)) + 1i * sigma * randn(size(kspace)));
    
    % 3. Apply Noise Correlation and Add to K-space
    % Reshape arrays for matrix multiplication
    [nx, n_readout, ny, n_coils] = size(kspace);
    noise_reshaped = reshape(noise_uncorrelated, [], n_coils);
    
    if cfg.params.noise_params.add_corrolated_noise
        % Apply the correlation matrix to the noise
        noise_correlated_reshaped = noise_reshaped * psi_sqrt;
    else

        noise_correlated_reshaped = noise_reshaped;
    end

    % Reshape correlated noise back to original dimensions
    noise_correlated = reshape(noise_correlated_reshaped, nx, n_readout, ny, n_coils);
    
    % Add the correlated noise to the original k-space
    kspace_noisy = kspace + noise_correlated;
    
end
