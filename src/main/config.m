function config = config()
    % Configuration for the synthetic MRI generation pipeline

    % --- Paths ---
    config.paths.patient_dir = '/data/tu_kilicb/derivatives/Run_03_from_demodata_hMRI_perfect_MPMs_custom_defaults/';
    config.paths.kspace_path = '/data/tu_kilicb/KP2BK/meas_MID00130_FID06609_t1w_kp_mtflash3d_v1s_0p65.dat';
    config.paths.im_mag_dir = '/data/tu_kilicb/KP2BK/20241028_27993.cb/DICOM_imports/MF27993_2x2_T1w_mag_S69/';
    config.paths.im_phase_dir = '/data/tu_kilicb/KP2BK/20241028_27993.cb/DICOM_imports/MF27993_2x2_T1w_phase_S70/';
    config.paths.output_dir  = '/data/tu_kilicb/synthetic_demodata_tissue_map_mask/';    
    config.paths.espirit_path = '/data/u_kilicb_software/ESPIRiT/';
    config.paths.mapvbvd_path = '../../scripts/mapVBVD.m';    %DO NOT MODIFY unless you know what you are doing
    config.paths.romeo_script_path = '../../scripts/romeo.jl';   

    % --- Sequence Parameters ---
    config.params.sequence = 'GRE';  %options are: 'GRE', 'MP2RAGE', 'FSE', 'FLAIR' 
    config.params.real_echo_times = [3,5.5,8,10.5,13,15.5]; % Can be obtained from the configuration file of the sequence.
    %config.params.TE_values = [3, 5.5, 8, 10.5, 13, 15.5];
    config.params.TE_values = [2.2];
    config.params.TR_values = [18.70];
    config.params.FA_values = [6, 20];
    config.params.TI_values = [1];  % For MP2RAGE each combination should be a cell representing first and second inversion recovery times. 
                                                                % MP2RAGE: TI_values = [ {10; 20}, {15; 20}, {15; 30}]
                                                                % FSE: TI_values = [10, 20, 30, 40]

   config.params.smap_params = struct('ncalib', 40,...                  % Size of the Auto Calibration Region (ACS) is taken as ncalib x calib (or if [nx, ny]
                                       'ksize', [5 5], ...              % Size of the interpolating kernel
                                       'eigThresh1', 0.1, ...           % It is used to only select eigenvectors with eigenvalues greater than eigThresh1 for the reconstruction matrix
                                       'eigThresh2', 0.8, ...           % It is used to only select eigenvectors with eigenvalues greater than eigThresh2 for the construction of sensititvity maps. If there is multiple eigenvectors with eigenvalue greater than eigThres2, resulting sens_map will be a linear combination of those. 
                                       'show_figures', 'none', ...      % % options: 'none','essential','detailed'
                                       'verbose', true, ...
                                       'smoothing_radius', 25);         % Radius of the gaussian smoothing kernel applied after ESPIRiT algorithm.                                              
                                                                
    
    config.params.kspace_params = struct('fft_shift' , true, ...        % Apply fftshift while taking FFT
                                         'normalize' , true, ...        % Normalize during FFT
                                         'verbose', true, ...
                                         'save_complex', true);         % Save magnitude and phase in seperate nifti files if true, otherwise just save the magnitude.
   
    config.params.noise_params =struct('SNR_dB', -1, ...                % Add noise based on SNR of resulting image. Set to -1 for no noise
                                        'add_correlated_noise',true,... % If true covariance matrix from noise scans are calculated and the added noise is correlated, if false awgn noise is added 
                                        'num_variations', 1);           % Change for generating multiple weighted images (with random noise distributions) )


    % --- Processing Options ---
    config.options.extract_brain = 'tissue_maps'; % 'bet', 'tissue_maps', or 'no'
    config.options.interp_method = 'trilinear';
    config.options.threshold_tmap = 0.2;    % Threshold for brain extraction using tissue maps. 
    config.options.threshold_bet = 0.5;     % Threshold used for BET. Fractional intensity threshold (0->1). Smaller values give larger brain outline estimates 
    config.options.signal_constant = 2000;   % Arbitrary signal constant representing the lumped multiplicative effects on the signal intensity    
    config.options.add_RF_bias = false; % Whether to reintroduce RF sensitivity bias to synthetic images. If MPMs are not corrected for bias, it's not recommended to apply it.
  

end