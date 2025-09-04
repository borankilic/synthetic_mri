function config = config()
    % Configuration for the synthetic MRI generation pipeline

    % --- Paths ---
    config.paths.patient_dir = '/data/tu_kilicb/KP2BK/20241028_27993.cb/mpm/';
    config.paths.kspace_path = '/data/tu_kilicb/KP2BK/meas_MID00130_FID06609_t1w_kp_mtflash3d_v1s_0p65.dat';
    config.paths.im_mag_dir = '/data/tu_kilicb/KP2BK/20241028_27993.cb/DICOM_imports/MF27993_2x2_T1w_mag_S69/';
    config.paths.im_phase_dir = '/data/tu_kilicb/KP2BK/20241028_27993.cb/DICOM_imports/MF27993_2x2_T1w_phase_S70/';
    config.paths.output_dir  = '/data/tu_kilicb/synthetic_imgs_generated_bet_mask/';    
    config.paths.espirit_path = '/data/u_kilicb_software/ESPIRiT/';
    config.paths.mapvbvd_path = '../../scripts/mapVBVD-main/';    %DO NOT MODIFY
    config.paths.romeo_script_path = '../../scripts/romeo.jl';    
    % --- Sequence Parameters ---
    config.params.real_echo_times = [3,5.5,8,10.5,13,15.5]; % Can be obtained from the configuration file of the sequence.
    config.params.TE_values = [3, 4, 5, 6, 7, 8, 9, 10];
    config.params.TR_values = [25];
    config.params.FA_values = [5, 10, 20, 30, 40];

    % --- Processing Options ---
    config.options.extract_brain = 'bet'; % 'bet', 'fsl', or 'no'
    config.options.interp_method = 'trilinear';
    config.options.signal_constant = 2000;
end