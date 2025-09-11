function [status, cmdout] = runRomeo(cfg)
% runRomeoWithMerge  Merge NIfTI files (if needed) and call romeo.jl
%
% [status, cmdout] = runRomeoWithMerge(magFolder, phaseFolder, echoTimes, B0_folder, romeoScriptPath)
%
% Inputs:
%   magFolder       - path to folder containing magnitude .nii files
%   phaseFolder     - path to folder containing phase .nii files
%   echoTimes       - numeric vector of echo times, e.g. [2.1,4.2,6.3]
%   B0_folder      - path for romeo output
%   romeoScriptPath - full path to romeo.jl (or executable romeo script)
%
% Requires: FSL installed and 'fslmerge' available in system PATH
    magFolder = cfg.paths.im_mag_dir;
    phaseFolder = cfg.paths.im_phase_dir;
    echoTimes = cfg.params.real_echo_times;
    output_path = cfg.paths.output_dir;
    mpm_dir = cfg.paths.patient_dir;
    romeoScriptPath = cfg.paths.romeo_script_path;


    
    B0_folder = fullfile(output_path, 'B0_map_unwrapped');
    % ensure output folder exists
    if ~exist(B0_folder, 'dir')
        mkdir(B0_folder);
    end

    % output merged file names
    mergedMag = fullfile(B0_folder, 'merged_mag.nii.gz');
    mergedPhase = fullfile(B0_folder, 'merged_phase.nii.gz');

    % merge magnitude if not existing
    if ~isfile(mergedMag)
        magFiles = dir(fullfile(magFolder, '*.nii'));
        if isempty(magFiles)
            error('No magnitude .nii files found in %s', magFolder);
        end
        magList = strjoin(arrayfun(@(f) fullfile(magFolder, f.name), magFiles, 'UniformOutput', false), ' ');
        cmdMagMerge = sprintf('FSL fslmerge -t "%s" %s', mergedMag, magList);
        fprintf('Merging magnitude images:\n%s\n', cmdMagMerge);
        system(cmdMagMerge);
    else
        fprintf('Merged magnitude image already exists. Skipping merge.\n');
    end

    % merge phase if not existing
    if ~isfile(mergedPhase)
        phaseFiles = dir(fullfile(phaseFolder, '*.nii'));
        if isempty(phaseFiles)
            error('No phase .nii files found in %s', phaseFolder);
        end
        phaseList = strjoin(arrayfun(@(f) fullfile(phaseFolder, f.name), phaseFiles, 'UniformOutput', false), ' ');
        cmdPhaseMerge = sprintf('FSL fslmerge -t "%s" %s', mergedPhase, phaseList);
        fprintf('Merging phase images:\n%s\n', cmdPhaseMerge);
        system(cmdPhaseMerge);
    else
        fprintf('Merged phase image already exists. Skipping merge.\n');
    end

    % format echo times
    etStr = ['[' strjoin(arrayfun(@(x) sprintf('%.6g', x), echoTimes, 'UniformOutput', false), ',') ']'];

    % build romeo command
    cmdRomeo = sprintf( './run_romeo.sh "%s" -p "%s" -m "%s" -t %s -B -o "%s"', ...
    romeoScriptPath, mergedPhase, mergedMag, etStr, B0_folder );
    
    if ~isfile(fullfile(B0_folder, "B0.nii"))
        fprintf('Running romeo:\n%s\n', cmdRomeo);
        [status, cmdout] = system(cmdRomeo);

        if status ~= 0
            warning('romeo exited with status %d. Output:\n%s', status, cmdout);
        else
            fprintf('romeo finished successfully.\n');
        end
    else
        fprintf('B0 map already exists, skipping calculation.\n')
    end
    


    %% Coregisteration of B0 maps onto synthetic images 
    % During th MPM creation, some rigid body transformation could have
    % been applied. Thus we coregister the B0 map using on of the weighted
    % magnitude images with synthetic images.
    magFiles = dir(fullfile(magFolder, '*.nii'));
    weighted_mag_filepath = fullfile(magFiles(1).folder, magFiles(1).name);
    mpmFiles = dir([mpm_dir 'Results/' '*R1.nii']);
    mpm_filepath = fullfile(mpmFiles(1).folder, mpmFiles(1).name);
    B0_file = fullfile(B0_folder, 'B0.nii');
    out_file = fullfile(B0_folder, 'B0_coreg.nii.gz');

    % Temporary transformation matrix file
    matfile = fullfile(B0_folder, 'transform_matrix.mat');

    % ensure geometry matches
    cmd_cpgeom = sprintf('FSL fslcpgeom "%s" "%s"', weighted_mag_filepath, B0_file);
    
    % estimate translation-only transform
    cmd1 = sprintf(['FSL flirt -in "%s" -ref "%s" -omat "%s" -dof 6 -cost normcorr -usesqform ' ...
                    '-searchrx 0 0 -searchry 0 0 -searchrz 0 0 -out "%s"'], ...
                   weighted_mag_filepath, mpm_filepath, matfile, fullfile(tempdir,'mag2synth.nii.gz'));
    
    % apply transform to B0 (header-matched)

    cmd2 = sprintf(['FSL flirt -in "%s" -applyxfm -init "%s" -ref "%s" -out "%s" -usesqform -interp trilinear'], ...
                   B0_file, matfile, mpm_filepath, out_file);


    if ~isfile(out_file)

        system(cmd_cpgeom);
        % Run rigid registration
        [status1, result1] = system(cmd1);
        if status1 ~= 0
            error('FLIRT registration failed:\n%s', result1);
        end
    
        % Apply transformation and resample
        [status2, result2] = system(cmd2);
        if status2 == 0
            fprintf('Registration and resampling complete. Output: %s\n', out_file);
        else
            error('FLIRT resampling failed:\n%s', result2);
        end
    else 
        fprintf('Registered B0 %s already exists. Skipping registeration ', out_file);
     end

end
