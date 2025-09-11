function register_b1_to_ref(data_struct , ref_path, output_dir,interp_method)
   
    
    try
        b1_map_path=data_struct.B1map.filepath;
        b1_ref_path=data_struct.B1ref.filepath;
    catch
        fprintf("No B1 map found. Skipping registeration.\n")
        return
    end
    [~, name, ext] = fileparts(b1_map_path);
    if strcmp(ext, '.gz') && endsWith(name, '.nii')
        name = extractBefore(name, '.nii');
    end
    output_path = fullfile(output_dir, [name '_registered.nii.gz']);
    if exist(output_path,'file')
        fprintf("Registered B1 map already exists.")
        return
    end
    matfile = tempname + ".mat";
    cmd1 = sprintf(['FSL flirt -in "%s" -ref "%s" -out "%s" -omat "%s" ' ...
                    '-interp "%s" -dof 6'], ...
                    b1_ref_path, ref_path, output_path, matfile, interp_method);
    cmd2 = sprintf(['FSL flirt -in "%s" -applyxfm -init "%s" -ref "%s" ' ...
                    '-out "%s" -interp "%s"'], ...
                    b1_map_path, matfile, ref_path, output_path, interp_method);
    if ~isfile(output_path)
        [status1, result1] = system(cmd1);
        if status1 ~= 0
            error('FLIRT registration failed:\n%s', result1);
        end
        [status2, result2] = system(cmd2);
        if status2 == 0
            fprintf('Registration and resampling complete. Output: %s\n', output_path);
        else
            error('FLIRT resampling failed:\n%s', result2);
        end
    else 
        fprintf('Registered B1 %s already exists. Skipping registeration ', output_path); 
    end
end
