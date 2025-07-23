function file_struct = load_mri_data(file_struct)
%LOAD_MRI_DATA Augment input struct with NIfTI data and header
%   file_struct = LOAD_MRI_DATA(file_struct)
%   Input: file_struct with modality fields (e.g., R1, R2, PD...) each having
%          a .filepath property pointing to a NIfTI file (string).
%   Output: same struct, where for each modality with a valid filepath:
%           .data       - image array (scaled),
%           .ref_header - NIfTI header struct from niftiinfo.

    modalities = fieldnames(file_struct);
    for i = 1:numel(modalities)
        mod = modalities{i};
        % Skip if no filepath field or empty
        if ~isfield(file_struct.(mod), 'filepath') || isempty(file_struct.(mod).filepath)
            fprintf('Skipping %s: no valid filepath provided.\n', mod);
            continue;
        end
        fp = file_struct.(mod).filepath;
        if ~exist(fp, 'file')
            warning('File not found for %s: %s. Skipping.\n', mod, fp);
            continue;
        end
        try
            % Read image and header
            nii_data = niftiread(fp);
            nii_info = niftiinfo(fp);
            % Apply scaling if available
            if isfield(nii_info, 'raw') && isfield(nii_info.raw, 'scl_slope')
                slope = nii_info.raw.scl_slope;
                intercept = nii_info.raw.scl_inter;
            else
                slope = 1;
                intercept = 0;
            end
            % Store in struct
            file_struct.(mod).data = double(nii_data) * slope + intercept;
            file_struct.(mod).ref_header = nii_info;
            fprintf('Loaded %s: %s\n', mod, fp);
        catch ME
            warning('Error loading %s from %s: %s', mod, fp, ME.message);
        end
    end
end
