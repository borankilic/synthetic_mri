function nii = load_nii(filename)
%LOAD_NII Robustly load a NIfTI file using available toolboxes.
%   nii = load_nii(filename)
%   Returns a struct with fields:
%       nii.img  -> image array
%       nii.hdr  -> header/metadata
%
%   Fallback order:
%   1) External load_nii (if available on path, but not this file)
%   2) SPM (spm_vol + spm_read_vols)
%   3) MATLAB built-in (niftiread + niftiinfo, R2017b+)

    % --- 1) Try external load_nii (not this function) ---
    external = which('load_nii');
    thisFile = [mfilename('fullpath') '.m'];
    if ~isempty(external) && ~strcmp(external, thisFile)
        try
            nii = load_nii(filename); % calls external
            fprintf('Loaded using external load_nii -> %s\n', filename);
            return;
        catch ME
            warning('External load_nii failed: %s\nFalling back...', E.message);
        end
    end

    % --- 2) Try SPM ---
    if exist('spm_vol','file') == 2 && exist('spm_read_vols','file') == 2
        try
            V = spm_vol(filename);
            img = spm_read_vols(V);
            nii.img = img;
            nii.hdr = V;

            % Add some fields to mimic NIfTI toolbox style
            nii.hdr.dime.datatype = V.dt(1);
            nii.hdr.dime.bitpix   = spm_type(V.dt(1), 'bits');
            nii.hdr.dime.cal_max  = max(img(:));
            nii.hdr.dime.cal_min  = min(img(:));

            fprintf('Loaded using SPM -> %s\n', filename);
            return;
        catch ME
            warning('SPM load failed: %s\nFalling back...', E.message);
        end
    end

    % --- 3) Try MATLAB built-in ---
    if exist('niftiread','file') == 2 && exist('niftiinfo','file') == 2
        try
            img  = niftiread(filename);
            info = niftiinfo(filename);
            nii.img = img;
            nii.hdr = info;

            fprintf('Loaded using MATLAB niftiread -> %s\n', filename);
            return;
        catch ME
            warning('niftiread/niftiinfo failed: %s\n', E.message);
        end
    end

    % --- If none worked ---
    error('Unable to read NIfTI file. Please install NIfTI toolbox, SPM, or MATLAB R2017b+.');

end
