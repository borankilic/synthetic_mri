function signal = compute_ernst_signal(data_struct, TE, TR, alpha, signal_constant, interp_method, extract_brain,)
%COMPUTE_ERNST_SIGNAL Compute MRI signal using Ernst equation
% signal = compute_ernst_signal(data_struct, TE, TR, alpha, signal_constant, interp_method)
% computes the MRI signal using the Ernst equation:
%
% S = K * PD * sin(α*B1) * (1 - exp(-TR*R1)) / (1 - cos(α*B1) * exp(-TR*R1)) * exp(-TE*R2)
%
% Inputs:
% - data_struct: struct containing PD, R1, R2, and optionally B1map or Registered_B1
% - TE: echo time (ms)
% - TR: repetition time (ms)
% - alpha: flip angle (degrees)
% - signal_constant: signal scaling constant
% - interp_method: interpolation method for B1 map resizing (default 'linear')

% Handle default argument
if nargin < 6
    interp_method = 'linear';
end

% Input validation
validateattributes(data_struct, {'struct'}, {});
validateattributes(TE, {'numeric'}, {'positive','scalar'});
validateattributes(TR, {'numeric'}, {'positive','scalar'});
validateattributes(alpha, {'numeric'}, {'positive','scalar'});
validateattributes(signal_constant, {'numeric'}, {'positive','scalar'});
validateattributes(interp_method, {'char','string'}, {});

% Extract and validate data
if strcmp(extract_brain, 'no'), PD = double(data_struct.PD.data); else, PD = double(data_struct.PD_brain.data); end
if strcmp(extract_brain, 'no'), R1 = double(data_struct.R1.data); else, R1 = double(data_struct.R1_brain.data); end
if strcmp(extract_brain, 'no'), R2 = double(data_struct.R2.data); else, R2 = double(data_struct.R2_brain.data); end

% Handle B1 correction
if isfield(data_struct, 'B1map_registered') && ~isempty(data_struct.B1map_registered)
    % Use pre-registered B1 map directly, no interpolation
    B1_scaled = double(data_struct.B1map_registered.data);
elseif isfield(data_struct, 'B1map') && ~isempty(data_struct.B1map)
    % Resize and interpolate B1 map to match PD dimensions
    B1_scaled = resize_b1_map(data_struct.B1map.data, PD, interp_method);
else
    % No B1 correction provided, assume 100%%
    B1_scaled = ones(size(PD)) * 100;
end

if isfield(data_struct, 'smap_registered') && ~isempty(data_struct.smap_registered)
    % Use pre-registered B1 map directly, no interpolation
    smap_scaled = double(data_struct.smap_registered.data);
    
elseif isfield(data_struct, 'smap') && ~isempty(data_struct.smap)
    % Resize and interpolate B1 map to match PD dimensions
    smap_scaled = resize_b1_map(data_struct.smap.data, PD, interp_method);
else
    % No B1 correction provided, assume 100%%
    smap_scaled = ones(size(PD)) * 100;
end

% Convert flip angle to radians and apply B1 correction
alpha_rad = deg2rad(alpha);
effective_alpha = alpha_rad .* (B1_scaled / 100);
smap_norm = (smap_scaled - min(smap_scaled(:))) / (max(smap_scaled(:)) - min(smap_scaled(:)));  %normalize to [0,1]




% Convert time units (ms to s) for rate constants
TR_s = TR / 1000;
TE_s = TE / 1000;

% Compute exponential terms
E1 = exp(-TR_s .* R1);
E2 = exp(-TE_s .* R2);

% Compute trigonometric terms
sin_alpha = sin(effective_alpha);
cos_alpha = cos(effective_alpha);

% Ernst equation components
numerator   = sin_alpha .* (1 - E1);
denominator = 1 - cos_alpha .* E1;

% Avoid division by zero
denominator = max(denominator, eps);

% Compute final signal 
%% IMPORTANT: WE HAVE REMOVED MULTIPLICATION WITH SMAP (WHICH WAS AN APPROXIMATION) BY DELETING ".* smap_norm" . to ACTIVATE IT JUST ADD THIS BACK TO THE LINE BELOW. 
%% WE NOW OBTAIN SMAPS USING ESPIRIT NOW FOR EACH COIL SEPERATELY
signal = signal_constant  .* PD .* (numerator ./ denominator) .* E2;

% Handle invalid values
signal(~isfinite(signal)) = 0;
signal(signal < 0) = 0;



end
