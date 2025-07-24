
if isempty(getenv('TOOLBOX_PATH'))
	error('Environment variable TOOLBOX_PATH is not set.');
end

addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));
kspace = rand(100,150);
sensitivities = bart('ecalib', kspace);