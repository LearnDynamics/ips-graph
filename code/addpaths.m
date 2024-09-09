%% Add code path to the environment
folder = fileparts(which(mfilename));
addpath(genpath(folder));
clear folder

if ~exist('chebfun')
    fprintf('Please download chebfun for improved performance\n');
end
