function paths_clean = getPathsNoObsNoise( pathObj, m )

%
% function paths_clean = getPathsNoObsNoise( pathObj, m )
%
% Returns the paths without the observation noise
%

% (c) Mauro Maggioni

if nargin<2
    M               = length(pathObj.paths);
    paths_clean     = cell(M,1);
    paths           = pathObj.paths;
    obs_noise       = pathObj.obs_noise;
    parfor m = 1:M
        paths_clean{m}  = paths{m}-obs_noise{m};
    end
else
    paths_clean  = pathObj.paths{m}-pathObj.obs_noise{m};
end


return