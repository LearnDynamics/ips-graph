function pathObj = get_paths( dyn_sys, M, varargin )

% OUT: pathObj is a structure with the following fields:
%   paths           : cell array of length M, with m-th cell representing the m-th path with an array of size N x dim x L
%   forcing_noise   : cell array of length M, with m-th cell representing the forcing noise on the m-th path with an array of size N x dim x L
%   obs_noise       : cell array of length M, with m-th cell representing the observation noise on the m-th path with an array of size N x dim x L

%% Input parser
p = inputParser;
addRequired(p, 'dyn_sys');
addRequired(p, 'M');
addOptional(p, 'saveON', 0);
addOptional(p, 'ParforProgressON', 0);
addOptional(p, 'loadON', 0);
addOptional(p, 'forcing_noise', []);

parse(p, dyn_sys, M, varargin{:});
saveON                  = p.Results.saveON;
loadON                  = p.Results.loadON;
forcing_noise           = p.Results.forcing_noise;
forcing_noise_provided  = ~isempty(forcing_noise);
ParforProgressON        = p.Results.ParforProgressON;

%% Check if we have local path file
if loadON
    Datadir = get_Data_dir( M, dyn_sys );
    if exist(Datadir, 'file')
        fprintf('Loading local trajectory file...')
        load(Datadir, 'paths','dyn_sys');
        return
    end
    
end

%% Generate paths
progressON      = 0;
rk4ON           = false;
paths           = cell(M, 1);
if ~forcing_noise_provided,    forcing_noise   = cell(M, 1); end
obs_noise       = cell(M, 1);

%if ParforProgressON;WaitMessage = parfor_wait(M, 'Waitbar', true);                                                             % MM: does not work with ThreadPools
%else; WaitMessage = 0;end

for i = 1:M
    if ischar( dyn_sys.initial )      %% MM:TBD: document this behavior for initial conditions
        X0                      = set_particle_initial_all_dim( dyn_sys.N, dyn_sys.d, dyn_sys.initial );                        %#ok<*PFBNS> 
    elseif iscell( dyn_sys.initial )
        X0                      = set_particle_initial_all_dim( dyn_sys.N, dyn_sys.d, squeeze(dyn_sys.initial{i}(:,:,1) ));     %#ok<*PFBNS> 
    else
        X0                      = set_particle_initial_all_dim( dyn_sys.N, dyn_sys.d, squeeze(dyn_sys.initial(i,:,:) )');        %#ok<*PFBNS> 
    end
    if ~forcing_noise_provided
        [xpath,forcing_noise{i}]= graph_forward_model( dyn_sys, X0, progressON, rk4ON);                                         % generate paths, save stoch. forcing term
    else
        xpath                   = graph_forward_model( dyn_sys, X0, progressON, rk4ON, forcing_noise{i});                       % generate paths, with the stoch. forcing term provided
    end
    if dyn_sys.obs_std>0
        obs_noise{i}            = dyn_sys.obs_std * randn(size(xpath));
        paths{i}                = xpath + obs_noise{i};
    else
        obs_noise{i}            = 0;
        paths{i}                = xpath;
    end
    
%    if ParforProgressON;WaitMessage.Send;end    
end

%if ParforProgressON;WaitMessage.Destroy;end

%% save data
if saveON        
    Datadir = get_Data_dir(M, dyn_sys, learning_setup);    
    if ~exist(dyn_sys.SAVE_DIR,'dir')
        mkdir(dyn_sys.SAVE_DIR); 
    end
    try
    save(Datadir, 'paths','dyn_sys');
    catch
        %MM save is buggy, if it fails let's not waste the computations done so far
    end
end

pathObj.paths           = paths;
pathObj.forcing_noise   = forcing_noise;
pathObj.obs_noise       = obs_noise;

end


function Datadir = get_Data_dir(M, dyn_sys, learning_setup)

% setting_name = ['M_', num2str(M), '_N_', num2str(I.N), '_steps_', num2str(I.L),'_ic',I.initial];

str_name    = sprintf('M_%i_N_%i_steps_%i_basis_case%i_d%i_obsnr', M, dyn_sys.N, dyn_sys.L, learning_setup.basis_case, dyn_sys.d);
str_name    = [str_name, num2str(dyn_sys.obs_std),'_ic',dyn_sys.initial, '_vis_', num2str(dyn_sys.viscosity)];

Datadir    = [dyn_sys.SAVE_DIR, 'all_path_', str_name, '.mat'];

end