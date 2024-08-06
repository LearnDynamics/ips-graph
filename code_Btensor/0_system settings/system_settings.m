function dyn_sys = system_settings( dyn_sys )

% settings of the interacting particle system on graph, and the integrator
%{
The IPS:   size(a) = N\times N;    size(c) = n\times 1
       \dot{X}^i = \sum_{j=1}^N a_{i,j} K(X^i,X^j) + sigma*dW^i;    i= 1,..., N
            K    = \sum_{k=1}^n c_k K_k
Goal: to learn a and c from data consisting of M trajectories

The system in vector form:
           \dot{X}=  a^\top GammaK c  + sigma*dW
%}

%% Settings of the Dynamical System and its integrator
if ~isfield(dyn_sys,'N')
    dyn_sys.N            = 6;                                                                                                   % number of agents
end
if ~isfield(dyn_sys,'d')
    dyn_sys.d            = 2;                                                                                                   % dim of state vector at each node
end
if ~isfield(dyn_sys,'dt')
    dyn_sys.dt           = 0.001;                                                                                               % step size of simulator
end
if ~isfield(dyn_sys,'T')
    dyn_sys.T            = 1;
end
if ~isfield(dyn_sys,'viscosity')
    dyn_sys.viscosity    = 0;                                                                                                   % viscosity (forcing term noise)
end

% settings for the observations
if ~isfield(dyn_sys,'obs_std')
    dyn_sys.obs_std      = 1e-3;                                                                                                % observation noise
end
if ~isfield(dyn_sys,'L')
    dyn_sys.L            = 10;                                                                                                  % number of observations in [t_0,T]
end
if ~isfield(dyn_sys,'tgrid')
    dyn_sys.tgrid        = linspace(0,dyn_sys.T,dyn_sys.L);
end

%% Graph
if ~isfield(dyn_sys,'sparsity')
    dyn_sys.sparsity     = 0.5;
end
if ~isfield(dyn_sys,'A')
    dyn_sys.A            = set_graph(dyn_sys.N, 'sparsity', dyn_sys.sparsity, 'plotON', 0);                            % generate influence graph
end

%% Initial condition
if ~isfield(dyn_sys,'initial')
    dyn_sys.initial      = 'N_0_5';  % or 'Unif_0_5', or 'DN_1_1_0.5'
end
if ~isfield(dyn_sys,'X0')
    dyn_sys.X0           = set_particle_initial_all_dim(dyn_sys.N, dyn_sys.d, dyn_sys.initial);   % Initial condition
end

%% Save directory
if ~isfield(dyn_sys,'SAVE_DIR')
    dyn_sys.SAVE_DIR = [getenv('HOME'),'/IPS_Graph_data/'];
    if ~exist(dyn_sys.SAVE_DIR,'dir'); mkdir(dyn_sys.SAVE_DIR); end
    dyn_sys.SAVE_DIR_fig = [getenv('HOME'),'/IPS_Graph_data/figures/'];
    if ~exist(dyn_sys.SAVE_DIR_fig,'dir'); mkdir(dyn_sys.SAVE_DIR_fig); end
end

%% Generate the folder for saving figures for paper
if ~isfield(dyn_sys,'PAPER_FIG_DIR')
    mydir  = pwd;
    idcs   = strfind(mydir,'/');
    newdir = mydir(1:idcs(end)-1);  % This is the main folder for the project
    dyn_sys.PAPER_FIG_DIR = [newdir, '/figures'];
end


% %% Generate the folder for saving figures for slides
% if ~isfield(dyn_sys,'SLIDE_FIG_DIR')
%     dyn_sys.SLIDE_FIG_DIR = [fileparts(PAPER_FIG_DIR), '/Slides/myfigures'];
% end
end