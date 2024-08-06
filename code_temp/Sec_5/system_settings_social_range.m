function I = system_settings()
% settings of the interacting particle system on graph, and the integrator
%{
The IPS:   size(a) = N\times N;    size(c) = n\times 1
       \dot{X}^i = \sum_{j=1}^N a_{i,j} K(X^i,X^j) + sigma*dW^i;    i= 1,..., N
            K    = \sum_{k=1}^n c_k K_k
Goal: to learn a and c from data consisting of M trajectories

The system in vector form:
           \dot{X}=  a^\top GammaK c  + sigma*dW
%}

%% setting of the DE system and its integrator
I.viscosity = 0;       % viscosity (noise)
I.N = 6;              % number of agents
I.d = 2;               % dim of state vectors
I.t0       = 0;
I.dt       = 0.001;
I.steps    = 10;
I.obs_std  = 1e-3;     % observation noise

%% Graph
sparsity = 0.5;
I.E = set_graph(I.N, 'sparsity', sparsity, 'plotON', 0);

%% Initial condition
% I.initial = 'DN_1_1_0.5';  % or 'Unif' *********
I.initial = 'N_0_5';  % or 'Unif'
% I.initial = 'Unif_0_5';

%% Choice of kernel
I.basis_case = 6;

%%
I = update_system_settings(I);

%% save directory
I.SAVE_DIR = [getenv('HOME'),'/IPS_Graph_data/'];
if ~exist(I.SAVE_DIR,'dir'); mkdir(I.SAVE_DIR); end
I.SAVE_DIR_fig = [getenv('HOME'),'/IPS_Graph_data/figures/'];
if ~exist(I.SAVE_DIR_fig,'dir'); mkdir(I.SAVE_DIR_fig); end

%% Generate the folder for saving figures for paper
mydir  = pwd;
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end)-1);  % This is the main folder for the project
I.PAPER_FIG_DIR = [newdir, '/Notes/figures'];

end



