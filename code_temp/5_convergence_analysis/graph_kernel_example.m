% Example of graph and kernel 

close all;  clear all;   clc;    addpaths;
rng(1)
%% System settings
% Specify the settings here
% In order to fix the figures for paper

I = system_settings(); % setting of the IPS and graph and its integrator
I.viscosity = 1e-4;    % viscosity (noise)
I.N         = 10;               % number of agents
I.d         = 1;               % dim of state vectors
I.t0        = 0;
I.dt        = 1e-3;     % time steps
I.steps     = 100;      % time steps
I.obs_std   = 1e-4;     % observation noise
I.E = set_graph(I.N, 'sparsity', 0.4, 'plotON', 0);
I.initial = 'Unif_0_5';
I.basis_case = 6; 

I = update_system_settings(I);