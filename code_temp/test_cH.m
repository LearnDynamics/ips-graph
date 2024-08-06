% Test the coercivity
close all;  clear all;   clc;    addpaths;   rng(10);

%% System settings
% Specify the settings here
% In order to fix the figures for paper

I = system_settings(); % setting of the IPS and graph and its integrator
I.viscosity = 0;    % viscosity (noise)
I.N = 6;               % number of agents
I.d = 4;               % dim of state vectors
I.t0       = 0;
I.dt       = 1e-3;     % time steps
I.steps    = 1;      % time steps
I.obs_std  = 0;     % observation noise
I.E = set_graph(I.N, 'sparsity', 0.4, 'plotON', 0);
I.initial = 'N_0_0.1';
I.basis_case = 5;

I = update_system_settings(I);
I.n = length(I.dict);       % Dimension of the hypothesis space
I.phi_kernel = get_kernel_from_c(I.c_true, I.dict);
I.TN = I.steps + 1;                                         % Total time steps
I.tgrid = I.t0:I.dt:(I.steps)*I.dt;
I.X0 = set_particle_initial_all_dim(I.N, I.d, I.initial);   % Initial condition
I.Z_true = get_Z_from_E_c(I.E, I.c_true);     % Z is the product of E and c
I.graph_norm  = norm(I.E,'fro');

%%
B = 500;
M = 10000;
[test_path, I] = get_M_path(I, M, 'ParforProgressON', 1,'saveON', 0, 'loadON', 0);

n = I.n;
d = I.d;
dict = I.dict;
initial = I.initial;
dict_mat = I.dict_mat;

all_cH = zeros(B, 1);
for b = 1:B
    c = randn(n, 1);
    test_kernel = get_kernel_from_c(c, dict);
    %     X0 = set_particle_initial_all_dim(2*M, d, initial);
    %     dif = X0(1:M, :) - X0(M+1:2*M, :);
    %     dis = sqrt(sum(dif.^2, 2));
    %     samples_KX = test_kernel(dis).*dif./dis;
    %     E_KX_2 = sum(samples_KX.^2, 'all')/M;
    %     E_KX = sum(samples_KX, 1)/M;
    %     var_KX = E_KX_2 - sum(E_KX.^2);
    %     norm_kernel_2 = c'*dict_mat*c;
    %     all_cH(b) = var_KX/norm_kernel_2;
    
    all_cH(b) = get_ratio(c, M, d, initial, dict, dict_mat);
    b
end
figure;
histogram(all_cH)


% function r = get_ratio(c, M, d, initial, dict, dict_mat)
% test_kernel = get_kernel_from_c(c, dict);
% X0 = set_particle_initial_all_dim(2*M, d, initial);
% dif = X0(1:M, :) - X0(M+1:2*M, :);
% dis = sqrt(sum(dif.^2, 2));
% samples_KX = test_kernel(dis).*dif./dis;
% E_KX_2 = sum(samples_KX.^2, 'all')/M;
% E_KX = sum(samples_KX, 1)/M;
% var_KX = E_KX_2 - sum(E_KX.^2);
% norm_kernel_2 = c'*dict_mat*c;
% r = var_KX/norm_kernel_2;
% end

