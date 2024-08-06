% Section 3.1 Typical example
% Use code based on B tensors

% Using Lenard Jones kernel, with rendudant basis
% See details below for the setting of kernel_type
% and the construction of the basis in learning_settings

% Todo:
% Finish the examples
% Show some interesting trajectories


%% path and clear

close all; clear all; addpaths;clc
rng(2);
% gcp;        % get current parallel pool



%% specify the graph and the desired initial conditions in the figures
options = 'spiral';
options = 'two_clusters';
% options = 'load';
switch options

    case 'spiral'
        % Set a initial condition
        X0                      = [0.5241    0.3062;
            0.3560    0.6185;
            0.9517    0.0557;
            0.8728    0.2870;
            0.6836    0.9951;
            0.6742    0.6197];

        % set a graph
        A = zeros(6, 6);
        A(3, 1) = 1;
        A(4, 3) = 1;
        A(6, 4) = 1;
        A(5, 6) = 1;
        A(2, 5) = 1;
        A(1, 2) = 1;

        % Observations and training data sizes
        M                       = 1000;                                                 % number of trajectories
        dyn_sys.L               = 50;      % number observations equispaced in time [0,T]

        % Settings for dynamics and learning
        dyn_sys.N               = 6;                                                                                                   % number of agents
        dyn_sys.d               = 2;                                                                                                    % dim of state vectors
        dyn_sys.viscosity       = 1e-3;  %  set 0 to test consistency                                       % stochastic force; viscosity (forcing noise in the dynamics)
        % dyn_sys.A               = set_graph(dyn_sys.N, 'sparsity', 0.4, 'plotON', 0);                                                   % create influence graph
        dyn_sys.A = A;
        dyn_sys.initial         = 'Unif_0_1';
        dyn_sys.dt              = 1e-4;
        dyn_sys.T               = dyn_sys.dt*(dyn_sys.L-1);
        dyn_sys.obs_std         = 1e-3;  %  set 0 to test consistency                                       % observation noise
        dyn_sys                 = system_settings( dyn_sys );                                                                           % settings of the IPS and graph and its integrator






    case 'two_clusters'
        X0 = [1.3001    0.1679;
    0.7814    0.4799;
    0.7454    0.5858;
    0.6136    1.1098;
    0.0988    1.2428;
    1.2959    0.6399];

        A =          [0         0         0         0    0.3635    0.7104;
                0.7144         0         0    0.9702         0         0;
                     0         0         0         0         0         0;
                0.6997    0.8523    0.9933         0    0.9316    0.7038;
                     0         0    0.1154    0.2425         0         0;
                     0    0.5231         0         0         0         0];
        % Observations and training data sizes
        M                       = 1000;                                                 % number of trajectories
        dyn_sys.L               = 50;      % number observations equispaced in time [0,T]

        % Settings for dynamics and learning
        dyn_sys.N               = 6;                                                                                                   % number of agents
        dyn_sys.d               = 2;                                                                                                    % dim of state vectors
        dyn_sys.viscosity       = 1e-3;  %  set 0 to test consistency                                       % stochastic force; viscosity (forcing noise in the dynamics)
        % dyn_sys.A               = set_graph(dyn_sys.N, 'sparsity', 0.4, 'plotON', 0);                                                   % create influence graph
        dyn_sys.A = A;
        dyn_sys.initial         = 'Unif_0_1.5';
        dyn_sys.dt              = 1e-4;
        dyn_sys.T               = dyn_sys.dt*(dyn_sys.L-1);
        dyn_sys.obs_std         = 1e-3;  %  set 0 to test consistency                                       % observation noise
        dyn_sys                 = system_settings( dyn_sys );
end







kernel_type             = 'typical_example_Lenard_Jones'; % 'randomSmoothFourierWithDecay'; % 6; %                  % see options in learning_settings.m

learning_setup            = learning_settings( kernel_type, dyn_sys);
learning_setup.Z_true     = get_Z_from_E_c( dyn_sys.A, learning_setup.c );       % Z is the product of A and c                                                      % Z is the product of E and c

dyn_sys.phi_kernel      = learning_setup.phi_kernel_cheb;                           % MM: logical inconsistency here and in the following because of this; learning_set.phi_kernel not used in any of what follows
dyn_sys.phi_kernel      = learning_setup.phi_kernel;                                % FL: learning_set.phi_kernel serves as the ground truth. --- It is a parametric inference. The chev is towards non-parametric inference
dyn_sys.n               = learning_setup.n;


%% a place to check the traj before we do the learning
true_para.A             = dyn_sys.A;
true_para.phi_kernel    = dyn_sys.phi_kernel;
% dyn_sys_temp = dyn_sys;
% dyn_sys_temp.dt = 1e-4;
% dyn_sys_temp.L = 100;
% 
% dyn_sys_temp.obs_std = 0;
% dyn_sys_temp.viscosity = 0;
% 
% dyn_sys_temp.tgrid = (1:dyn_sys_temp.L)*dyn_sys_temp.dt;
% true_path_temp               = generate_pred_path(true_para,  dyn_sys_temp, X0);
% %
% figure; subplot(121); hold on; view(3);grid on
% for i = 1:dyn_sys.N
%     plot3(squeeze(true_path_temp(i, 1, :)), squeeze(true_path_temp(i, 2, :)), dyn_sys_temp.tgrid, 'linewidth', 3)
% end
% 
% subplot(122)
% imagesc(dyn_sys_temp.A')
%% Generate training and testing paths
% M_test                      = 50;
fprintf('\nGenerating trajectories, M = %i...',M);tic
trainingPathsObj            = get_paths( dyn_sys, M,        'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
fprintf('done (%.2f sec).',toc);



% update the exploration measure rho and the basis matrix using paths
learning_setup                = update_dict_mat( learning_setup, trainingPathsObj.paths);
% compute the exploration measure rho and the basis matrix using all paths
% if learning_setup.basis_changed ==1
%     learning_setup.Z_true_orig    = learning_setup.Z_true;        % just in case a change of basis is applied when dict_mat is singular
%     learning_setup.Z_true         = get_Z_from_E_c( dyn_sys.A, learning_setup.c );
% end
learning_setup.kernel_norm    = sqrt( learning_setup.c'*learning_setup.dict_mat*learning_setup.c );                                     % MM: this is the estimated kernel norm, assuming the kernel is in the hypothesis space



%% ALS and ORALS

estALS = learn_kernel_graph_ALS( trainingPathsObj.paths, dyn_sys, learning_setup, 'niter', 10, 'normalizeON', 1, ...
    'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
    'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
    'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );

matFactorize  = 'ALS';    % %  'SVD', 'ALS', 'SVD+ORALS'
estORALS = learn_kernel_graph_ORALS_B(trainingPathsObj.paths, dyn_sys, matFactorize, ...
    learning_setup, 'plotON', 0, 'reg_method', 'lsqminnorm');



%% Generate a sample path
pred_path_ALS       = generate_pred_path(estALS,    dyn_sys, X0);
pred_path_ORALS     = generate_pred_path(estORALS,  dyn_sys, X0);
true_path           = generate_pred_path(true_para,  dyn_sys, X0);



%% Ensemble trajectory error
N = 100;
pred_path_ALS_all = cell(N, 1);
pred_path_ORALS_all = cell(N, 1);
true_path_all = cell(N, 1);

for i = 1:N
    % i
    pred_path_ALS_all{i}       = generate_pred_path(estALS,    dyn_sys, X0);
    pred_path_ORALS_all{i}     = generate_pred_path(estORALS,  dyn_sys, X0);
    true_path_all{i}           = generate_pred_path(true_para,  dyn_sys, X0);
end



%%
ErrorALS(1) = kernel_err(estALS.c, learning_setup);
ErrorALS(2) = graph_err(estALS.A, dyn_sys);
ErrorALS(3) = traj_err( {true_path}, {pred_path_ALS} );
[ALS_ensemble_traj_error_mean, ~, ~, ALS_ensemble_traj_error_all] = traj_err(true_path_all, pred_path_ALS_all);
ErrorALS(4) = ALS_ensemble_traj_error_mean;
ErrorALS(5) = std(ALS_ensemble_traj_error_all);

ErrorORALS(1) = kernel_err(estORALS.c, learning_setup);
ErrorORALS(2) = graph_err(estORALS.A, dyn_sys);
ErrorORALS(3) = traj_err( {true_path}, {pred_path_ORALS} );

[ORALS_ensemble_traj_error_mean, ~, ~, ORALS_ensemble_traj_error_all] = traj_err(true_path_all, pred_path_ORALS_all);
ErrorORALS(4) = ORALS_ensemble_traj_error_mean;
ErrorORALS(5) = std(ORALS_ensemble_traj_error_all);

% sqrt(sum((true_path - pred_path_ORALS).^2*dyn_sys.dt, 'all'));

ErrorALS =ErrorALS';
ErrorORALS =ErrorORALS';
clear ErrorTable; 
ErrorType   = ["Kernel";"Graph";"Path in fig";"Path Ensemble mean";"Path Ensemble SD"];
ErrorTable  = table(ErrorType);
ErrorTable = addvars(ErrorTable,ErrorALS); 
ErrorTable = addvars(ErrorTable,ErrorORALS);  

fprintf('\n');
disp(ErrorTable)






%%
plot_typical_example




%%
% ErrorALS(1) = kernel_err(estALS.c, learning_setup);
% ErrorALS(2) = graph_err(estALS.A, dyn_sys);
% 
% 
% ErrorORALS(1) = kernel_err(estORALS.c, learning_setup);
% ErrorORALS(2) = graph_err(estORALS.A, dyn_sys);
% 
% 
% 
% % sqrt(sum((true_path - pred_path_ORALS).^2*dyn_sys.dt, 'all'));
% 
% ErrorALS =ErrorALS';
% ErrorORALS =ErrorORALS';
% 
% clear ErrorTable; 
% ErrorType   = ["Kernel";"Graph";"Ensemble Paths"];
% ErrorTable  = table(ErrorType);
% ErrorTable = addvars(ErrorTable,ErrorALS); 
% ErrorTable = addvars(ErrorTable,ErrorORALS);  
% 
% fprintf('\n');
% disp(ErrorTable)


