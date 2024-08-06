% Section 4.3 Multitype kernel
% Use code based on B tensors
% Use 3 way ALS to learn multitype kernel

% Quanjun Lang (c) 2024
close all; clear all; addpaths;clc

rng(20);

% gcp;        % get current parallel pool

%% Observations and training data sizes
M                       = 400;                                                 % number of trajectories
dyn_sys.L               = 50;      % number observations equispaced in time [0,T]
runALS                  = true;
runORALS                = true;   % using B tensor and long matrix
runORALS_normalMat      = true;   % without using B tensor, use normal matrices 

% Settings for dynamics and learning

dyn_sys.N               = 10;                                                                                                   % number of agents
dyn_sys.d               = 2;                                                                                                    % dim of state vectors
dyn_sys.viscosity       = 1e-3;  %  set 0 to test consistency                                       % stochastic force; viscosity (forcing noise in the dynamics)
dyn_sys.A               = set_graph(dyn_sys.N, 'sparsity', 0.4, 'plotON', 0);                                                   % create influence graph
dyn_sys.initial         = 'Unif_0_5';
dyn_sys.dt              = 1e-3;
dyn_sys.T               = dyn_sys.dt*(dyn_sys.L-1);
dyn_sys.obs_std         = 1e-3;  %  set 0 to test consistency                                       % observation noise
dyn_sys                 = system_settings( dyn_sys );                                                                           % settings of the IPS and graph and its integrator



M_test = 10;
dt_test= 1e-2;
L_test = 500;

%% Generate data with 2 kernels and learning using ALS-3 and ALS   %%%%%%%%

% Generate data
kernel_type                 = 'multitype';  %'multi_type';
basis_choice                = 'spline_long_short';
learning_set                = learning_settings( kernel_type, dyn_sys, basis_choice);    
dyn_sys.n                   = learning_set.n; 
dyn_sys.phi_kernel_choices  = learning_set.phi_kernel_choices;
dyn_sys.kernel_idx          = learning_set.kernel_idx;

fprintf('\nGenerating trajectories, M = %i...',M);tic
trainingPathsObj            = get_paths_multitype( dyn_sys, M,        'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
fprintf('done (%.2f sec).',toc);
learning_set                = update_dict_mat( learning_set, trainingPathsObj.paths, 'rho', 0);                                          % compute the exploration measure rho and the basis matrix using all paths

% Learning using Three ways ALS
Q = learning_set.num_kernel_choices;

estALS_3_multi     = learn_kernel_graph_ALS_multitype( trainingPathsObj.paths, dyn_sys, learning_set, 'niter', 10, 'normalizeON', 1, ...
    'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_u', 1e-3, ...
    'stop_thres_relDelta_v', 1e-3, 'plotON', 1, 'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
    'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] ,'K_means_on_v', 1);% 1e-2*norm(I.A,'fro') );

learning_set.sorted_id          = sort_id(learning_set.kernel_idx)';
estALS_3_multi.sorted_id              = sort_id(estALS_3_multi.kernel_idx)';
estALS_3_multi.sorted_multi_kernel    = get_multi_type_kernel_from_id(estALS_3_multi.coef_mat, estALS_3_multi.sorted_id, Q, learning_set.dict);

% Learning using Two ways ALS
estALS_2_multi = learn_kernel_graph_ALS( trainingPathsObj.paths, dyn_sys, learning_set, 'niter', 10, 'normalizeON', 1, ...
    'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
    'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
    'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );

% Generating testing trajectory
dyn_sys_test            = dyn_sys;
dyn_sys_test.L          = L_test;
dyn_sys_test.dt         = dt_test;
dyn_sys_test.viscosity  = 0;
dyn_sys_test.obs_std    = 0;
dyn_sys_test.T          = dyn_sys_test.dt*(dyn_sys_test.L-1);
testingPathsObj_multi         = get_paths_multitype( dyn_sys_test, M_test,   'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
initial_val         = testingPathsObj_multi.initial;

% Generating prediction trajectory
estALS_3_multi.pred_path  = generate_pred_path_multitype(estALS_3_multi, dyn_sys_test, initial_val);
estALS_2_multi.pred_path    = generate_pred_path(estALS_2_multi, dyn_sys_test, initial_val);

% Trajectory prediction error
[estALS_3_multi.traj_err_rel, ~, estALS_3_multi.traj_err, estALS_3_multi.traj_err_rel_seq]   = traj_err( testingPathsObj_multi.paths, estALS_3_multi.pred_path );
[estALS_2_multi.traj_err_rel, ~, estALS_2_multi.traj_err, estALS_2_multi.traj_err_rel_seq]     = traj_err( testingPathsObj_multi.paths, estALS_2_multi.pred_path );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data with 1 kernel and learning using ALS-3 and ALS   %%%%%%%%
kernel_type             = 'multitype';  %'multi_type';
basis_choice            = 'fake_multitype';
learning_set            = learning_settings( kernel_type, dyn_sys, basis_choice);    

dyn_sys.n                   = learning_set.n; 
dyn_sys.phi_kernel_choices  = learning_set.phi_kernel_choices;
dyn_sys.kernel_idx          = learning_set.kernel_idx;

fprintf('\nGenerating trajectories, M = %i...',M);tic
trainingPathsObj            = get_paths_multitype( dyn_sys, M,        'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );

fprintf('done (%.2f sec).',toc);

learning_set                    = update_dict_mat( learning_set, trainingPathsObj.paths, 'rho', 0);                                          % compute the exploration measure rho and the basis matrix using all paths


% Learning using Three ways alternating least squares
Q = learning_set.num_kernel_choices;

estALS_3_single     = learn_kernel_graph_ALS_multitype( trainingPathsObj.paths, dyn_sys, learning_set, 'niter', 10, 'normalizeON', 1, ...
    'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_u', 1e-3, ...
    'stop_thres_relDelta_v', 1e-3, 'plotON', 1, 'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
    'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] ,'K_means_on_v', 1);% 1e-2*norm(I.A,'fro') );

learning_set.sorted_id          = sort_id(learning_set.kernel_idx)';
estALS_3_single.sorted_id              = sort_id(estALS_3_single.kernel_idx)';
estALS_3_single.sorted_multi_kernel    = get_multi_type_kernel_from_id(estALS_3_single.coef_mat, estALS_3_single.sorted_id, Q, learning_set.dict);

estALS_2_single = learn_kernel_graph_ALS( trainingPathsObj.paths, dyn_sys, learning_set, 'niter', 10, 'normalizeON', 1, ...
    'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
    'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
    'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );


% Generating testing trajectory
dyn_sys_test            = dyn_sys;
dyn_sys_test.L          = L_test;
dyn_sys_test.dt         = dt_test;
dyn_sys_test.viscosity  = 0;
dyn_sys_test.obs_std    = 0;
dyn_sys_test.T          = dyn_sys_test.dt*(dyn_sys_test.L-1);
testingPathsObj_single         = get_paths_multitype( dyn_sys_test, M_test,   'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );

initial_val         = testingPathsObj_single.initial;
estALS_3_single.pred_path  = generate_pred_path_multitype(estALS_3_single, dyn_sys_test, initial_val);
estALS_2_single.pred_path    = generate_pred_path(estALS_2_single, dyn_sys_test, initial_val);

[estALS_3_single.traj_err_rel, ~, estALS_3_single.traj_err, estALS_3_single.traj_err_rel_seq]   = traj_err( testingPathsObj_single.paths, estALS_3_single.pred_path );
[estALS_2_single.traj_err_rel, ~, estALS_2_single.traj_err, estALS_2_single.traj_err_rel_seq]   = traj_err( testingPathsObj_single.paths, estALS_2_single.pred_path );

%% Print out the error table

clear ErrorTable; 

Q_true_1 = zeros(4, 1);
Q_true_2 = zeros(4, 1);


Q_true_1(1) = estALS_2_single.traj_err_rel;
Q_true_1(2) = std(estALS_2_single.traj_err_rel_seq);
Q_true_1(3) = estALS_3_single.traj_err_rel;
Q_true_1(4) = std(estALS_3_single.traj_err_rel_seq);


Q_true_2(1) = estALS_2_multi.traj_err_rel;
Q_true_2(2) = std(estALS_2_multi.traj_err_rel_seq);
Q_true_2(3) = estALS_3_multi.traj_err_rel;
Q_true_2(4) = std(estALS_3_multi.traj_err_rel_seq);


ErrorType   = ["Q est 1 , traj erro mean"; "Q est 1 , traj error std"; "Q est 2, traj erro mean"; "Q est 2, traj error std"; ];
ErrorTable  = table(ErrorType);
ErrorTable = addvars(ErrorTable,Q_true_1); 
ErrorTable = addvars(ErrorTable,Q_true_2);  

fprintf('\n');
disp(ErrorTable)


% estALS_2_multi.traj_err
% estALS_3_multi.traj_err
% 
% estALS_2_multi.traj_err_rel
% estALS_3_multi.traj_err_rel
% 
% 
% estALS_2_single.traj_err
% estALS_3_single.traj_err
% 
% estALS_2_single.traj_err_rel
% estALS_3_single.traj_err_rel


%% plot some trajectoryies


figure;hold on;
ind     = 5;
axis    = 1;
plot(squeeze(estALS_2_single.pred_path{ind}(:, axis, :))'    , '.'      , 'color', 'red', 'LineWidth', 2)
plot(squeeze(estALS_3_single.pred_path{ind}(:, axis, :))' , '.'       , 'color', 'black', 'LineWidth', 2)
plot(squeeze(testingPathsObj_single.paths{ind}(:, axis, :))', 'color', 'blue', 'LineWidth', 4)

figure;hold on;
ind     = 5;
axis    = 1;
plot(squeeze(estALS_2_multi.pred_path{ind}(:, axis, :))'    , '.'      , 'color', 'red', 'LineWidth', 2)
plot(squeeze(estALS_3_multi.pred_path{ind}(:, axis, :))' , '.'       , 'color', 'black', 'LineWidth', 2)
plot(squeeze(testingPathsObj_multi.paths{ind}(:, axis, :))', 'color', 'blue', 'LineWidth', 4)





%% PLot for 3 way ALS on single kernel trajectory

% Section_4_3_helper_ALS3_on_single_kernel_traj

%% Plot for the general result, comparing Kmeans and not using Kmeans
% Section_4_3_helper_compare_kmeans