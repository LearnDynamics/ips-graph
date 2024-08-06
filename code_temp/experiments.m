% Learn the adjacency matrix and kernel parameter of an Interacting Particle System on graph.
close all; clear all; addpaths;

gcp;

%% Observations and training data sizes
experimentOption = 4;
experiment_setups

%% Generate training and testing paths
M_test                      = 50;
Z_true                      = get_Z_from_E_c(dyn_sys.A, learning_set.c);     % Z is the product of E and c
fprintf('\nGenerating trajectories, M = %i...',M);tic
trainingPathsObj            = get_paths( dyn_sys, M,        'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
testingPathsObj             = get_paths( dyn_sys, M_test,   'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
fprintf('done (%.2f sec).\n',toc);

% update the exploration measure rho and the basis matrix using paths
learning_set                = update_dict_mat( learning_set, trainingPathsObj.paths );                                          % compute the exploration measure rho and the basis matrix using all paths
learning_set.kernel_norm    = sqrt( learning_set.c'*learning_set.dict_mat*learning_set.c );                                     % MM: this is the estimated kernel norm, assuming the kernel is in the hypothesis space

%%
if false %MM fix this withcorrect initial conditions and forcing terms
    dyn_sys_smalldt             = dyn_sys;
    dyn_sys_smalldt.dt          = dyn_sys.dt/10;
    dyn_sys_smalldt.initial     = testingPathsObj.paths;
    testing_paths_smalldt       = get_paths( dyn_sys_smalldt, M_test, 'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
    fprintf('\n traj err with smaller dt: %f',traj_err(testing_paths,testing_paths_smalldt));                                                                                  % MM: sanity test: compare with simulations with smaller dt
    dyn_sys_samedt             = dyn_sys;
    dyn_sys_samedt.dt          = dyn_sys.dt;
    dyn_sys_samedt.initial     = testingPathsObj.paths;
    testing_paths_samedt       = get_paths( dyn_sys_samedt, M_test, 'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
    fprintf('\n traj err with same dt: %f',traj_err(testing_paths,testing_paths_samedt));                                                                                  % MM: sanity test: compare with simulations with smaller dt
    dyn_sys_largedt             = dyn_sys;
    dyn_sys_largedt.dt          = dyn_sys.dt*10;
    dyn_sys_largedt.initial     = testingPathsObj.paths;
    testing_paths_largedt       = get_paths( dyn_sys_largedt, M_test, 'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
    fprintf('\n traj err with larger dt: %f',traj_err(testing_paths,testing_paths_largedt));                                                                                  % MM: sanity test: compare with simulations with smaller dt
end

%% Learning using ORALS
[estORSVD.A, estORSVD.c, estORALS.A, estORALS.c, estORALS.Z, ~, condA_orals,~,estORALS.time] = ...
    learn_kernel_graph_ORALS(trainingPathsObj.paths, dyn_sys, learning_set, 'plotON', 0, 'reg_method', 'pinv');                 % reg_methods: ID, RKHS, None, lsqminnorm

%% Learning using alternating least squares
fprintf('\nALS...'); tic;
estALS     = learn_kernel_graph_ALS( trainingPathsObj.paths, dyn_sys, learning_set, 'niter', 10, 'normalizeON', 1, ...
    'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
    'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
    'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );
fprintf('done (in %f sec., %d iter).',toc,size(estALS.stats.c_hat,2));

%% Compute some error statistics
estORSVD.A_err  = graph_err ( estORSVD.A, dyn_sys );
estORSVD.c_err  = kernel_err( estORSVD.c, learning_set );
estORSVD.Z_err  = get_Z_error( get_Z_from_E_c( estORSVD.A, estORSVD.c ), dyn_sys, learning_set );
estORALS.A_err  = graph_err ( estORALS.A, dyn_sys );
estORALS.c_err  = kernel_err( estORALS.c, learning_set );
estORALS.Z_err  = get_Z_error( get_Z_from_E_c( estORALS.A, estORALS.c ), dyn_sys, learning_set );
estALS.A_err    = graph_err ( estALS.A, dyn_sys );
estALS.c_err    = kernel_err( estALS.c, learning_set );
estALS.Z_err    = get_Z_error( get_Z_from_E_c( estALS.A, estALS.c ), dyn_sys, learning_set );

try
    estALS.ALS_error_kernel    = zeros(size(estALS.stats.c_hat,2),1);
    estALS.ALS_error_graph     = zeros(size(estALS.stats.c_hat,2),1);

    for q = 1:size(estALS.stats.c_hat,2)
        estALS.ALS_error_kernel(q) = kernel_err( estALS.stats.c_hat(:,q), learning_set );
        estALS.ALS_error_graph(q)  = graph_err ( estALS.stats.A_hat(:,:,q), dyn_sys );
    end
catch
end

[estORSVD.pathTestErr,~,estORSVD.meanL2traj]    = getPathTestError( dyn_sys, learning_set, estORSVD.A, estORSVD.c, testingPathsObj );
[estORALS.pathTestErr,~,estORALS.meanL2traj]    = getPathTestError( dyn_sys, learning_set, estORALS.A, estORALS.c, testingPathsObj );
[estALS.pathTestErr,~,estALS.meanL2traj]        = getPathTestError( dyn_sys, learning_set, estALS.A, estALS.c, testingPathsObj );

%% Report estimation errors
ErrorType   = ["Kernel";"Graph";"Tensor Z";"Paths";"Time"];
ErrorALS    = [ estALS.c_err; estALS.A_err; estALS.Z_err; estALS.pathTestErr;estALS.time.ALS];
ErrorORALS  = [ estORALS.c_err; estORALS.A_err; estORALS.Z_err; estORALS.pathTestErr; estORALS.time];
ErrorORSVD  = [ estORSVD.c_err; estORSVD.A_err; estORSVD.Z_err; estORSVD.pathTestErr; estORALS.time];
ErrorTable  = table(ErrorType,ErrorALS,ErrorORALS,ErrorORSVD);

ErrorTable

%% Plots
main_plots