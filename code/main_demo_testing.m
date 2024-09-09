% Learn the weight/adjacency matrix and kernel parameter of an Interacting Particle System on graph.
% @copyright: Quanjun Lang, Mauro Maggioni, Fei Lu, Xiong Wang

close all; clear all; addpaths;

%% Observations and training data sizes
M                       = 128;                                                                                                   % number of trajectories
dyn_sys.L               = 16;                                                                                                    % number observations equispaced in time [0,T]
runALS                  = true;
runORALS                = true;                                                                                                 % using B tensor and long matrix                                                                                               % without using B tensor, use normal matrices 

%% Settings for dynamics and learning

dyn_sys.N               = 6;                                                                                                    % number of agents
dyn_sys.d               = 2;                                                                                                    % dim of state vectors
dyn_sys.viscosity       = 1e-3;  %  set 0 to test consistency                                                                   % stochastic force; viscosity (forcing noise in the dynamics)
dyn_sys.A               = set_graph(dyn_sys.N, 'sparsity', 0.4, 'plotON', 0);                                                   % create influence graph
dyn_sys.initial         = 'Unif_0_3';
dyn_sys.dt              = 1e-3;
dyn_sys.T               = dyn_sys.dt*(dyn_sys.L-1);
dyn_sys.obs_std         = 1e-2;  %  set 0 to test consistency                                                                   % observation noise
dyn_sys                 = system_settings( dyn_sys );                                                                           % settings of the IPS and graph and its integrator


kernel_type             =  6; % 'randomSmoothFourierWithDecay'; % 6; %                                                          % see options in learning_settings.m
n =  16;                                                                                                                        %  it is used in learning_settings

learning_setup          = learning_settings( kernel_type, dyn_sys, struct('n',n) );
learning_setup.Z_true   = get_Z_from_E_c( dyn_sys.A, learning_setup.c );                                                        % Z is the product of A and c

dyn_sys.phi_kernel      = learning_setup.phi_kernel;                           
dyn_sys.n               = learning_setup.n; 

%% Demo (1 traj) generate particles
% Sanity check.
% If the kernel is bad, it will generate divergence particle trajectory.

progressON  = true;
plotON      = false;
RIPcompON   = false;

if  ~plotON
    xpath = graph_forward_model( dyn_sys, dyn_sys.X0, progressON );
else % Plot trajectories and the graph
    I1 = dyn_sys;
    I1.steps = 1e4;
    xpath = graph_forward_model(I1, I1.X0, progressON);

    if plotON &&  I1.d == 2
        graph_plot_motion(xpath, I1, plotON);
        plot_graph(I1.A);
        figure;fplot(I1.phi_kernel, [0, 5]); title('True interaction kernel')
    end

    if sum(abs(xpath(:, :, end)), 'all') > 1e10
        disp('trajectory is divergent. Kernel is bad'); return;
    end
end

%% Generate training and testing paths
M_test                      = 50;
fprintf('\nGenerating trajectories, M = %i...',M);tic
trainingPathsObj            = get_paths( dyn_sys, M,        'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
testingPathsObj             = get_paths( dyn_sys, M_test,   'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
fprintf('done (%.2f sec).',toc);

% update the exploration measure rho and the basis matrix using paths
learning_setup                  = update_dict_mat( learning_setup, trainingPathsObj.paths);                                       % compute the exploration measure rho and the basis matrix using all paths
if learning_setup.basis_changed ==1
    learning_setup.Z_true_orig    = learning_setup.Z_true;                                                                      % just in case a change of basis is applied when dict_mat is singular
    learning_setup.Z_true         = get_Z_from_E_c( dyn_sys.A, learning_setup.c );
end

learning_setup.kernel_norm    = sqrt( learning_setup.c'*learning_setup.dict_mat*learning_setup.c );                             % this is the estimated kernel norm, assuming the kernel is in the hypothesis space

%% Learning using ORALS
matFactorize  = 'ALS';    % %  'SVD', 'ALS', 'SVD+ORALS'
if runORALS 
   fprintf('\nRunning ORALS...');
   estORALS = learn_kernel_graph_ORALS_B(trainingPathsObj.paths, dyn_sys, matFactorize, learning_setup,...
          'plotON', 0, 'reg_method', 'lsqminnorm');         % reg_methods: ID, RKHS, None, lsqminnorm,pinv, pinvreg    
end

%% Learning using Alternating Least Squares
if runALS
    fprintf('\nRunning ALS...');
    estALS     = learn_kernel_graph_ALS( trainingPathsObj.paths, dyn_sys, learning_setup, 'niter', 10, 'normalizeON', 1, ...
        'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
        'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
        'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );
end


%% Compute some error statistics

if runORALS    % ORALS output includes ORSVD estimator
    estORALS.A_err           = graph_err ( estORALS.A, dyn_sys );
    estORALS.c_err           = kernel_err( estORALS.c, learning_setup );
    estORALS.Z_err           = get_Z_error( get_Z_from_E_c( estORALS.A, estORALS.c ), dyn_sys, learning_setup );
    [estORALS.pathTestErr,estORALS.estTestPathObj,estORALS.meanL2traj]    = getPathTestError( dyn_sys, learning_setup, estORALS.A, estORALS.c, testingPathsObj );
    ErrorORALS               = [ estORALS.c_err; estORALS.A_err; estORALS.Z_err; estORALS.pathTestErr; estORALS.time_ORALS];

    % SVD estimator error: it is estimated in ALS
    estORSVD.A_err           = graph_err ( estORALS.Esvd, dyn_sys );
    estORSVD.c_err           = kernel_err( estORALS.csvd, learning_setup );
    estORSVD.Z_err           = get_Z_error( get_Z_from_E_c( estORALS.Esvd, estORALS.csvd ), dyn_sys, learning_setup );
    estORSVD.time_ORSVD      = 1; % estORALS.time_build_B_tensor+estORALS.time_solveOR ;
    [estORSVD.pathTestErr,estORSVD.estTestPathObj,estORSVD.meanL2traj]    = getPathTestError( dyn_sys, learning_setup, estORALS.Esvd, estORALS.csvd, testingPathsObj );
    ErrorORSVD      = [ estORSVD.c_err; estORSVD.A_err; estORSVD.Z_err; estORSVD.pathTestErr; estORSVD.time_ORSVD];
end

if runALS
    estALS.A_err    = graph_err ( estALS.A, dyn_sys );
    estALS.c_err    = kernel_err( estALS.c, learning_setup );
    estALS.Z_err    = get_Z_error( get_Z_from_E_c( estALS.A, estALS.c ), dyn_sys, learning_setup );

    try
        estALS.ALS_error_kernel    = zeros(size(estALS.stats.c_hat,2),1);
        estALS.ALS_error_graph     = zeros(size(estALS.stats.c_hat,2),1);

        for q = 1:size(estALS.stats.c_hat,2)
            estALS.ALS_error_kernel(q) = kernel_err( estALS.stats.c_hat(:,q), learning_setup );
            estALS.ALS_error_graph(q)  = graph_err ( estALS.stats.A_hat(:,:,q), dyn_sys );
        end
    catch
    end
    [estALS.pathTestErr,estALS.estTestPathObj,estALS.meanL2traj]        = getPathTestError( dyn_sys, learning_setup, estALS.A, estALS.c, testingPathsObj );
    ErrorALS    = [ estALS.c_err; estALS.A_err; estALS.Z_err; estALS.pathTestErr;estALS.time.ALS];
end

%% Report estimation errors
clear ErrorTable; 
ErrorType   = ["Kernel";"Graph";"Tensor Z";"Paths";"Time"];
ErrorTable  = table(ErrorType);
if runALS;      ErrorTable = addvars(ErrorTable,ErrorALS); end
if runORALS;    ErrorTable = addvars(ErrorTable,ErrorORALS);  end   % ErrorTable = addvars(ErrorTable,ErrorORSVD); 

fprintf('\n');
disp(ErrorTable)

%% Plots
main_plots

