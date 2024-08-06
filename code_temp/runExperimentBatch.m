function runExperimentBatch( expSetup, paramStructs )

gcp;

running_tic = tic;
warning off;

%% Conduct one experiment for each combination of values of the parameter
for p = 1:length(paramStructs)
    p_idx = expSetup.p_perm( p );
    curFileName = sprintf('%s/LIKSonGraphs_largeTest_%i.mat',expSetup.saveDir,p_idx);
    if exist( curFileName, 'file')
        try
            tmp = load(curFileName,'curFileName');
            if numel(fieldnames(tmp))==0 || ~COMPLETE_PARTIAL_RESULTS
                fprintf('\nSkipping %d as found results from previous run in %s',p_idx,expSetup.saveDir);
                continue;
            end
        catch
            fprintf('\nFile %d seems corrupt: will re-run this experiment.',p_idx);
        end
    end
    save(curFileName,'curFileName','-v7.3');                                                                                    % create the file so other processes do not work on this example
    for pp = 1:expSetup.ntrials
        fprintf('\nRun %d.%d/%d : M=%d, L=%d, N=%d, n=%d, visc=%.2e, obs_std=%.2e, A_sp=%.2e', ...
            p,pp,length(paramStructs),paramStructs{p_idx}.M,paramStructs{p_idx}.L,paramStructs{p_idx}.N,paramStructs{p_idx}.n,paramStructs{p_idx}.viscosity,paramStructs{p_idx}.obs_std,paramStructs{p_idx}.A_sparsity);
        % Set up the experiment
        M                       = paramStructs{p_idx}.M;
        dyn_sys                 = paramStructs{p_idx};
        dyn_sys.A               = set_graph(dyn_sys.N, 'sparsity', dyn_sys.A_sparsity, 'plotON', 0);                            % create influence graph
        dyn_sys.initial         = 'Unif_0_5';
        dyn_sys.dt              = 1e-2;
        dyn_sys.T               = dyn_sys.dt*dyn_sys.L;
        dyn_sys                 = system_settings( dyn_sys );                                                                   % settings of the IPS and graph and its integrator

        learning_set            = learning_settings( 'randomSmoothFourierWithDecay',dyn_sys,struct('n',paramStructs{p_idx}.n) );
        learning_set.Z_true     = get_Z_from_E_c( dyn_sys.A, learning_set.c );                                                  % Z is the product of E and c

        dyn_sys.phi_kernel      = learning_set.phi_kernel_cheb;

        % Generate trajectories
        trainingPathsObj            = get_paths( dyn_sys, M,        'ParforProgressON', 0,'saveON', 0, 'loadON', 0 );
        testingPathsObj             = get_paths( dyn_sys, expSetup.M_test,   'ParforProgressON', 0,'saveON', 0, 'loadON', 0 );

        learning_set                = update_dict_mat( learning_set, trainingPathsObj.paths );                                  % compute the exploration measure rho and the basis matrix using all paths
        learning_set.kernel_norm    = sqrt( learning_set.c'*learning_set.dict_mat*learning_set.c );                             % this is the estimated kernel norm, assuming the kernel is in the hypothesis space

        %% Learning using ORALS
        matFactorize  = 'ALS';    % %  'SVD', 'ALS', 'SVD+ORALS'
        if isfield(expSetup,'runORALS_normaMat') && expSetup.runORALS_normalMat
            fprintf('\n\tRunning ORALS...')
            [estORALS0.Esvd, estORALS0.csvd, estORALS0.A, estORALS0.c, estORALS0.Z, condA_orals,~,estORALS0.time] = ...
                learn_kernel_graph_ORALS(trainingPathsObj.paths, dyn_sys, learning_set, 'plotON', 0, 'reg_method', 'lsqminnorm');         % reg_methods: ID, RKHS, None, lsqminnorm,pinv
        end
        if isfield(expSetup,'runORALS') && expSetup.runORALS
            estORALS = learn_kernel_graph_ORALS_B(trainingPathsObj.paths, dyn_sys, matFactorize, learning_set,...
                'plotON', 0, 'reg_method', 'lsqminnorm');         % reg_methods: ID, RKHS, None, lsqminnorm,pinv, pinvreg
        end

        %% Learning using alternating least squares
        if expSetup.runALS
            fprintf('\n\tRunning ALS...')
            estALS     = learn_kernel_graph_ALS( trainingPathsObj.paths, dyn_sys, learning_set, 'niter', 10, 'normalizeON', 1, ...
                'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
                'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
                'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );
            if ~expSetup.saveALSstats,   estALS=rmfield(estALS,'stats'); end
        end

        %% Compute some error statistics
        if isfield(expSetup,'runORALS_normaMat') && expSetup.runORALS_normalMat
            estORALS0.A_err         = graph_err ( estORALS0.A, dyn_sys );
            estORALS0.kernel_err    = kernel_err( estORALS0.c, learning_set );
            estORALS0.Z_err         = get_Z_error( get_Z_from_E_c( estORALS0.A, estORALS0.c ), dyn_sys, learning_set );
            [estORALS0.pathTestErr,estORALS0.estTestPathObj,estORALS0.meanL2traj]    = getPathTestError( dyn_sys, learning_set, estORALS0.A, estORALS0.c, testingPathsObj );
            ErrorORALS0             = [ estORALS0.c_err; estORALS0.A_err; estORALS0.Z_err; estORALS0.pathTestErr; estORALS0.time];
        end
        if isfield(expSetup,'runORALS') && expSetup.runORALS    % ORALS output includes ORSVD estimator
            estORALS.A_err          = graph_err ( estORALS.A, dyn_sys );
            estORALS.kernel_err     = kernel_err( estORALS.c, learning_set );
            estORALS.Z_err          = get_Z_error( get_Z_from_E_c( estORALS.A, estORALS.c ), dyn_sys, learning_set );
            [estORALS.pathTestErr,estORALS.estTestPathObj,estORALS.meanL2traj]    = getPathTestError( dyn_sys, learning_set, estORALS.A, estORALS.c, testingPathsObj );
            ErrorORALS              = [ estORALS.c_err; estORALS.A_err; estORALS.Z_err; estORALS.pathTestErr; estORALS.time];

            % SVD estimator error: it is estamated in ALS
            estORSVD.A_err          = graph_err ( estORALS.Esvd, dyn_sys );
            estORSVD.kernel_err     = kernel_err( estORALS.csvd, learning_set );
            estORSVD.Z_err          = get_Z_error( get_Z_from_E_c( estORALS.Esvd, estORALS.csvd ), dyn_sys, learning_set );
            estORSVD.time_ORSVD     = 1; % estORALS.time_build_B_tensor+estORALS.time_solveOR ;
            [estORSVD.pathTestErr,estORSVD.estTestPathObj,estORSVD.meanL2traj]    = getPathTestError( dyn_sys, learning_set, estORALS.Esvd, estORALS.csvd, testingPathsObj );
            ErrorORSVD              = [ estORSVD.c_err; estORSVD.A_err; estORSVD.Z_err; estORSVD.pathTestErr; estORSVD.time_ORSVD];
        end

        if expSetup.runALS
            estALS.A_err        = graph_err ( estALS.A, dyn_sys );
            estALS.kernel_err   = kernel_err( estALS.c, learning_set );
            estALS.Z_err        = get_Z_error( get_Z_from_E_c( estALS.A, estALS.c ), dyn_sys, learning_set );

            try
                estALS.ALS_error_kernel    = zeros(size(estALS.stats.c_hat,2),1);
                estALS.ALS_error_graph     = zeros(size(estALS.stats.c_hat,2),1);

                for q = 1:size(estALS.stats.c_hat,2)
                    estALS.ALS_error_kernel(q) = kernel_err( estALS.stats.c_hat(:,q), learning_set );
                    estALS.ALS_error_graph(q)  = graph_err ( estALS.stats.A_hat(:,:,q), dyn_sys );
                end
            catch
            end
            [estALS.pathTestErr,~,estALS.meanL2traj]    = getPathTestError( dyn_sys, learning_set, estALS.A, estALS.c, testingPathsObj );
            ErrorALS                                    = [ estALS.kernel_err; estALS.A_err; estALS.Z_err; estALS.pathTestErr;estALS.time.ALS];
        end

        %% Report estimation errors
        ErrorType   = ["Kernel";"Graph";"Tensor Z";"Paths";"Time"];
        ErrorTable  = table(ErrorType);
        if expSetup.runALS,                                                         ErrorTable = addvars(ErrorTable,ErrorALS); end
        if isfield(expSetup,'runORALS') && expSetup.runORALS,                       ErrorTable = addvars(ErrorTable,ErrorORALS); end
        if isfield(expSetup,'runORALS_normaMat') && expSetup.runORALS_normalMat,    ErrorTable = addvars(ErrorTable,ErrorORALS0); end

        fprintf('\n');
        disp(ErrorTable)

        %% Save information about this run
        cur_run.dyn_sys{pp}      = dyn_sys;
        cur_run.learning_set{pp} = learning_set;
        if expSetup.runALS,      cur_run.estALS{pp}      = estALS;   cur_run.ErrorALS{pp}    = ErrorALS;     end
        if isfield(expSetup,'runORALS') && expSetup.runORALS,    cur_run.estORALS{pp}    = estORALS; cur_run.ErrorORALS{pp}  = ErrorORALS;   end
        if isfield(expSetup,'runORALS_normaMat') && expSetup.runORALS_normalMat,   cur_run.estORALS0{pp}    = estORALS0; cur_run.ErrorORALS0{pp}  = ErrorORALS0;  end
    end
    cur_run.params  = paramStructs{p_idx};

    save(curFileName,'cur_run','p_idx','-v7.3');
    fprintf('\nEstimated time remaining: %e sec.',toc(running_tic)/p*(length(paramStructs)-p));
    close all;
end

return