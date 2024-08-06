%% test for regularization 

close all; clear all; addpaths;
rng(10); gcp;

SAVE_DIR = [getenv('HOME'),'/DataAnalyses/IPS_graph/'];
if ~exist(SAVE_DIR,'dir'); mkdir(SAVE_DIR); end
data_saving_folder = SAVE_DIR;


All_reg_methods = {'None', 'pinv', 'lsqminnorm', 'ID', 'RKHS'};
Mseq = 2.^(6:1:11);
nsimu    = 10; 

L = 6; N = 20; 
% kernel_type = 'randomSmoothFourierWithDecay';  kernel_type_str= 'randomFourier';
% kernel_type = 6; kernel_type_str = 'LJ6_'; %     % see options in learning_settings.m
kernel_type = 1; kernel_type_str = 'LJ1_'; %    % see options in learning_settings.m
kernel_type = 3; kernel_type_str = 'FS3_'; %    % see options in learning_settings.m


runALS                  = true;
runORALS                = true;   % using B tensor and long matrix

n_regSeq = length(All_reg_methods); 
nMseq    = length(Mseq); 


data_filename = [data_saving_folder,kernel_type_str, sprintf('regu_test_N%i_L%i.mat',N,L)];

if ~exist(data_filename,'file')
    % Generate a large dataset of M trajectories, then randomly sample a given size for multiple times.
    error_ALS   = zeros(5,nMseq,nsimu,n_regSeq);    %  [ c_err; A_err; Z_err; pathTestErr; time];
    error_ORALS = zeros(5,nMseq,nsimu,n_regSeq);
    condA_ORALS  = zeros(N,nMseq,nsimu,n_regSeq);   % condition number 
    condL_ORALS  = zeros(N,nMseq,nsimu,n_regSeq);
    cond_ALS    = zeros(2,nMseq,nsimu,n_regSeq);    %  only for kernel Est, since the lsqnonneg plays an important role for graph estimation

    M                       = 2*Mseq(end);    % number of trajectories
    dyn_sys.L               = L;      % number observations equispaced in time [0,T]

    % Settings for dynamics and learning
    dyn_sys.N               = N;                                                                                                   % number of agents
    dyn_sys.d               = 2;                                                                                                    % dim of state vectors
    dyn_sys.viscosity       = 1e-3;  %  set 0 to test consistency                                       % stochastic force; viscosity (forcing noise in the dynamics)
    dyn_sys.A               = set_graph(dyn_sys.N, 'sparsity', 0.4, 'plotON', 0);                                                   % create influence graph
    dyn_sys.initial         = 'Unif_0_5';
    dyn_sys.dt              = 1e-3;
    dyn_sys.T               = dyn_sys.dt*(dyn_sys.L-1);
    dyn_sys.obs_std         = 1e-2;  %  set 0 to test consistency                                       % observation noise
    dyn_sys                 = system_settings( dyn_sys );                                                                           % settings of the IPS and graph and its integrator


    % kernel_type             =  6; % 'randomSmoothFourierWithDecay'; % 6; %                  % see options in learning_settings.m
    n =  16;   %  it is determined in learning_settings, where kernel_type determines n
    learning_set            = learning_settings( kernel_type, dyn_sys, struct('n',n) );
    learning_set.Z_true     = get_Z_from_E_c( dyn_sys.A, learning_set.c );                                                          % Z is the product of E and c

    dyn_sys.phi_kernel      = learning_set.phi_kernel_cheb;                           % MM: logical inconsistency here and in the following because of this; learning_set.phi_kernel not used in any of what follows
    dyn_sys.phi_kernel      = learning_set.phi_kernel;                                % FL: learning_set.phi_kernel serves as the ground truth. --- It is a parametric inference. The chev is towards non-parametric inference
    dyn_sys.n               = learning_set.n;

    progressON  = false;    plotON      = false;

    % Generate training and testing paths
    M_test                      = 50;
    Z_true                      = get_Z_from_E_c(dyn_sys.A, learning_set.c);     % Z is the product of E and c
    fprintf('\nGenerating trajectories, M = %i...',M);tic
    trainingPathsObj            = get_paths( dyn_sys, M,        'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
    testingPathsObj             = get_paths( dyn_sys, M_test,   'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
    fprintf('done (%.2f sec).',toc);

    % update the exploration measure rho and the basis matrix using paths
    learning_set                = update_dict_mat( learning_set, trainingPathsObj.paths );                                          % compute the exploration measure rho and the basis matrix using all paths
    learning_set.kernel_norm    = sqrt( learning_set.c'*learning_set.dict_mat*learning_set.c );                                     % MM: this is the estimated kernel norm, assuming the kernel is in the hypothesis space


    for  ii = 1:nMseq
        M_train = Mseq(ii);
        for j= 1:n_regSeq
            regu_method= All_reg_methods{j};
            parfor nn =1:nsimu
                m_indx = randperm(M,M_train);
                %% Learning using ORALS
                % for kernel_type= 6: ID is bad; others are good: lsqminnorm, RKHS, None, pinv
                % for kernel_type= 'randomSmoothFourierWithDecay'; all good: ID, RKHS, None, pinv; but lsqminnorm is slower

                matFactorize  = 'ALS';    % %  'SVD', 'ALS', 'SVD+ORALS'
                if runORALS
                    estORALS = learn_kernel_graph_ORALS_B(trainingPathsObj.paths(m_indx), dyn_sys, matFactorize, learning_set,...
                        'plotON', 0, 'reg_method', regu_method);         % reg_methods: ID, RKHS, None, lsqminnorm,pinv, pinvreg
                    estORALS.A_err           = graph_err ( estORALS.A, dyn_sys );
                    estORALS.c_err           = kernel_err( estORALS.c, learning_set );
                    estORALS.Z_err           = get_Z_error( get_Z_from_E_c( estORALS.A, estORALS.c ), dyn_sys, learning_set );
                    [estORALS.pathTestErr,estORALS.estTestPathObj,estORALS.meanL2traj]    = getPathTestError( dyn_sys, learning_set, estORALS.A, estORALS.c, testingPathsObj );
                    ErrorORALS               = [ estORALS.c_err; estORALS.A_err; estORALS.Z_err; estORALS.pathTestErr; estORALS.time_ORALS];
                    error_ORALS(:,ii,nn,j)     = ErrorORALS;
                    condA_ORALS(:,ii,nn,j)    = estORALS.condA_orals';
                    condL_ORALS(:,ii,nn,j)    =  estORALS.condL';
                end
                %% Learning using alternating least squares
                if runALS
                    estALS     = learn_kernel_graph_ALS( trainingPathsObj.paths(m_indx), dyn_sys, learning_set, 'niter', 10, 'normalizeON', 1, ...
                        'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
                        'reg_methodK', regu_method, 'reg_methodA', 'lsqnonneg', ...
                        'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );

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
                    [estALS.pathTestErr,estALS.estTestPathObj,estALS.meanL2traj]        = getPathTestError( dyn_sys, learning_set, estALS.A, estALS.c, testingPathsObj );
                    ErrorALS    = [ estALS.c_err; estALS.A_err; estALS.Z_err; estALS.pathTestErr;estALS.time.ALS];
                    error_ALS(:,ii,nn,j)   = ErrorALS;
                    cond_ALS(:,ii,nn,j)    = estALS.condA_kernel([1,2]);
                end
            end
        end
        save(data_filename,'error_ALS','error_ORALS','cond_ALS','condA_ORALS','condL_ORALS','N','Mseq');
    end

end %% end if data 

%% present the results 
load(data_filename);
fig_dir = [SAVE_DIR,'figures/']; 
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
fig_str = [ fig_dir,'fig_regu',kernel_type_str,sprintf('_N%i_L%i',N,L)]; 
 
% the regularizers are helpful when the problem is ill-posed, e.g., when M is small or when the probem in the large sample limit is ill-conditioned 
    % % % error_ALS   = zeros(5,nMseq,nsimu,n_regSeq);  % [ c_err; A_err; Z_err; pathTestErr; time];
    % % % error_ORALS = zeros(5,nMseq,nsimu,n_regSeq);
    % % % condA_ORALS  = zeros(N,nMseq,nsimu,n_regSeq);
    % % % cond_ALS    = zeros(2,nMseq,nsimu,n_regSeq); 

% coef estimators 
     h = figure; % plot the error of regularized estimators 
    err_c_orals= squeeze(error_ORALS(1,1,:,:));  % nMseq x nsimu
    err_c_als   = squeeze(error_ALS(1,1,:,:));  % nMseq x nsimu
     
    subplot(1,2,1); boxplot(err_c_orals); ax = gca; ax.XTickLabel= All_reg_methods;  ylabel('ORALS: kernel error'); 
    subplot(1,2,2); boxplot(err_c_als); ax = gca; ax.XTickLabel= All_reg_methods; ylabel('ALS: kernel error'); 
    % sgtitle(['Error in c estimator, ',sprintf('M=%i',Mseq(1))])
    % exportgraphics(h, figname,'BackgroundColor','none')
    figname = [fig_str,sprintf('M%i_c',Mseq(1))];
    set_positionFontsAll; 


    figure; % plot the error of regularized estimators 
    err_c_orals= squeeze(error_ORALS(1,end,:,:));  % nMseq x nsimu
    err_c_als   = squeeze(error_ALS(1,end,:,:));  % nMseq x nsimu
     
    subplot(1,2,1); boxplot(err_c_orals); ax = gca; ax.XTickLabel= All_reg_methods; ylabel('ORALS: kernel error');
    subplot(1,2,2); boxplot(err_c_als); ax = gca; ax.XTickLabel= All_reg_methods; ylabel('ALS: kernel error');
   % sgtitle(['Error in c estimator, ',sprintf('M=%i',Mseq(2)) ])
    figname = [fig_str,sprintf('M%i_c',Mseq(end))];
    set_positionFontsAll; 

% A estimators 
          figure; % plot the error of regularized estimators 
    err_A_orals= squeeze(error_ORALS(2,1,:,:));  % nMseq x nsimu
    err_A_als   = squeeze(error_ALS(2,1,:,:));  % nMseq x nsimu
     
    subplot(1,2,1); boxplot(err_A_orals); ax = gca; ax.XTickLabel= All_reg_methods; ylabel('ORALS: graph error'); 
    subplot(1,2,2); boxplot(err_A_als); ax = gca; ax.XTickLabel= All_reg_methods;   ylabel('ALS: graph error');
    % sgtitle(['Error in A estimator, ',sprintf('M=%i',Mseq(1))])
    figname = [fig_str,sprintf('M%i',Mseq(1)),'_A'];
    set_positionFontsAll; 

           figure; % plot the error of regularized estimators 
    err_A_orals= squeeze(error_ORALS(2,end,:,:));  % nMseq x nsimu
    err_A_als   = squeeze(error_ALS(2,end,:,:));  % nMseq x nsimu
     
    subplot(1,2,1); boxplot(err_A_orals); ax = gca; ax.XTickLabel= All_reg_methods; ylabel('ORALS: graph error');
    subplot(1,2,2); boxplot(err_A_als); ax = gca; ax.XTickLabel= All_reg_methods;   ylabel('ALS: graph error');
    % sgtitle(['Error in A estimator, ',sprintf('M=%i',Mseq(end))])
    figname = [fig_str,sprintf('M%i_A',Mseq(end))];
    set_positionFontsAll; 

 % condition number    
    figure; % plot the condition number 
    mean_condA_orals= squeeze(mean(condA_ORALS(:,:,:,1),1)); % N,nMseq,nsimu,n_regSeq ==>>  % nMseq x nsimu
    mean_cond_als   = squeeze(mean(cond_ALS(2,:,:,1),1));  % nMseq x nsimu
     
    ymin = min(mean_condA_orals(end,:))*0.9; ymax = min(max(mean_condA_orals(1,:)),max(mean_condA_orals(2,:))*10);
    subplot(1,2,1); boxplot(mean_condA_orals');   ax = gca; ax.YAxis.Scale ="log";ylim([ymin,ymax]); ax.XTickLabel= Mseq;
                       ylabel('Condition number ORALS');  ax.XLabel.String ='Sample size M'; 
    subplot(1,2,2); boxplot(mean_cond_als');      ax = gca; ax.YAxis.Scale ="log";  ax.XTickLabel= Mseq; 
                       ylabel('Condition number: ALS -kernel'); ax.XLabel.String ='Sample size M'; 
    figname = [fig_str,sprintf('M%i_cond',Mseq(1))];
    set_positionFontsAll; 

