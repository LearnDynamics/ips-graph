% Convergence of graph and kernel error

% We wish to show the convergence as M increases

% There are two important threholds for the number M*d*steps
% (N-1)+n is the requirement for ALS
% (N-1)*n is the requirement for ORALS

% But we need also consider d and steps
% The mesh grid for the M is designed as following, so that 
% M*d*steps would cover the range ( (N-1)*n, Large)
% Hence M should be in the range ( (N-1)*n/(d*steps), Large)

% Taking into consideration the d and steps.
% the two important threshold become
% ((N-1)+n)/(d*steps)
% ((N-1)*n)/(d*steps)

%% Test Convergence in M for both ALS and ORALS
loadON = 0;     % Load previous results or not
saveON = 0;     % Save 

% generate data of multple trajectories
num0 = (dyn_sys.N-1+learning_setup.n)/dyn_sys.d/dyn_sys.L; 
num1 = (dyn_sys.N-1)*learning_setup.n/dyn_sys.d/dyn_sys.L;

Mseq_small = ceil(linspace(num0/4, num1*4, 10));
Mseq_large = ceil(10.^linspace(log10(num1/2), log10(5000), 10));
Mseq = sort([Mseq_small, Mseq_large],'ascend');

MM = Mseq(end)*10;                                                                                                              % Generate more data and sample from these data when testing convergence

fprintf('\nThe sequence of M is');disp(Mseq);

str_name = sprintf('conv_in_M_N%i_n%i_M%i_obsnr_',dyn_sys.N,learning_setup.n,Mseq(end));
str_name = [str_name, num2str(dyn_sys.obs_std),'regu'];
data_filename = [dyn_sys.SAVE_DIR,str_name,'.mat'];

if ~exist(data_filename,'file') || ~loadON
    
    fprintf('Generating trajectories, M = %i ...', MM);tic
    all_xpath = get_paths(dyn_sys, MM, 'ParforProgressON', 1,'saveON', 0,'loadON', 0);
    fprintf('done (%.2f sec).\n',toc);
    
    test_num = 1; % multiple tests  --- seems not necessary
    
    Mseq_length = length(Mseq);
    
    error_k_orals   = zeros(Mseq_length, test_num);
    error_k_orsvd   = zeros(Mseq_length, test_num);
    error_k_als     = zeros(Mseq_length, test_num);
    
    error_g_orals   = zeros(Mseq_length, test_num);
    error_g_orsvd   = zeros(Mseq_length, test_num);
    error_g_als     = zeros(Mseq_length, test_num);
    
    error_Z_orals   = zeros(Mseq_length, test_num);
    error_Z_orsvd   = zeros(Mseq_length, test_num);
    error_Z_als     = zeros(Mseq_length, test_num);
    % error_Z_or      = zeros(Mseq_length, test_num);
    
    c_ORALS_seq     = cell(Mseq_length, test_num);
    c_ORSVD_seq     = cell(Mseq_length, test_num);
    c_ALS_seq       = cell(Mseq_length, test_num);
    
    E_ORALS_seq     = cell(Mseq_length, test_num);
    E_ORSVD_seq     = cell(Mseq_length, test_num);
    E_ALS_seq       = cell(Mseq_length, test_num);
    
    learning_setup    = update_dict_mat(learning_setup, all_xpath.paths(1:MM));   
    
    
    %%

    estALS_all   = cell(Mseq_length,test_num);
    estORALS_all = cell(Mseq_length,test_num);
    for b = 1:test_num     % TBD: re-organize to run in parallel  - parfor is already used in ORALS, no need to use it here.
        path_id   = randperm(MM, Mseq(end));
        test_path = all_xpath.paths(path_id); % sample from a large collection of path
        for i = 1:Mseq_length
            fprintf('\n M-sequence:  %i out of %i : \n',i, Mseq_length);
            kernel_norm = sqrt(learning_setup.c'*learning_setup.dict_mat*learning_setup.c);
            graph_norm  = norm(dyn_sys.A,'fro');
            Z_true = get_Z_from_E_c(dyn_sys.A, learning_setup.c);

            M = Mseq(i);
            % [E_ORSVD, c_ORSVD, E_ORALS, c_ORALS, Z_OR,~,condA_orals,~,time_ORALS(i,b)] = learn_kernel_graph_ORALS(test_path(1:M), dyn_sys, learning_set, 'reg_method', 'pinv');

            % ORALS
            estORALS = learn_kernel_graph_ORALS_B(test_path(1:M), dyn_sys, 'ALS', learning_setup,...
                'plotON', 0, 'reg_method', 'lsqminnorm');         % reg_methods: ID, RKHS, None, lsqminnorm,pinv, pinvreg
            estORALS_all{i,b} = estORALS;
            E_ORALS = estORALS.A;     c_ORALS = estORALS.c;    time_ORALS(i,b) = estORALS.time_ORALS;
            E_ORSVD = estORALS.Esvd;  c_ORSVD = estORALS.csvd;

            % ALS
            % estALS  = learn_kernel_graph_ALS(test_path(1:M), dyn_sys, learning_set, 'niter', 50, 'stop_thres', dyn_sys.obs_std, 'normalizeON', 1, 'plotON', 0,'reg_methodK', 'pinv', 'reg_methodA', 'lsqnonneg','A_sparsity_thres', 0,'A_sparsity_prior',1);
            estALS     = learn_kernel_graph_ALS( test_path(1:M), dyn_sys, learning_setup, 'niter', 10, 'normalizeON', 1, ...
                'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
                'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
                'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );
            estALS_all{i,b} = estALS;
            E_ALS = estALS.A;         c_ALS = estALS.c;    time_ALS(i,b) = estALS.time.ALS;

            c_ALS_seq{i, b}     = c_ALS;
            c_ORALS_seq{i, b}   = c_ORALS;
            c_ORSVD_seq{i, b}   = c_ORSVD;

            E_ALS_seq{i, b}     = E_ALS;
            E_ORALS_seq{i, b}   = E_ORALS;
            E_ORSVD_seq{i, b}   = E_ORSVD;

            error_k_orals(i, b) = kernel_err( c_ORALS, learning_setup );
            error_k_orsvd(i, b) = kernel_err( c_ORSVD, learning_setup );
            error_k_als(i, b)   = kernel_err( c_ALS, learning_setup );

            error_g_orals(i, b) = graph_err( E_ORALS, dyn_sys );
            error_g_orsvd(i, b) = graph_err( E_ORSVD, dyn_sys );
            error_g_als(i, b)   = graph_err( E_ALS, dyn_sys );

            error_Z_als(i, b)   = get_Z_error( get_Z_from_E_c(E_ALS, c_ALS), dyn_sys, learning_setup );
            error_Z_orsvd(i, b) = get_Z_error( get_Z_from_E_c(E_ORSVD, c_ORSVD), dyn_sys, learning_setup );
            error_Z_orals(i, b) = get_Z_error( get_Z_from_E_c(E_ORALS, c_ORALS), dyn_sys, learning_setup );
            % error_Z_or(i, b)    = get_Z_error( Z_OR, dyn_sys, learning_set );
        end
    end

    try
    save(data_filename, 'error_k_orals', 'error_k_orsvd', 'error_k_als', 'error_g_orsvd','error_g_orals','error_g_als','error_Z_orsvd','error_Z_orals','error_Z_als','Mseq');% 'error_Z_or',
    catch %MM : save is buggy and fails, let's continue so computation is not wasted.
    end
else
    if loadON
        load(data_filename, 'error_k_orals', 'error_k_orsvd', 'error_k_als', 'error_g_orsvd','error_g_orals','error_g_als','error_Z_orsvd','error_Z_orals','error_Z_als','Mseq'); % 'error_Z_or',
    end
end

%% Plot graph and kernel estimation errors as a function of M
bigFig; 
loglog(Mseq,error_g_orals(:,1),'-o','color',[0.3010 0.7450 0.9330],'linewidth',2); hold on;
loglog(Mseq,error_k_orals(:,1),'-ob','linewidth',2); hold on;
loglog(Mseq,error_g_orsvd(:,1),'-o','color',[0.125 0.629 0.6940],'linewidth',2); hold on;
loglog(Mseq,error_k_orsvd(:,1),'-og','color',[0 0.4 0],'linewidth',2); hold on;
loglog(Mseq,error_g_als(:,1),'-o','color',[0.9290 0.6940 0.1250],'linewidth',2); hold on;
loglog(Mseq,error_k_als(:,1),'-or','linewidth',2); hold on;

xline(num1)
xline(num0)
xlabel('Sample size M'); ylabel('Relative L2 Error');
legend('Graph: ORALS','Kernel: ORALS','Graph: ORSVD','Kernel: ORSVD','Graph: ALS','Kernel: ALS');
legend('location','best');
title('Relative error as sample size increases')
axis tight; grid on;

figname = [dyn_sys.SAVE_DIR_fig,str_name];
set_positionFontsAll;

%% Plot Z estimation errors as a function of M
bigFig;hold on;
plot(log10(Mseq), log10(error_Z_als(:, 1)), '-x', 'LineWidth', 2, 'MarkerSize', 12);
plot(log10(Mseq), log10(error_Z_orsvd(:, 1)), '-o', 'LineWidth', 2, 'MarkerSize', 12);
plot(log10(Mseq), log10(error_Z_orals(:, 1)), '-.o', 'LineWidth', 2, 'MarkerSize', 12);
% plot(log10(Mseq), log10(error_Z_or(:, 1)), '--o', 'LineWidth', 2, 'MarkerSize', 12);
legend('ALS', 'ORSVD', 'ORALS', 'OR')
xlabel('log_{10} M')
ylabel('log_{10} Z error')
title('Compare the estimation of Z')
axis tight;grid on;
set_positionFontsAll;

%% Plot computation time as a function of M
bigFig;hold on;
loglog(Mseq, time_ALS(:, 1), '-or', 'LineWidth', 2, 'MarkerSize', 12);
loglog(Mseq, time_ORALS(:, 1), '-ob', 'LineWidth', 2, 'MarkerSize', 12);
legend({'ALS', 'ORALS'})
xlabel('log_{10} M')
ylabel('log_{10} computation time')
title('Compare computation')
axis tight;grid on;
set_positionFontsAll;


