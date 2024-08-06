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
saveON = 1;     % Save

% generate data of multple trajectories
num0 = (I.N-1+I.n)/I.d/I.steps;
num1 = (I.N-1)*I.n/I.d/I.steps;

% Mseq consists of two parts.
% samll range and large range
% Each has 10 steps. merge together to compute
% Analyze separately
L = 8;
Mseq_small = ceil(10.^linspace(log10(num0/4), log10(num0*4), L));
Mseq_large = ceil(10.^linspace(log10(num1/2), log10(2000), L));
Mseq = [Mseq_small, Mseq_large];

MM = Mseq(end)*3;  % Generate more data and sample from these data when testing convergence

fprintf('\nThe sequence of M is');disp(Mseq);

str_name = sprintf('conv_in_M_N%i_n%i_M%i_obsnr_',I.N,I.n,Mseq(end));
str_name = [str_name, num2str(I.obs_std),'regu'];
data_filename = [I.SAVE_DIR,str_name,'.mat'];

if ~exist(data_filename,'file') || ~loadON
    
    fprintf('Generating trajectories, M = %i...', MM);tic
    all_xpath = get_M_path(I, MM, 'ParforProgressON', 1,'saveON', 1,'loadON', 0);
    fprintf('done (%.2f sec).\n',toc);
    
    %% FL: to compute the exploration measure rho and the basis matrix here: using all paths
    I = update_dict_mat(I, all_xpath); % This is it
    
    test_num = 1; % multiple tests  --- seems not necessary
    Mseq_length = length(Mseq);
    
    kernel_norm = sqrt(I.c_true'*I.dict_mat*I.c_true);
    graph_norm  = norm(I.E,'fro');
    
    
    all_c.ORALS_seq = cell(Mseq_length, test_num);
    all_c.ORSVD_seq = cell(Mseq_length, test_num);
    all_c.ALS_seq = cell(Mseq_length, test_num);
    
    all_E.ORALS_seq = cell(Mseq_length, test_num);
    all_E.ORSVD_seq = cell(Mseq_length, test_num);
    all_E.ALS_seq = cell(Mseq_length, test_num);

    regu = 'RKHS'; % 'None'; %'lsqminnorm'; 'ID','RKHS'
    
    %%
    for b = 1:test_num
        % TBD: re-organize to run in parallel  - parfor is already used in ORALS, no need to use it here.
        % OK
        path_id   = randperm(MM, Mseq(end));
        test_path = all_xpath(path_id); % sample from a large collection of path
        for i = 1:Mseq_length
            fprintf('\n M-sequence:  %i out of %i : \n',i, Mseq_length);
            
            M = Mseq(i);
            [E_ORSVD, c_ORSVD, E_ORALS, c_ORALS, Z_OR, ~, condA_orals] = learn_kernel_graph_ORALS(test_path(1:M), I, 'reg_method', regu);
            [E_ALS, c_ALS]  = learn_kerne_graph_ALS(test_path(1:M), I, 'niter', 50, 'normalizeON', 1, 'plotON', 0,'reg_method', regu);
            
            all_c.ALS_seq{i, b} = c_ALS;
            all_c.ORALS_seq{i, b} = c_ORALS;
            all_c.ORSVD_seq{i, b} = c_ORSVD;
            
            all_E.ALS_seq{i, b} = E_ALS;
            all_E.ORALS_seq{i, b} = E_ORALS;
            all_E.ORSVD_seq{i, b} = E_ORSVD;
        end
    end
    save(data_filename);
else
    if loadON
        load(data_filename);
    end
end


%% Compute error
for i = 1:Mseq_length
    for b = 1:test_num
        error.k_orals(i, b) = kernel_err(all_c.ORALS_seq{i, b}, I);
        error.k_orsvd(i, b) = kernel_err(all_c.ORSVD_seq{i, b}, I);
        error.k_als(i, b) = kernel_err(all_c.ALS_seq{i, b}, I);
        
        error.g_orals(i, b) = graph_err(all_E.ORALS_seq{i, b}, I);
        error.g_orsvd(i, b) = graph_err(all_E.ORSVD_seq{i, b}, I);
        error.g_als(i, b) = graph_err(all_E.ALS_seq{i, b}, I);
        
        Z_ALS   = get_Z_from_E_c(E_ALS, c_ALS);
        Z_ORSVD = get_Z_from_E_c(E_ORSVD, c_ORSVD);
        Z_ORALS = get_Z_from_E_c(E_ORALS, c_ORALS);
        
        
        error.Z_als(i, b)   = get_Z_error(Z_ALS, I);
        error.Z_orsvd(i, b) = get_Z_error(Z_ORSVD, I);
        error.Z_orals(i, b) = get_Z_error(Z_ORALS, I);
        error.Z_or(i, b)    = get_Z_error(Z_OR, I);
    end
end

%% For small range of M
ind = 1:L;
ttl = 'Relative error as sample size increases (Large range)';
draw_graph_kernel_M(Mseq, ind, I, error, ttl);
% figname = [I.SAVE_DIR_fig,str_name,regu, 'small'];
% set_positionFontsAll;

%% For Large range of M
ind = L+1:2*L;
ttl = 'Relative error as sample size increases (Large range)';
draw_graph_kernel_M(Mseq, ind, I, error, ttl);
% figname = [I.SAVE_DIR_fig,str_name, regu,'large'];
% set_positionFontsAll;

%% All range of M
[~, ind] = sort(Mseq);
ttl = 'Relative error as sample size increases (All)';
draw_graph_kernel_M(Mseq, ind, I, error, ttl);
figname = [I.SAVE_DIR_fig,str_name, regu];
set_positionFontsAll;

%%
%
%
% %%
% figure;hold on;
% plot(log10(Mseq), log10(error_Z_als(:, 1)), '-x', 'LineWidth', 2, 'MarkerSize', 12);
% plot(log10(Mseq), log10(error_Z_orsvd(:, 1)), '-o', 'LineWidth', 2, 'MarkerSize', 12);
% plot(log10(Mseq), log10(error_Z_orals(:, 1)), '-.o', 'LineWidth', 2, 'MarkerSize', 12);
% plot(log10(Mseq), log10(error_Z_or(:, 1)), '--o', 'LineWidth', 2, 'MarkerSize', 12);
% legend('ALS', 'ORSVD', 'ORALS', 'OR')
% xlabel('log_{10} M')
% ylabel('log_{10} Z error')
% title('Compare the estimation of Z')
%
% % str_name = sprintf('Zconv_in_M_N%i_n%i_M%i',I.N,I.n,Mseq(end));
% % figname = [I.SAVE_DIR_fig,str_name];
% % set_positionFontsAll;
%
%




function draw_graph_kernel_M(Mseq, ind, I, error, ttl)
num0 = (I.N-1+I.n)/I.d/I.steps;
num1 = (I.N-1)*I.n/I.d/I.steps;

figure; % plot error of ORALS in M
loglog(Mseq(ind),error.g_orals(ind,1),'-x','linewidth',1); hold on;
loglog(Mseq(ind),error.k_orals(ind,1),'-o','linewidth',1); hold on;
loglog(Mseq(ind),error.g_orsvd(ind,1),'--x','linewidth',1.5); hold on;
loglog(Mseq(ind),error.k_orsvd(ind,1),'--o','linewidth',1.5); hold on;
loglog(Mseq(ind),error.g_als(ind,1),'-.x','linewidth',1.5); hold on;
loglog(Mseq(ind),error.k_als(ind,1),'-.o','linewidth',1.5); hold on;
grid on
xline(num1)
text(num1,I.obs_std/1.3,'num1')
xline(num0)
text(num0,I.obs_std/1.3,'num0')
yline(I.obs_std);
text(min(Mseq(ind)),I.obs_std*1.3,'obs std')
ylim([I.obs_std/2, 10])
xlabel('Sample size M'); ylabel('Relative L2 Error');
legend('Graph: ORALS','Kernel: ORALS','Graph: ORSVD','Kernel: ORSVD','Graph: ALS','Kernel: ALS');
legend('location','best');
title(ttl)

end