% Convergence of graph and kernel error in sample size M

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
loadON = 1;     % Load previous results or not
saveON = 1;     % Save

% generate data of multple trajectories
num0 = (I.N-1+I.n); % /I.d; % /I.steps;
num1 = (I.N-1)*I.n; %/I.d; %/I.steps;

% Mseq consists of two parts.
% samll range and large range; merge together to compute
% Analyze separately
L = 8;   % length of Mseq in the large range; 
regu = 'RKHS'; % 'None'; %'lsqminnorm'; 'ID','RKHS'

Mseq_small = ceil(10.^linspace(log10(num0/2), log10(num1*4), L));
Mseq_large = ceil(10.^linspace(log10(num1*8), log10(10000), L));
Mseq = [Mseq_small, Mseq_large];

MM = Mseq(end)*3;  % Generate more data and sample from these data when testing convergence
fprintf('\nThe sequence of M is');disp(Mseq);

str_name      = sprintf('conv_in_M_N%i_n%i_M%i_L%i_obsnr',I.N,I.n,Mseq(end),I.steps);
str_viscosity = ['_visc',num2str(I.viscosity)]; 
str_name      = [str_name, num2str(I.obs_std),str_viscosity,'_ic',I.initial];

est_filename = [I.SAVE_DIR,str_name,'regu',regu,'.mat']; % saves estimator for each regu; 
rip_filename = [I.SAVE_DIR,str_name,'_rip','.mat'];      % 



%% compute the estimator; the different regularization method 
if ~exist(est_filename,'file') || ~loadON
    
    % generate data: and update I: exploration measure, basis matrix 
    fprintf('Generating trajectories, M = %i, and update I ...', MM);tic    % generate data if does not exist
    [all_xpath,I] = get_M_path(I, MM, 'ParforProgressON', 1,'saveON', 1,'loadON', 1); 
    fprintf('done (%.2f sec).\n',toc);
    
   % compute the RIP and condition number if file does not exist
    plotRIP = 1; 
    [RIP_seq,COND_seq]= compute_RIP_cond(rip_filename,all_xpath,I,Mseq,num0,num1,plotRIP); 

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
    comput_time = zeros(2,Mseq_length, test_num); 

    % compute the estimators 
    for b = 1:test_num                 % About parallel - parfor is used in ORALS & ALS, thus not here.         
        path_id   = randperm(MM, Mseq(end));
        test_path = all_xpath(path_id); % sample from a large collection of path
        for i = 1:Mseq_length
            fprintf('\n M-sequence:  %i out of %i : \n',i, Mseq_length);
            
            M = Mseq(i);
            [E_ORSVD, c_ORSVD, E_ORALS, c_ORALS, Z_OR, ~, condA_orals, ~, timeORALS] = learn_kernel_graph_ORALS(test_path(1:M), I, 'reg_method', regu);
            [E_ALS, c_ALS,~,~,timeALS]  = learn_kerne_graph_ALS(test_path(1:M), I, 'niter', 50, 'normalizeON', 1, 'plotON', 0,'reg_method', regu);
            
            all_c.ALS_seq{i, b} = c_ALS;
            all_c.ORALS_seq{i, b} = c_ORALS;
            all_c.ORSVD_seq{i, b} = c_ORSVD;
            
            all_E.ALS_seq{i, b} = E_ALS;
            all_E.ORALS_seq{i, b} = E_ORALS;
            all_E.ORSVD_seq{i, b} = E_ORSVD;

            comput_time(1,i,b) = timeORALS;
            comput_time(2,i,b) = timeALS; 
        end
    end
    save(est_filename,'all_E','all_c','condA_orals','graph_norm','kernel_norm','I','comput_time'); 
%     save(data_filename);
else
    load(est_filename,'all_E','all_c','condA_orals','graph_norm','kernel_norm','I','comput_time');
end


%% Compute error
for i = 1:Mseq_length
    for b = 1:test_num
        E_ALS   = all_E.ALS_seq{i, b}; 
        E_ORSVD = all_E.ORSVD_seq{i, b}; 
        E_ORALS = all_E.ORALS_seq{i, b};

        c_ALS   = all_c.ALS_seq{i, b}; 
        c_ORALS = all_c.ORALS_seq{i, b};
        c_ORSVD = all_c.ORSVD_seq{i, b};

        error.k_orals(i, b) = kernel_err(c_ORALS, I);
        error.k_orsvd(i, b) = kernel_err(c_ORSVD, I);
        error.k_als(i, b) = kernel_err(c_ALS, I);
        
        error.g_orals(i, b) = graph_err(E_ORALS, I);
        error.g_orsvd(i, b) = graph_err(E_ORSVD, I);
        error.g_als(i, b) = graph_err(E_ALS, I);
        

        Z_ALS   = get_Z_from_E_c(E_ALS, c_ALS);
        Z_ORSVD = get_Z_from_E_c(E_ORSVD, c_ORSVD);
        Z_ORALS = get_Z_from_E_c(E_ORALS, c_ORALS);
       
        error.Z_als(i, b)   = get_Z_error(Z_ALS, I);
        error.Z_orsvd(i, b) = get_Z_error(Z_ORSVD, I);
        error.Z_orals(i, b) = get_Z_error(Z_ORALS, I);
        error.Z_or(i, b)    = get_Z_error(Z_OR, I);
    end
end

% %% For small range of M
% ind = 1:L;
% ttl = 'Relative error as sample size increases (Large range)';
% draw_graph_kernel_M(Mseq, ind, I, error, ttl);
% % figname = [I.SAVE_DIR_fig,str_name,regu, 'small'];
% % set_positionFontsAll;
% 
% %% For Large range of M
% ind = L+1:2*L;
% ttl = 'Relative error as sample size increases (Large range)';
% draw_graph_kernel_M(Mseq, ind, I, error, ttl);
% % figname = [I.SAVE_DIR_fig,str_name, regu,'large'];
% % set_positionFontsAll;

%% All range of M
[~, ind] = sort(Mseq);
ttl = 'Relative error as sample size increases';
draw_graph_kernel_M(Mseq, ind, I, error, ttl);
figname = [I.SAVE_DIR_fig,str_name,'regu',regu,];
set_positionFontsAll;

plot_computTime(comput_time,Mseq,figname); 



% <<<<<<< Updated upstream
% =======
% figure;
% subplot(121);
% loglog(Mseq(ind), RIP_seq, '-', 'LineWidth', 3)
% title('Change of RIP ratio in ALS')
% yline(0.5)
% text(10,0.5,'\delta=0.5')
% ylim([min(0.48, min(RIP_seq)-0.05), max(RIP_seq)+0.1])
% xline(num0,'k','n0')
% xline(num1,'k','n1')
% % text(num0,10,'num0')
% % plot(log10(Mseq(ind)), log10(RIP_seq(ind)), '-')
% subplot(122);
% loglog(Mseq(ind), COND_seq(ind), '-', 'LineWidth', 3)
% title('Change of Conditional number in ORALS')
% xline(num1)
% text(num1,10,'n1')
% % xline(num0)
% % text(num0, 10,'num0')
% % plot(log10(Mseq(ind)), log10(COND_seq(ind)), '-')
% >>>>>>> Stashed changes
%%
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




