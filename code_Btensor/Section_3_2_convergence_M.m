% Convergence with respect to M 
close all;  clear all;   clc;    addpaths;   rng(10);

%% System settings
% Specify the settings here to fix the figures for paper
I = [];
I = system_settings(I); % setting of the IPS and graph and its integrator
I.viscosity = 1e-4;    % viscosity (noise)
I.N         = 6;               % number of agents
I.d         = 2;               % dim of state vectors
I.t0       = 0;
I.dt       = 1e-3;     % time step size
I.steps    = 1;        % time steps
I.obs_std  = 1e-4;       % 1e-4;     % observation noise
I.A = set_graph(I.N, 'sparsity', 0.4, 'plotON', 0);
I.initial = 'Unif_0_2';
% I.basis_case = 1;
% I = update_system_settings(I);         % this part to be changed 

PAPER_FIG_DIR = I.PAPER_FIG_DIR;        % saving folder for figures in paper

kernel_type             = 'typical_example_Lenard_Jones';  % 6: for parametric inference
                              % 6 'typical_example_Lenard_Jones': reduandant-basis; 
                              % 'randomSmoothFourierWithDecay': many but good basis                % see options in learning_settings.m
n =  16;   %  it is determined in learning_settings, where kernel_type determines n

learning_setup            = learning_settings( kernel_type, I, struct('n',n) );
learning_setup.Z_true     = get_Z_from_E_c( I.A, learning_setup.c );       % Z is the product of A and c                                                      % Z is the product of E and c

I.phi_kernel      = learning_setup.phi_kernel_cheb;                           % MM: logical inconsistency here and in the following because of this; learning_set.phi_kernel not used in any of what follows
%I.phi_kernel      = learning_setup.phi_kernel;                                % FL: learning_set.phi_kernel serves as the ground truth. --- It is a parametric inference. The chev is towards non-parametric inference

I.n               = learning_setup.n;


%% Test Convergence in M for both ALS and ORALS
L_M = 10;
M_seq = floor(10.^linspace(1, 4.5, L_M));

MM = M_seq(end)*2;       % Generate more data and sample from them
test_num = 100;


%% Debug test setting
% L_M = 5;
% M_seq = floor(10.^linspace(1, 2, L_M));
% 
% MM = M_seq(end)*2;       % Generate more data and sample from them
% test_num = 5;

%%

regu = 'lsqminnorm';       % reg_methods: ID, RKHS, None, lsqminnorm,pinv, pinvreg

loadON = 0;     % Load previous results or not
saveON = 1;     % Save

fprintf('\nThe sequence of M is');
for i = 1:L_M;fprintf('\t%d', M_seq(i));  end;   fprintf('\n');

str_name      = sprintf('conv_in_M_N%i_kernelType%s_d%i_L%i_vis%i_obsstd%i', I.N,kernel_type, I.d, I.steps, I.viscosity, I.obs_std);
str_name      = [str_name,'_ic',I.initial, '_regu_', regu, 'test_num_', num2str(test_num)];


est_filename = [I.SAVE_DIR,str_name,'regu',regu,'.mat']; % saves estimator for each regu;
rip_filename = [I.SAVE_DIR,str_name,'_rip','.mat'];      %

fprintf(['File name is : ', str_name, '\n']);


%% compute the estimator; the different regularization method
if ~exist(est_filename, 'file') || ~loadON

    %% % Generate data with largest M
    fprintf('Generating trajectories, M = %i, and update I ...', MM);tic
    pathObj   = get_paths(I, MM, 'ParforProgressON', 1,'saveON', 0,'loadON', 0);   % saveON = 1 is not recommended, since data can be efficiently generated
    all_xpath = pathObj.paths;
    fprintf('done (%.2f sec).\n',toc);

    % update the exploration measure rho and the basis matrix using paths
    learning_setup                = update_dict_mat( learning_setup, all_xpath);                                          % compute the exploration measure rho and the basis matrix using all paths
    if learning_setup.basis_changed ==1
        learning_setup.Z_true_orig    = learning_setup.Z_true;        % just in case a change of basis is applied when dict_mat is singular
        learning_setup.Z_true         = get_Z_from_E_c( dyn_sys.A, learning_setup.c );
    end

    learning_setup.kernel_norm    = sqrt( learning_setup.c'*learning_setup.dict_mat*learning_setup.c );


    all_c.ORALS_seq = cell(L_M, test_num);
    all_c.ORSVD_seq = cell(L_M, test_num);
    all_c.ALS_seq   = cell(L_M, test_num);  

    all_E.ORALS_seq = cell(L_M, test_num);
    all_E.ORSVD_seq = cell(L_M, test_num);
    all_E.ALS_seq   = cell(L_M, test_num); 
    time_ORALS = zeros(L_M, test_num);
    time_ALS   = zeros(L_M, test_num);

    for i = 1:L_M
        fprintf('\nM sequence:  %i out of %i : \n', i, L_M);
        M = M_seq(i);

        for b = 1:test_num
            path_id   = randperm(MM, M);
            test_path = all_xpath(path_id);         % sample from a large collection of path

            % ORALS
%            [E_ORSVD, c_ORSVD, E_ORALS, c_ORALS, Z_OR, ~, ~, ~, ~] = learn_kernel_graph_ORALS(test_path, dyn_sys, 'reg_method', regu);
            estORALS = learn_kernel_graph_ORALS_B(test_path, I, 'ALS', learning_setup,...
                'plotON', 0, 'reg_method',regu);         % reg_methods: ID, RKHS, None, lsqminnorm,pinv, pinvreg
            E_ORALS = estORALS.A;     c_ORALS = estORALS.c;    time_ORALS(i,b) = estORALS.time_ORALS;
            E_ORSVD = estORALS.Esvd;  c_ORSVD = estORALS.csvd;

            % ALS
%            [E_ALS, c_ALS, ~, ~, ~]  = learn_kerne_graph_ALS(test_path, dyn_sys, 'niter', 100, 'normalizeON', 1, 'plotON', 0,'reg_method', regu);
            % estALS  = learn_kernel_graph_ALS(test_path(1:M), dyn_sys, learning_set, 'niter', 50, 'stop_thres', dyn_sys.obs_std, 'normalizeON', 1, 'plotON', 0,'reg_methodK', 'pinv', 'reg_methodA', 'lsqnonneg','A_sparsity_thres', 0,'A_sparsity_prior',1);
            estALS     = learn_kernel_graph_ALS( test_path, I, learning_setup, 'niter', 10, 'normalizeON', 1, ...
                'stop_thres_testpaths', I.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
                'reg_methodK', regu, 'reg_methodA', 'lsqnonneg', ...
                'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );
            E_ALS = estALS.A;         c_ALS = estALS.c;    time_ALS(i,b) = estALS.time.ALS;


            all_c.ALS_seq{i, b} = c_ALS;
            all_c.ORALS_seq{i, b} = c_ORALS;
            all_c.ORSVD_seq{i, b} = c_ORSVD;

            all_E.ALS_seq{i, b} = E_ALS;
            all_E.ORALS_seq{i, b} = E_ORALS;
            all_E.ORSVD_seq{i, b} = E_ORSVD;
        end
    end
    save(est_filename,'all_E','all_c', 'I',"time_ORALS","learning_setup","time_ALS");
else
    load(est_filename,'all_E','all_c', 'I',"time_ORALS","learning_setup","time_ALS");
end

%% Compute error
for i = 1:size(all_E.ALS_seq,1)
    for b = 1:size(all_E.ALS_seq,2)
        E_ALS   = all_E.ALS_seq{i, b};
        E_ORSVD = all_E.ORSVD_seq{i, b};
        E_ORALS = all_E.ORALS_seq{i, b};

        c_ALS   = all_c.ALS_seq{i, b};
        c_ORALS = all_c.ORALS_seq{i, b};
        c_ORSVD = all_c.ORSVD_seq{i, b};

        error.k_orals(i, b) = kernel_err(c_ORALS, learning_setup);
        error.k_orsvd(i, b) = kernel_err(c_ORSVD,learning_setup );
        error.k_als(i, b)   = kernel_err(c_ALS, learning_setup);

        error.g_orals(i, b) = graph_err(E_ORALS, I);
        error.g_orsvd(i, b) = graph_err(E_ORSVD, I);
        error.g_als(i, b)   = graph_err(E_ALS, I);

    end
end

%% plot
% figure;
% mksz = 13;
k_y_max = max([error.k_orals; error.k_als], [], 'all');
k_y_min = min([error.k_orals; error.k_als], [], 'all');

g_y_min = min([error.g_orals; error.g_als], [], 'all');
g_y_max = max([error.g_orals; error.g_als], [], 'all');

% Boxplot
warning('off', 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout')

y_min = min(g_y_min, k_y_min);
y_max = max(g_y_max, k_y_max);

M_ALS   = I.n+I.N^2;      str_M_ALS= 'N^2+n'; 
M_ORALS = I.n*(I.N)^2;    str_M_ORALS= 'N^2\times n';

% M_ALS   = I.n+I.N;      str_M_ALS= 'N+n'; 
% M_ORALS = I.n*(I.N);    str_M_ORALS= 'N\times n';

critical_M_on = 0; 

ff = figure;
% subplot(141);
tiledlayout(1, 4, 'TileSpacing', 'normal', 'padding', 'compact');
nexttile;
boxplot(error.g_als',  log10(M_seq));  xlabel('M');ylabel('Graph error');title('ALS');  grid on
hold on;%plot(1:L_M, trimmean(error.g_als, 95, 2))
ylim([g_y_min, g_y_max]); 
ax = gca;ax.YAxis.Scale = "log";
if critical_M_on==1
    xline(log10(M_ALS));     xline(log10(M_ORALS)); 
    text(log10(M_ALS)-1.6, 0.03, str_M_ALS)
    text(log10(M_ORALS)+0.2, 0.03, str_M_ORALS)
end
set(gca, 'xtick',1:9/7:10, 'XTickLabels', {'$10^1$', '', '$10^2$', '', '$10^3$', '', '$10^4$', ''} , 'TickLabelInterpreter', 'latex')
% set(gca, 'xtick',1:18/7:10, 'XTickLabels', compose('$10^%d$',0:3), 'TickLabelInterpreter', 'latex')



% subplot(142);
nexttile;
boxplot(error.g_orals', M_seq); xlabel('M');% ylabel('Graph error');
title('ORALS');grid on
hold on;
plot(1:L_M, trimmean(error.g_orals, 95, 2))
ylim([g_y_min, g_y_max]); 
ax = gca;ax.YAxis.Scale = "log";
if critical_M_on==1
    xline(log10(M_ALS));     xline(log10(M_ORALS)); 
    text(log10(M_ALS)-1.6, 0.03, str_M_ALS)
    text(log10(M_ORALS)+0.2, 0.03, str_M_ORALS)
end
set(gca, 'xtick',1:9/7:10, 'XTickLabels', {'$10^1$', '', '$10^2$', '', '$10^3$', '', '$10^4$', ''} , 'TickLabelInterpreter', 'latex')


% subplot(143);
nexttile;
boxplot(error.k_als',  M_seq);  xlabel('M');ylabel('Kernel error');title('ALS');  grid on
hold on;
plot(1:L_M, trimmean(error.k_als, 95, 2))
ylim([k_y_min, k_y_max]); 
ax = gca;ax.YAxis.Scale = "log";
if critical_M_on==1
    xline(log10(M_ALS));     xline(log10(M_ORALS)); 
    text(log10(M_ALS)-1.6, 0.03, str_M_ALS)
    text(log10(M_ORALS)+0.2, 0.03, str_M_ORALS)
end
set(gca, 'xtick',1:9/7:10, 'XTickLabels', {'$10^1$', '', '$10^2$', '', '$10^3$', '', '$10^4$', ''} , 'TickLabelInterpreter', 'latex')


%subplot(144);
nexttile;
boxplot(error.k_orals', M_seq); xlabel('M');% ylabel('Kernel error');
title('ORALS');grid on
hold on;
plot(1:L_M, trimmean(error.k_orals, 95, 2))
ylim([k_y_min, k_y_max]); 
ax = gca;ax.YAxis.Scale = "log";
if critical_M_on==1
    xline(log10(M_ALS));     xline(log10(M_ORALS)); 
    text(log10(M_ALS)-1.6, 0.03, str_M_ALS)
    text(log10(M_ORALS)+0.2, 0.03, str_M_ORALS)
end
set(gca, 'xtick',1:9/7:10, 'XTickLabels', {'$10^1$', '', '$10^2$', '', '$10^3$', '', '$10^4$', ''} , 'TickLabelInterpreter', 'latex')

set(gcf,'Position',[100 100 1000 250])
% set(findall(gcf,'-property','FontSize'),'FontSize', 13)                                                                        MM: This is buggy, it completely changes the scale of the Y-axis in some plots (e.g. the last one): 
set_positionFontsAll;

% set_positionFontsAll;
% tightfig
% set(ff,'PaperSize',[400, 100])

figname = [PAPER_FIG_DIR,str_name, '.pdf'];
saveas(gcf, figname);

