%% Convergence with respect to observation noise 
% @copyright: Quanjun Lang, Mauro Maggioni, Fei Lu
close all;  clear all;   clc;    addpaths;   rng(10);

%% System settings
% Specify the settings here to fix the figures for paper
I = [];
I = system_settings(I); % setting of the IPS and graph and its integrator
I.viscosity = 0;    % viscosity (noise)
I.N         = 6;               % number of agents
I.d         = 2;               % dim of state vectors
I.t0       = 0;
I.dt       = 1e-3;     % time step size
I.steps    = 1;        % time steps
I.obs_std  = 1e-4;       % 1e-4;     % observation noise
I.A = set_graph(I.N, 'sparsity', 0.4, 'plotON', 0);
I.initial = 'Unif_0_2';


PAPER_FIG_DIR = I.PAPER_FIG_DIR;        % saving folder for figures in paper

kernel_type             = 6;  % 6: for parametric inference
                              % 6 'typical_example_Lenard_Jones': reduandant-basis; 
                              % 'randomSmoothFourierWithDecay': many but good basis                % see options in learning_settings.m
n =  16;   %  it is determined in learning_settings, where kernel_type determines n

learning_setup            = learning_settings( kernel_type, I, struct('n',n) );
learning_setup.Z_true     = get_Z_from_E_c( I.A, learning_setup.c );       % Z is the product of A and c                                                      % Z is the product of E and c

% I.phi_kernel      = learning_setup.phi_kernel_cheb;                           
I.phi_kernel      = learning_setup.phi_kernel;                                

I.n               = learning_setup.n;


%% Test Convergence in observation noise for both ALS and ORALS
L = 7;obs_nse_seq = 10.^linspace(0, -6, L);
M = 1000;
MM = M*5;       % Generate more data and sample from them
test_num = 100;

% debug settings
% L = 3;obs_nse_seq = 10.^linspace(0, -6, L);
% M = 10;
% MM = M*5;       % Generate more data and sample from them
% test_num = 5;



regu = 'lsqminnorm';       % reg_methods: ID, RKHS, None, lsqminnorm,pinv, pinvreg
 fprintf('regu Method for ALS and ORALS = %s,   ...\n ',regu);

loadON = 0;     % Load previous results or not
saveON = 1;     % Save

fprintf('\nThe sequence of observation noise is');
for i = 1:L;fprintf('\t%d', obs_nse_seq(i));end;fprintf('\n');

str_name      = sprintf('conv_in_obsnse_N%i_kernelType%i_d%i_M%i_L%i_vis%i_', I.N, kernel_type, I.d, M, I.steps, I.viscosity);
str_name      = [str_name,'_ic',I.initial, '_regu_', regu, '_nserange_', num2str(min(obs_nse_seq)), '_', num2str(max(obs_nse_seq)), '_L', num2str(L), 'test_num', num2str(test_num)];

est_filename = [I.SAVE_DIR,str_name,'regu',regu,'.mat']; % saves estimator for each regu;
rip_filename = [I.SAVE_DIR,str_name,'_rip','.mat'];      %

fprintf('\n The file name is : \n');
disp(str_name);
%% compute the estimator; the different regularization method
if ~exist(est_filename,'file') || ~loadON
    
    all_c.ORALS_seq = cell(L, test_num);
    all_c.ORSVD_seq = cell(L, test_num);
    all_c.ALS_seq   = cell(L, test_num);
    
    all_E.ORALS_seq = cell(L, test_num);
    all_E.ORSVD_seq = cell(L, test_num);
    all_E.ALS_seq   = cell(L, test_num);
    
    time_ORALS = zeros(L, test_num);
    time_ALS   = zeros(L, test_num);

    % compute the estimators
    for i = 1:L
        fprintf('\nObservation Noise sequence:  %i out of %i : \n', i, L);
        % generate data
        I.obs_std = obs_nse_seq(i);
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

        
        for b = 1:test_num
            path_id   = randperm(MM, M);
            test_path = all_xpath(path_id);         % sample from a large collection of path

            % ORALS
            estORALS = learn_kernel_graph_ORALS_B(test_path, I, 'ALS', learning_setup,...
                'plotON', 0, 'reg_method',regu);         % reg_methods: ID, RKHS, None, lsqminnorm,pinv, pinvreg
            E_ORALS = estORALS.A;     c_ORALS = estORALS.c;    time_ORALS(i,b) = estORALS.time_ORALS;
            E_ORSVD = estORALS.Esvd;  c_ORSVD = estORALS.csvd;

            % ALS
            estALS     = learn_kernel_graph_ALS( test_path, I, learning_setup, 'niter', 10, 'normalizeON', 1, ...
                'stop_thres_testpaths', I.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 0, ...
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
    % Compute error
    for i = 1:L
        for b = 1:test_num
            E_ALS   = all_E.ALS_seq{i, b};
            E_ORSVD = all_E.ORSVD_seq{i, b};
            E_ORALS = all_E.ORALS_seq{i, b};

            c_ALS   = all_c.ALS_seq{i, b};
            c_ORALS = all_c.ORALS_seq{i, b};
            c_ORSVD = all_c.ORSVD_seq{i, b};

            error.k_orals(i, b) = kernel_err(c_ORALS, learning_setup);
            error.k_orsvd(i, b) = kernel_err(c_ORSVD, learning_setup);
            error.k_als(i, b) = kernel_err(c_ALS, learning_setup);

            error.g_orals(i, b) = graph_err(E_ORALS, I);
            error.g_orsvd(i, b) = graph_err(E_ORSVD, I);
            error.g_als(i, b) = graph_err(E_ALS, I);
        end
    end

    save(est_filename,'all_E','all_c', 'I','error',"obs_nse_seq",'learning_setup');
else
    load(est_filename,'all_E','all_c', 'I','error',"obs_nse_seq",'learning_setup');
end


%% Boxplot

k_y_max = max([error.k_orals; error.k_als], [], 'all');
k_y_min = min([error.k_orals; error.k_als], [], 'all');

g_y_min = min([error.g_orals; error.g_als], [], 'all');
g_y_max = max([error.g_orals; error.g_als], [], 'all');
% 


figure; 
tiledlayout(1, 4, 'TileSpacing', 'normal', 'padding', 'compact');
nexttile;
% subplot(141);
boxplot(error.g_als',  obs_nse_seq);  xlabel('Observation Noise');ylabel('Graph error');title('ALS');  grid on;
hold on;plot(L:-1:1, mean(error.g_als, 2))
ylim([g_y_min, g_y_max]); 
ax = gca;ax.YAxis.Scale = "log"; 
set(gca, 'xdir', 'reverse', 'xtick',1:L, 'Xticklabel', {'$10^{-6}$', '', '$10^{-4}$', '', '$10^{-2}$', '', '$10^{0}$'}, 'TickLabelInterpreter', 'latex');


% subplot(142);
nexttile;
boxplot(error.g_orals', obs_nse_seq); xlabel('Observation Noise');ylabel('Graph error');title('ORALS');grid on
hold on;plot(L:-1:1, mean(error.g_orals, 2))
ylim([g_y_min, g_y_max]); ax = gca;ax.YAxis.Scale = "log";
set(gca, 'xdir', 'reverse', 'xtick',1:L, 'Xticklabel', {'$10^{-6}$', '', '$10^{-4}$', '', '$10^{-2}$', '', '$10^{0}$'}, 'TickLabelInterpreter', 'latex');


% subplot(143);
nexttile;
boxplot(error.k_als',  obs_nse_seq);  xlabel('Observation Noise');ylabel('Kernel error');title('ALS');  grid on
hold on;plot(L:-1:1, mean(error.k_als, 2))
ylim([k_y_min, k_y_max]); ax = gca;ax.YAxis.Scale = "log";  
set(gca, 'xdir', 'reverse', 'xtick',1:L, 'Xticklabel', {'$10^{-6}$', '', '$10^{-4}$', '', '$10^{-2}$', '', '$10^{0}$'}, 'TickLabelInterpreter', 'latex');


% subplot(144);
nexttile;
boxplot(error.k_orals', obs_nse_seq); xlabel('Observation Noise');ylabel('Kernel error');title('ORALS');grid on
hold on;plot(L:-1:1, mean(error.k_orals, 2))
ylim([k_y_min, k_y_max]); ax = gca;ax.YAxis.Scale = "log";
set(gca, 'xdir', 'reverse', 'xtick',1:L, 'Xticklabel', {'$10^{-6}$', '', '$10^{-4}$', '', '$10^{-2}$', '', '$10^{0}$'}, 'TickLabelInterpreter', 'latex');

%set(gcf,'Position',[100 100 1200 350])
 set(gcf,'Position',[100 100 800 250])
set(findall(gcf,'-property','FontSize'),'FontSize', 13)
% tightfig
 set_positionFontsAll;

%% save figure to paper figure folder
figname = [PAPER_FIG_DIR, '/convergence_obs_std.pdf'];
saveas(gcf, figname);