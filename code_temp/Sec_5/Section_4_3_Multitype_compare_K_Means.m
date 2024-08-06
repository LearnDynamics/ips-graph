% Section 4.3 Multitype kernel
% Use code based on B tensors
% Use 3 way ALS to learn multitype kernel

% Quanjun Lang (c) 2024



close all; clear all; addpaths;clc

rng(20);

% gcp;        % get current parallel pool

%% Observations and training data sizes
M                       = 400;                                                 % number of trajectories
dyn_sys.L               = 20;      % number observations equispaced in time [0,T]
runALS                  = true;
runORALS                = true;   % using B tensor and long matrix
runORALS_normalMat      = true;   % without using B tensor, use normal matrices 

% Settings for dynamics and learning

dyn_sys.N               = 10;                                                                                                   % number of agents
dyn_sys.d               = 2;                                                                                                    % dim of state vectors
dyn_sys.viscosity       = 1e-3;  %  set 0 to test consistency                                       % stochastic force; viscosity (forcing noise in the dynamics)
dyn_sys.A               = set_graph(dyn_sys.N, 'sparsity', 0.2, 'plotON', 0);                                                   % create influence graph
dyn_sys.initial         = 'Unif_0_10';
dyn_sys.dt              = 1e-3;
dyn_sys.T               = dyn_sys.dt*(dyn_sys.L-1);
dyn_sys.obs_std         = 1e-3;  %  set 0 to test consistency                                       % observation noise
dyn_sys                 = system_settings( dyn_sys );                                                                           % settings of the IPS and graph and its integrator


kernel_type             =  'multitype';  %'multi_type';
basis_choice            = 'spline_long_short';
learning_set            = learning_settings( kernel_type, dyn_sys, basis_choice);    

dyn_sys.n               = learning_set.n; 
dyn_sys.phi_kernel_choices = learning_set.phi_kernel_choices;
dyn_sys.kernel_idx = learning_set.kernel_idx;
%% Demo (1 traj) generate particles
% Sanity check.
% If the kernel is wild, it will generate divergence particle trajectory.

progressON  = true;
plotON      = false;
RIPcompON   = false;

if  ~plotON
    xpath = graph_forward_model_multitype( dyn_sys, dyn_sys.X0, progressON );
elseif plotON ==1 % Plot trajectories and the graph
    I1 = dyn_sys;
    I1.L = 1000;
    I1.dt = 0.02;
    xpath = graph_forward_model_multitype(I1, I1.X0, progressON);

    if plotON==1 &&  I1.d == 2
        graph_plot_motion(xpath, I1, plotON);
        plot_graph(I1.A);
        % figure;fplot(I1.phi_kernel, [0, 5]); title('True interaction kernel')
    end

    if sum(abs(xpath(:, :, end)), 'all') > 1e10
        disp('trajectory is divergent. Kernel is bad'); return;
    end
end

%% Generate training and testing paths
% M_test                      = 2;
% Z_true                      = get_Z_from_E_c(dyn_sys.A, learning_set.c);     % Z is the product of E and c
fprintf('\nGenerating trajectories, M = %i...',M);tic
trainingPathsObj            = get_paths_multitype( dyn_sys, M,        'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
% testingPathsObj             = get_paths_multitype( dyn_sys, M_test,   'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
fprintf('done (%.2f sec).',toc);

% update the exploration measure rho and the basis matrix using paths
learning_set                    = update_dict_mat( learning_set, trainingPathsObj.paths );                                          % compute the exploration measure rho and the basis matrix using all paths
% learning_set.kernel_norm_seq    = zeros(learning_set.num_kernel_choices, 1);
% for i = 1:learning_set.num_kernel_choices
%     learning_set.kernel_norm_seq(i) = sqrt( learning_set.c_choices{i}'*learning_set.dict_mat*learning_set.c_choices{i} );                                     % MM: this is the estimated kernel norm, assuming the kernel is in the hypothesis space
% end

%%

%% Learning using Three ways alternating least squares
K_means_on_v = 1;
estALS_3_Kmeans     = learn_kernel_graph_ALS_multitype( trainingPathsObj.paths, dyn_sys, learning_set, 'niter', 10, 'normalizeON', 1, ...
    'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_u', 1e-3, ...
    'stop_thres_relDelta_v', 1e-3, 'plotON', 1, 'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
    'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] ,'K_means_on_v', K_means_on_v);% 1e-2*norm(I.A,'fro') );

%%
estALS_3_Kmeans.coef_err = zeros(estALS_3_Kmeans.ALS_n_iter, 1);
for i = 1:estALS_3_Kmeans.ALS_n_iter
    estALS_3_Kmeans.coef_err(i) = kernel_multi_err(estALS_3_Kmeans.stats.coef_mat(:, :, i), learning_set.coef_mat,  learning_set);
end
% estALS_3_Kmeans.coef_err = squeeze(sqrt(sum((estALS_3_Kmeans.stats.coef_mat - learning_set.coef_mat).^2, [1,2])));
estALS_3_Kmeans.Ahat_err = squeeze(sqrt(sum((estALS_3_Kmeans.stats.A_hat - dyn_sys.A).^2, [1,2])));


%% Learning using Three ways alternating least squares
K_means_on_v = 0;
estALS_3     = learn_kernel_graph_ALS_multitype( trainingPathsObj.paths, dyn_sys, learning_set, 'niter', 10, 'normalizeON', 1, ...
    'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_u', 1e-3, ...
    'stop_thres_relDelta_v', 1e-3, 'plotON', 1, 'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
    'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] ,'K_means_on_v', K_means_on_v);% 1e-2*norm(I.A,'fro') );

%%
estALS_3.coef_err = zeros(estALS_3.ALS_n_iter, 1);
for i = 1:estALS_3.ALS_n_iter
    estALS_3.coef_err(i) = kernel_multi_err(estALS_3.stats.coef_mat(:, :, i), learning_set.coef_mat,  learning_set);
end
% estALS_3.coef_err = squeeze(sqrt(sum((estALS_3.stats.coef_mat - learning_set.coef_mat).^2, [1,2])));
estALS_3.Ahat_err = squeeze(sqrt(sum((estALS_3.stats.A_hat - dyn_sys.A).^2, [1,2])));


%% 



%%
fprintf('True coefficient matrix')
learning_set.coef_mat

fprintf('Est coefficient matrix, No K means')
estALS_3.u*estALS_3.v'

fprintf('Est coefficient matrix, Use K means')
estALS_3_Kmeans.u*estALS_3_Kmeans.v'


estALS_3_Kmeans.sorted_id = sort_id(estALS_3_Kmeans.kernel_idx)';
estALS_3.sorted_id = sort_id(estALS_3.kernel_idx)';
learning_set.sorted_id = sort_id(learning_set.kernel_idx)';

fprintf('Kernel Classification difference, NO Kmeans')
sum(abs(estALS_3.sorted_id - learning_set.sorted_id))

fprintf('Kernel Classification difference, Kmeans')
sum(abs(estALS_3_Kmeans.sorted_id - learning_set.sorted_id))

%% Plot Q kernels based on the following rules

Q = learning_set.num_kernel_choices;
true_multi_kernel = get_multi_type_kernel_from_id(learning_set.coef_mat, learning_set.sorted_id, Q, learning_set.dict);
ALS_Kmeans_multi_kernel = get_multi_type_kernel_from_id(estALS_3_Kmeans.coef_mat, estALS_3_Kmeans.sorted_id, Q, learning_set.dict);
ALS_No_Kmeans_multi_kernel = get_multi_type_kernel_from_id(estALS_3.coef_mat, estALS_3.sorted_id, Q, learning_set.dict);






c       = colororder;
blue    = c(1, :);
red     = c(2, :);
yellow   = c(3, :);






figure;
tiledlayout(1, Q+1, 'TileSpacing', 'compact', 'padding', 'compact')

%%%%%%%%%%%%%%%%%%%%%%%%%   Error decay with iteration   %%%%%%%%%%%%%%%%%%
nexttile; hold on;grid on;
plot(0:estALS_3.ALS_n_iter-1, log10(estALS_3_Kmeans.coef_err), 'color', blue, 'LineStyle','-','LineWidth',4, ...
    'DisplayName','Using K means, kernel error');
plot(0:estALS_3.ALS_n_iter-1, log10(estALS_3_Kmeans.Ahat_err), 'color', red, 'LineStyle','-','LineWidth',4, ...
    'DisplayName','Using K means, $\textbf{a}$ error');
plot(0:estALS_3.ALS_n_iter-1, log10(estALS_3.coef_err), 'color', blue, 'LineStyle',':','LineWidth',4, ...
    'DisplayName','No K means, kernel error');
plot(0:estALS_3.ALS_n_iter-1, log10(estALS_3.Ahat_err), 'color', red, 'LineStyle',':','LineWidth',4, ...
    'DisplayName','No K means, $\textbf{a}$ error');
title('Error decay with iterations')
legend('Interpreter','Latex')
xlim([0, estALS_3.ALS_n_iter-1])
ylim([-2, 2])
xlabel('Iterations')
ylabel('log_{10} Error')
xticks(0:estALS_3.ALS_n_iter-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



right_threshold = 5;

for i = 1:Q
    nexttile; hold on;grid on;
    % L = length(learning_set.rho)*learning_set.rho_dx;
    % Num = 
    rgrid = learning_set.rho_dx:learning_set.rho_dx:right_threshold;
    L = length(rgrid);


    % nexttile; hold on;grid on;
    % L = length(learning_set.rho)*learning_set.rho_dx;
    % rgrid = learning_set.rho_dx:learning_set.rho_dx:L;

    yyaxis right
    plt_4 = area(rgrid, learning_set.rho(1:L),'FaceAlpha',0.2, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');
    ylabel('\rho')

    yyaxis left
    plt_1 = plot(rgrid, true_multi_kernel{i}(rgrid),'-', 'linewidth',4,'DisplayName','True', 'color', blue);
    plt_2 = plot(rgrid, ALS_Kmeans_multi_kernel{i}(rgrid),'-o','linewidth', 2, 'MarkerSize', 4,'DisplayName','ALS with K means', 'color', red);
    plt_3 = plot(rgrid, ALS_No_Kmeans_multi_kernel{i}(rgrid),'-*','linewidth', 2, 'MarkerSize', 2,'DisplayName','ALS without K means', 'color', yellow);


    xlabel('r')
    ylabel(['\phi_' num2str(i), '(r)'])
    if i == 1;legend([plt_1, plt_2,plt_3,plt_4], 'Location', 'southeast');end
    
    xlim([learning_set.rho_dx, 5])
    title(['Kernel Type ', num2str(i)])
    ax = gca;
    ax.YAxis(2).Color = [.7 .7 .7];
    ax.YAxis(1).Color = [.1 .1 .1];
    ax.XAxis.Color = [.1 .1 .1];
    

    % set(gca, 'Children', flipud(get(gca, 'Children')) )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'Position',[100 100 1300 300])
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

figname = [dyn_sys.PAPER_FIG_DIR, '/multitype_error_decay_kernel_est.pdf'];
saveas(gcf, figname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
