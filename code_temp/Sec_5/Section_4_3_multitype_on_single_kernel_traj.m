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



kernel_type             = 'multitype';  %'multi_type';
basis_choice            = 'fake_multitype';
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

estALS_3_Kmeans.sorted_id = sort_id(estALS_3_Kmeans.kernel_idx)';
% estALS_3.sorted_id = sort_id(estALS_3.kernel_idx)';
learning_set.sorted_id = sort_id(learning_set.kernel_idx)';

% fprintf('Kernel Classification difference, Kmeans')
% sum(abs(estALS_3.sorted_id - learning_set.sorted_id))

fprintf('Kernel Classification difference, NO Kmeans')
sum(abs(estALS_3_Kmeans.sorted_id - learning_set.sorted_id))

%% Plot Q kernels based on the following rules

c       = colororder;
blue    = c(1, :);
red     = c(2, :);
yellow   = c(3, :);


Q = learning_set.num_kernel_choices;
true_multi_kernel = get_multi_type_kernel_from_id(learning_set.coef_mat, learning_set.sorted_id, Q, learning_set.dict);
ALS_Kmeans_multi_kernel = get_multi_type_kernel_from_id(estALS_3_Kmeans.coef_mat, estALS_3_Kmeans.sorted_id, Q, learning_set.dict);
% ALS_No_Kmeans_multi_kernel = get_multi_type_kernel_from_id(estALS_3.coef_mat, estALS_3.sorted_id, Q, learning_set.dict);


figure;
tiledlayout(1, Q, 'TileSpacing', 'compact', 'padding', 'compact')

right_threshold = 5;

for i = 1:Q
    nexttile; hold on;grid on;
    % L = length(learning_set.rho)*learning_set.rho_dx;
    % Num = 
    rgrid = learning_set.rho_dx:learning_set.rho_dx:right_threshold;
    L = length(rgrid);

    yyaxis right
    plt_4 = area(rgrid, learning_set.rho(1:L),'FaceAlpha',0.2, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');
    ylabel('\rho')

    yyaxis left
    plt_1 = plot(rgrid, true_multi_kernel{1}(rgrid),'-', 'linewidth',4,'DisplayName','True', 'color', blue);
    plt_2 = plot(rgrid, ALS_Kmeans_multi_kernel{i}(rgrid),'-o','linewidth', 2, 'MarkerSize', 4,'DisplayName','ALS of 2 types', 'color', red);

    xlabel('r')
    ylabel(['\phi_' num2str(i), '(r)'])
    if i == 1;legend([plt_1, plt_2, plt_4], 'Location', 'southeast');end
    
    xlim([learning_set.rho_dx, 5])
    title(['Kernel Type ', num2str(i)])
    ax = gca;
    ax.YAxis(2).Color = [.7 .7 .7];
    ax.YAxis(1).Color = [.1 .1 .1];
    ax.XAxis.Color = [.1 .1 .1];
    
    % set(gca, 'Children', flipud(get(gca, 'Children')) )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'Position',[100 100 700 300])
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

figname = [dyn_sys.PAPER_FIG_DIR, '/multitype_single_traj.pdf'];
saveas(gcf, figname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
