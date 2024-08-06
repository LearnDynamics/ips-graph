% This is an example of Kuramoto model.

clc
close all
clear all
addpaths;
rng(1)

%% specify X0 and A for generating interesting gif
X0 = [-1.5071;
        1.6298;
        0.7083;
       -0.7223;
        0.9200;
       -0.1001;
        0.2441;
        0.2873;
       -1.3940;
       -0.5521];


A =   [       0    0.3517         0         0    0.5032    0.4639         0         0         0         0;
        0.5754         0    0.6167    0.2981         0         0         0         0         0    0.6547;
             0         0         0         0         0         0    0.1980         0    0.4234         0;
        0.0768         0         0         0         0         0         0         0         0         0;
             0         0         0         0         0    0.5278    0.3771    0.8198         0    0.5280;
             0    0.7266    0.7867    0.4645         0         0    0.9048         0         0         0;
        0.8143         0         0         0         0         0         0         0    0.8902         0;
             0         0         0         0    0.1332         0         0         0         0    0.5409;
             0         0         0         0         0    0.7115         0    0.5242         0         0;
             0    0.5902    0.0274    0.8339    0.8538         0         0    0.2303    0.1680         0];


% all_M = [2, 4, 8, 16, 32, 64, 128, 256];
all_M = [8, 64, 512];
test_num = length(all_M);
batch_num = 20;

batch_num = 3;

%% load settings
% Observations and training data sizes
M                       = 1000;                                                 % number of trajectories
dyn_sys.L               = 100;      % number observations equispaced in time [0,T]

% Settings for dynamics and learning
dyn_sys.N               = 10;                                                                                                   % number of agents
dyn_sys.d               = 1;                                                                                                    % dim of state vectors
dyn_sys.viscosity       = 1e-4;  %  set 0 to test consistency                                       % stochastic force; viscosity (forcing noise in the dynamics)
% dyn_sys.A               = set_graph(dyn_sys.N, 'sparsity', 0.4, 'plotON', 0);                                                   % create influence graph
dyn_sys.A               = A;
dyn_sys.initial         = 'Unif_-2_2';
dyn_sys.dt              = 1e-3;
dyn_sys.T               = dyn_sys.dt*(dyn_sys.L-1);
dyn_sys.obs_std         = 1e-3;  %  set 0 to test consistency                                       % observation noise
dyn_sys                 = system_settings( dyn_sys );

kernel_type             = 'Kuramoto'; % 'randomSmoothFourierWithDecay'; % 6; %                  % see options in learning_settings.m
% basis_choice            = 'trigonometric';
basis_choice            = 'inHypo_Space';

learning_setup            = learning_settings( kernel_type, dyn_sys, basis_choice);
learning_setup.Z_true     = get_Z_from_E_c( dyn_sys.A, learning_setup.c );       % Z is the product of A and c                                                      % Z is the product of E and c

% dyn_sys.phi_kernel      = learning_setup.phi_kernel_cheb;                           % MM: logical inconsistency here and in the following because of this; learning_set.phi_kernel not used in any of what follows
dyn_sys.phi_kernel      = learning_setup.phi_kernel;                                % FL: learning_set.phi_kernel serves as the ground truth. --- It is a parametric inference. The chev is towards non-parametric inference
dyn_sys.n               = learning_setup.n;

true_para.A             = dyn_sys.A;
true_para.phi_kernel    = dyn_sys.phi_kernel;
% true_path           = generate_pred_path(true_para,  dyn_sys, X0);

%% Plot traj, graph and gif
plotON = 1;

if plotON ==1       % Plot trajectories and the graph
    dyn_sys_temp = dyn_sys;
    dyn_sys_temp.dt = 0.01;
    dyn_sys_temp.L = 10000;

    dyn_sys_temp.obs_std = 0;
    dyn_sys_temp.viscosity = 1e-4;

    dyn_sys_temp.tgrid = (1:dyn_sys_temp.L)*dyn_sys_temp.dt;
    true_path_temp               = generate_pred_path(true_para,  dyn_sys_temp, X0);


    %% Traj and graph
    figure;
    tiledlayout(1,2,'TileSpacing','tight','Padding','compact');
    nexttile()
    temp_tgrid = dyn_sys_temp.tgrid;
    plot(temp_tgrid, squeeze(true_path_temp)', 'linewidth', 2)
    text(60, 50, '6, 7, 9')
    text(60, 15, '2, 4, 5')
    text(60, -5, '1, 3, 8, 10')
    xlabel('Time t')
    ylabel('angel X')
    title('trajectories')

    nexttile();
    colormap(fliplr(gray(15)')')
    imagesc(A')
    colorbar
    xticks(1:10)
    yticks(1:10)
    title('Graph adjacency matrix a')

    set(gcf,'Position',[100 100 600 300])
    set(findall(gcf,'-property','FontSize'),'FontSize',13)
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    figname = [dyn_sys.PAPER_FIG_DIR, '/kuramoto_traj_graph_1d.pdf'];
    saveas(gcf, figname);

    %% Gif of particle spining on a circle
    GIF_on = 1;
    if GIF_on;graph_plot_motion_kuramoto(true_path_temp, dyn_sys_temp, 1);end
end

%% Within Hypothesis space

% rng(19)
rng(30)


basis_choice            = 'inHypo_Space';
learning_setup          = learning_settings( kernel_type, dyn_sys, basis_choice);
dyn_sys.phi_kernel      = learning_setup.phi_kernel;                                % FL: learning_set.phi_kernel serves as the ground truth. --- It is a parametric inference. The chev is towards non-parametric inference
dyn_sys.n               = learning_setup.n;

K_in = cell(test_num, batch_num);

for i = 1:test_num
    for j = 1:batch_num
        K_in{i, j} = Kuramoto_main(dyn_sys, learning_setup, all_M(i), 'RKHS_plain');
    end
end


% Not in Hypothesis space
basis_choice            = 'trigonometric';
learning_setup          = learning_settings( kernel_type, dyn_sys, basis_choice);
dyn_sys.phi_kernel      = learning_setup.phi_kernel;                                % FL: learning_set.phi_kernel serves as the ground truth. --- It is a parametric inference. The chev is towards non-parametric inference
dyn_sys.n               = learning_setup.n;

K_out = cell(test_num, batch_num);

for i = 1:test_num
    for j = 1:batch_num
        K_out{i, j} = Kuramoto_main(dyn_sys, learning_setup, all_M(i), 'RKHS_plain');
    end
end

% Plot est kernels
% Kuramoto_plot_est_kernels(K_in, learning_setup, all_M)
% Kuramoto_plot_est_kernels(K_out, learning_setup, all_M)

%%
phi_est.in_H.ALS            = cell(test_num, 1);
phi_est.in_H.ORALS          = cell(test_num, 1);
phi_est.in_H.all_rgrid      = cell(test_num, 1);

phi_est.out_H.ALS           = cell(test_num, 1);
phi_est.out_H.ORALS         = cell(test_num, 1);
phi_est.out_H.all_rgrid     = cell(test_num, 1);

for i = 1:test_num
    phi_est.in_H.all_rgrid{i}   = K_in{i, 1}.learning_setup.rho_bin_edges(1:end);
    phi_est.in_H.ALS{i}         = zeros(batch_num, length(phi_est.in_H.all_rgrid{i}));
    phi_est.in_H.ORALS{i}        = zeros(batch_num, length(phi_est.in_H.all_rgrid{i}));

    for j = 1:batch_num
        phi_est.in_H.ALS{i}(j, :)   = K_in{i, j}.estALS.phi_kernel(phi_est.in_H.all_rgrid{i});
        phi_est.in_H.ORALS{i}(j, :) = K_in{i, j}.estORALS.phi_kernel(phi_est.in_H.all_rgrid{i});
    end

    phi_est.out_H.all_rgrid{i}   = K_out{i, 1}.learning_setup.rho_bin_edges(1:end);
    phi_est.out_H.ALS{i}         = zeros(batch_num, length(phi_est.out_H.all_rgrid{i}));
    phi_est.out_H.ALS{i}        = zeros(batch_num, length(phi_est.out_H.all_rgrid{i}));

    for j = 1:batch_num
        phi_est.out_H.ALS{i}(j, :)   = K_out{i, j}.estALS.phi_kernel(phi_est.out_H.all_rgrid{i});
        phi_est.out_H.ORALS{i}(j, :) = K_out{i, j}.estORALS.phi_kernel(phi_est.out_H.all_rgrid{i});
    end
end


%%

c       = colororder;
blue    = c(1, :);
red     = c(2, :);
yellow   = c(3, :);


% figure;
figure;tiledlayout(2, 4, 'TileSpacing', 'tight', 'Padding','compact');

%%%%%%%%%%%%%%%%%%%%% Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile();
colormap(fliplr(gray(15)')')
imagesc(A')
c = colorbar;
c.Location = 'westoutside';
xticks(1:10)
yticks(1:10)
title('Adjacency matrix $\mathbf{a}$', 'Interpreter','latex')



%%%%%%%%%%%%%%%%%%%%% Kernel est in Hypo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
in_ALS_stats = cell(test_num, 1);
in_ORALS_stats = cell(test_num, 1);
for i = 1:test_num
    nexttile;hold on;grid on;
    rho     = K_in{i}.learning_setup.rho;
    dx      = K_in{i}.learning_setup.rho_dx;


    in_ALS_stats{i}.mean = mean(phi_est.in_H.ALS{i})';
    in_ALS_stats{i}.std = std(phi_est.in_H.ALS{i})';

    % Some ORALS gives NAN
    ind = isnan(sum(phi_est.in_H.ORALS{i}, 2));
    phi_est.in_H.ORALS{i}(ind, :) = [];

    in_ORALS_stats{i}.mean = mean(phi_est.in_H.ORALS{i})';
    in_ORALS_stats{i}.std = std(phi_est.in_H.ORALS{i})';


    plt_1 = plot(phi_est.in_H.all_rgrid{i}, learning_setup.phi_kernel(phi_est.in_H.all_rgrid{i}),'-', 'linewidth',4,'DisplayName','True');
    plt_2 = plot(phi_est.in_H.all_rgrid{i}, in_ALS_stats{i}.mean, '-', 'color', red, 'linewidth',2,'DisplayName','ALS');
    plt_3 = plot(phi_est.in_H.all_rgrid{i}, in_ORALS_stats{i}.mean,'-', 'color', yellow, 'linewidth',2,'DisplayName','ORALS');
    plt_4 = area(phi_est.in_H.all_rgrid{i}(2:end), rho,'FaceAlpha',0.4, 'EdgeAlpha',0.1, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');

    % plot(phi_est.in_H.all_rgrid{i}, K_in{i}.estALS.phi_kernel(phi_est.in_H.all_rgrid{i}), '.', 'color', red, 'linewidth',2,'DisplayName','ALS');
    % plot(phi_est.in_H.all_rgrid{i}, K_in{i}.estORALS.phi_kernel(phi_est.in_H.all_rgrid{i}), '.', 'color', yellow, 'linewidth',2,'DisplayName','ALS');

    fill([phi_est.in_H.all_rgrid{i}, fliplr(phi_est.in_H.all_rgrid{i})], [in_ORALS_stats{i}.mean - in_ORALS_stats{i}.std; flip(in_ORALS_stats{i}.mean + in_ORALS_stats{i}.std)], yellow, 'FaceAlpha',0.3, 'EdgeAlpha',0);
    fill([phi_est.in_H.all_rgrid{i}, fliplr(phi_est.in_H.all_rgrid{i})], [in_ALS_stats{i}.mean - in_ALS_stats{i}.std; flip(in_ALS_stats{i}.mean + in_ALS_stats{i}.std)], red, 'FaceAlpha',0.3, 'EdgeAlpha',0);

    xlim([dx, phi_est.in_H.all_rgrid{i}(end)])
    ylim([-2.5, 2.5])
    title(['$M = ', num2str(all_M(i)), ', \mathcal{H}_\phi$'], 'Interpreter','latex')

    if i == 3
        legend([plt_1, plt_2,plt_3,plt_4], 'Location', 'southwest');
        ylabel('$\phi(r)$', 'Interpreter','latex')
    end
end


    
%%%%%%%%%%%%%%%%%%%%% Traj %%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile;
temp_tgrid = dyn_sys_temp.tgrid;
    plot(temp_tgrid, squeeze(true_path_temp)', 'linewidth', 2)
    text(60, 50, '6, 7, 9')
    text(60, 15, '2, 4, 5')
    text(60, -5, '1, 3, 8, 10')
    xlabel('Time $t$', 'Interpreter','latex')
    ylabel('Angel $X$', 'Interpreter','latex')
    title('Trajectories', 'Interpreter','latex')



%%%%%%%%%%%%%%%%%%%%% Kernel est outside Hypo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_ALS_stats = cell(test_num, 1);
out_ORALS_stats = cell(test_num, 1);
for i = 1:test_num
    rho     = K_out{i}.learning_setup.rho;
    dx      = K_out{i}.learning_setup.rho_dx;

    nexttile;hold on;grid on;
    out_ALS_stats{i}.mean = mean(phi_est.out_H.ALS{i})';
    out_ALS_stats{i}.std = std(phi_est.out_H.ALS{i})';

    % Some ORALS gives NAN
    ind = isnan(sum(phi_est.out_H.ORALS{i}, 2));
    phi_est.out_H.ORALS{i}(ind, :) = [];

    out_ORALS_stats{i}.mean = mean(phi_est.out_H.ORALS{i})';
    out_ORALS_stats{i}.std = std(phi_est.out_H.ORALS{i})';


    plt_1 = plot(phi_est.out_H.all_rgrid{i}, learning_setup.phi_kernel(phi_est.out_H.all_rgrid{i}),'-', 'linewidth',4,'DisplayName','True');
    plt_2 = plot(phi_est.out_H.all_rgrid{i}, out_ALS_stats{i}.mean, '-', 'color', red, 'linewidth',2,'DisplayName','ALS');
    plt_3 = plot(phi_est.out_H.all_rgrid{i}, out_ORALS_stats{i}.mean,'-', 'color', yellow, 'linewidth',2,'DisplayName','ORALS');
    plt_4 = area(phi_est.out_H.all_rgrid{i}(2:end), rho,'FaceAlpha',0.4, 'EdgeAlpha',0.1, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');

    % plot(phi_est.out_H.all_rgrid{i}, K_out{i}.estALS.phi_kernel(phi_est.out_H.all_rgrid{i}), '.', 'color', red, 'linewidth',2,'DisplayName','ALS');
    % plot(phi_est.out_H.all_rgrid{i}, K_out{i}.estORALS.phi_kernel(phi_est.out_H.all_rgrid{i}), '.', 'color', yellow, 'linewidth',2,'DisplayName','ALS');


    fill([phi_est.out_H.all_rgrid{i}, fliplr(phi_est.out_H.all_rgrid{i})], [out_ORALS_stats{i}.mean - out_ORALS_stats{i}.std; flip(out_ORALS_stats{i}.mean + out_ORALS_stats{i}.std)], yellow, 'FaceAlpha',0.3, 'EdgeAlpha',0);
    fill([phi_est.out_H.all_rgrid{i}, fliplr(phi_est.out_H.all_rgrid{i})], [out_ALS_stats{i}.mean - out_ALS_stats{i}.std; flip(out_ALS_stats{i}.mean + out_ALS_stats{i}.std)], red, 'FaceAlpha',0.3, 'EdgeAlpha',0);

    xlim([dx, phi_est.out_H.all_rgrid{i}(end)])
    ylim([-2.5, 2.5])
    title(['$M = ', num2str(all_M(i)), ', \mathcal{H}$'], 'Interpreter','latex')

    if i == 3
        legend([plt_1, plt_2,plt_3,plt_4], 'Location', 'southwest');
        ylabel('$\phi(r)$', 'Interpreter','latex')
    end
end







%%%%%%%%%%%%%%%%%%%%%%% Set font size and paperposition %%%%%%%%%%%%%%%%%%%
set(gcf,'Position',[100 100 1050 500])
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figname = [dyn_sys.PAPER_FIG_DIR, '/Kuramoto_all.pdf'];
saveas(gcf, figname);


%% Giant figure
% 
% figure;tiledlayout(2, 4, 'TileSpacing', 'tight', 'Padding','compact');
% 
% %%%%%%%%%%%%%%%%%%%%% Graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nexttile();
% colormap(fliplr(gray(15)')')
% imagesc(A')
% c = colorbar;
% c.Location = 'westoutside';
% xticks(1:10)
% yticks(1:10)
% title('Adjacency matrix $\mathbf{a}$', 'Interpreter','latex')
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%% Kernel est in Hypo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:test_num
% 
%     nexttile;hold on;grid on;
% 
%     rho     = K_in{i}.learning_setup.rho;
%     dx      = K_in{i}.learning_setup.rho_dx;
%     rgrid   = K_in{i}.learning_setup.rho_bin_edges(2:end);
% 
%     plt_1 = plot(rgrid, learning_setup.phi_kernel(rgrid),'-', 'linewidth',4,'DisplayName','True');
%     plt_2 = plot(rgrid, K_in{i}.estALS.phi_kernel(rgrid),'-o','linewidth',2,'DisplayName','ALS');
%     plt_3 = plot(rgrid, K_in{i}.estORALS.phi_kernel(rgrid),'-*','linewidth',2,'DisplayName','ORALS');
%     plt_4 = area(rgrid, rho,'FaceAlpha',0.4, 'EdgeAlpha',0.1, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');
% 
%     % xlabel('r')
%     xlim([dx, rgrid(end)])
%     ylim([-2.5, 2.5])
%     title(['$M = ', num2str(all_M(i)), ', \mathcal{H}_\phi$'], 'Interpreter','latex')
% 
%     if i == 1
%         legend([plt_1, plt_2,plt_3,plt_4], 'Location', 'southwest');
%         ylabel('$\phi(r)$', 'Interpreter','latex')
%     end
% end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%% Traj %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nexttile;
% temp_tgrid = dyn_sys_temp.tgrid;
%     plot(temp_tgrid, squeeze(true_path_temp)', 'linewidth', 2)
%     text(60, 50, '6, 7, 9')
%     text(60, 15, '2, 4, 5')
%     text(60, -5, '1, 3, 8, 10')
%     xlabel('Time $t$', 'Interpreter','latex')
%     ylabel('Angel $X$', 'Interpreter','latex')
%     title('Trajectories', 'Interpreter','latex')
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%% Kernel est outside Hypo %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:test_num
% 
%     nexttile;hold on;grid on;
% 
%     rho     = K_out{i}.learning_setup.rho;
%     dx      = K_out{i}.learning_setup.rho_dx;
%     rgrid   = K_out{i}.learning_setup.rho_bin_edges(2:end);
% 
%     plt_1 = plot(rgrid, learning_setup.phi_kernel(rgrid),'-', 'linewidth',4,'DisplayName','True');
%     plt_2 = plot(rgrid, K_out{i}.estALS.phi_kernel(rgrid),'-o','linewidth',2,'DisplayName','ALS');
%     plt_3 = plot(rgrid, K_out{i}.estORALS.phi_kernel(rgrid),'-*','linewidth',2,'DisplayName','ORALS');
%     plt_4 = area(rgrid, rho,'FaceAlpha',0.4, 'EdgeAlpha',0.1, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');
% 
%     xlabel('$r$', 'Interpreter','latex')
%     xlim([dx, rgrid(end)])
%     ylim([-2.5, 2.5])
%     title(['$M = ', num2str(all_M(i)), ', \mathcal{H}$'], 'Interpreter','latex')
% 
%     if i == 1
%         legend([plt_1, plt_2,plt_3,plt_4], 'Location', 'southwest');
%         ylabel('$\phi(r)$', 'Interpreter','latex')
%     end
% end
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%% Set font size and paperposition %%%%%%%%%%%%%%%%%%%
% set(gcf,'Position',[100 100 1050 500])
% set(findall(gcf,'-property','FontSize'),'FontSize',15)
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% figname = [dyn_sys.PAPER_FIG_DIR, '/Kuramoto_all.pdf'];
% % saveas(gcf, figname);

