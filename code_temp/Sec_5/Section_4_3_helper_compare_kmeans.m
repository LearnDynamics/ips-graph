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
% saveas(gcf, figname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%