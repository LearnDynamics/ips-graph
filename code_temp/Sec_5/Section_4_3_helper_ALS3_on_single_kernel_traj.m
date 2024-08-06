
c       = colororder;
blue    = c(1, :);
red     = c(2, :);
yellow   = c(3, :);


Q = learning_set.num_kernel_choices;
true_multi_kernel = get_multi_type_kernel_from_id(learning_set.coef_mat, learning_set.sorted_id, Q, learning_set.dict);
ALS_Kmeans_multi_kernel = get_multi_type_kernel_from_id(estALS_3_single.coef_mat, estALS_3_single.sorted_id, Q, learning_set.dict);
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
% saveas(gcf, figname);