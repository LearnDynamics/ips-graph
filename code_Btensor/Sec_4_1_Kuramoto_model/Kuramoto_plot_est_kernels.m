function Kuramoto_plot_est_kernels(K, learning_setup, all_M)

test_num = length(all_M);

figure;tiledlayout(1, test_num, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:test_num
    
    nexttile;hold on;grid on;

    rho     = K{i}.learning_setup.rho;
    dx      = K{i}.learning_setup.rho_dx;
    rgrid   = K{i}.learning_setup.rho_bin_edges(2:end);

    plt_1 = plot(rgrid, learning_setup.phi_kernel(rgrid),'-', 'linewidth',4,'DisplayName','True');
    plt_2 = plot(rgrid, K{i}.estALS.phi_kernel(rgrid),'-o','linewidth',2,'DisplayName','ALS');
    plt_3 = plot(rgrid, K{i}.estORALS.phi_kernel(rgrid),'-*','linewidth',2,'DisplayName','ORALS');
    plt_4 = area(rgrid, rho,'FaceAlpha',0.4, 'EdgeAlpha',0.1, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');

    xlabel('r')
    xlim([dx, rgrid(end)])
    ylim([-2.5, 2.5])
    title(['M = ', num2str(all_M(i))])

    if i == 1
        legend([plt_1, plt_2,plt_3,plt_4], 'Location', 'southwest');
        ylabel('\phi(r)', 'Interpreter','latex')
    end

end

%%%%%%%%%%%%%%%%%%%%%%% Set font size and paperposition %%%%%%%%%%%%%%%%%%%
set(gcf,'Position',[100 100 1000 350])
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% figname = [dyn_sys.PAPER_FIG_DIR, '/Kuramoto_all.pdf'];
% saveas(gcf, figname);




end