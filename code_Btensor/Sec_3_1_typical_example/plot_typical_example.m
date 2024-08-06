%% plot graph
figure;
tiledlayout(1, 4, 'TileSpacing', 'normal', 'padding', 'normal', 'InnerPosition', [0.04, 0.11, 0.91, 0.75])
% This command is very efficient in setting a tight figutre
% You can customize the margin in a tiledlayout.

nexttile;
G = digraph(dyn_sys.A,'omitselfloops');
plot(G,'Layout','force', 'EdgeLabel', round(G.Edges.Weight,2),'linewidth',1,'markersize',10)
% title('True Graph', 'Interpreter', 'latex')
title('True Graph')
nexttile;
imagesc(dyn_sys.A')
colormap(fliplr(gray(15)')')
% colormap default
% colorbar('Ticks',colorbarticks)
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
xticks(1:6)
yticks(1:6)
colorbar
title('True a')
% title('Adjacency Matrix $$\mathbf{a}$$', 'Interpreter', 'latex')


nexttile;
imagesc(abs(estALS.A - dyn_sys.A)');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
barH = colorbar; barH.Ruler.Exponent = -3;
title('|a - a_{ALS}|')
% title('$$|\mathbf{a} - \mathbf{a}_{ALS}|$$', 'Interpreter', 'latex')

nexttile;
imagesc(abs(estORALS.A - dyn_sys.A)');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
barH = colorbar; barH.Ruler.Exponent = -3;
% title('$$|\mathbf{a} - \mathbf{a}_{ORALS}|$$', 'Interpreter', 'latex')
title('|a - a_{ORALS}|')
% set_positionFontsAll;
set(gcf,'Position',[100 100 1300 320])
set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

figname = [dyn_sys.PAPER_FIG_DIR, '/graph_true_est_',options, '.pdf'];
saveas(gcf, figname);

% figname = [SLIDE_FIG_DIR, '/graph_true_est_spiral.pdf'];
% saveas(gcf, figname);

%% kernel plot



figure;
tiledlayout(2, 4, 'TileSpacing', 'compact', 'padding', 'compact')

%%%%%%%%%%%%%%%%%%%%%%%% true traj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(1, [2, 1]);
hold on; view(3);grid on
for i = 1:dyn_sys.N
    plot3(squeeze(true_path(i, 1, :)), squeeze(true_path(i, 2, :)), dyn_sys.tgrid, 'linewidth', 3)
    if i == 4
        text(true_path(i, 1, 1)-0.09, true_path(i, 2, 1)+0.1, 0, ['X_', num2str(i)],'HorizontalAlignment','left','FontSize',10);
    elseif i == 3
        text(true_path(i, 1, 1)-0.1, true_path(i, 2, 1), 0, ['X_', num2str(i)],'HorizontalAlignment','left','FontSize',10);
    else
        text(true_path(i, 1, 1)+0.05, true_path(i, 2, 1), 0, ['X_', num2str(i)],'HorizontalAlignment','left','FontSize',10);
    end
end
xlabel('Coordinate 1', 'Rotation',27)
ylabel('Coordinate 2', 'Rotation',-40)
zlabel('Time')
title('True trajectory')

phi_kernel_ALS      = estALS.phi_kernel;
phi_kernel_ORALS    = estORALS.phi_kernel;


% figure;
% tiledlayout(1, 2, 'TileSpacing', 'normal', 'padding', 'compact')


%%%%%%%%%%%%%%%% kernel estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(2, [2, 1])
hold on;grid on;
L = length(learning_setup.rho)*learning_setup.rho_dx;
rgrid = learning_setup.rho_dx:learning_setup.rho_dx:L;
plt_1 = plot(rgrid, dyn_sys.phi_kernel(rgrid),'-', 'linewidth',4,'DisplayName','True');
plt_2 = plot(rgrid, phi_kernel_ALS(rgrid),'-o','linewidth', 2, 'MarkerSize', 4,'DisplayName','ALS');
plt_3 = plot(rgrid, phi_kernel_ORALS(rgrid),'-*','linewidth', 2, 'MarkerSize', 2,'DisplayName','ORALS');
ylim([-170, 22])
% set(gca, 'YScale', 'log')
rectangle('position',[0.8, -10, 0.5, 22], 'LineWidth', 1)
xlabel('r')
ylabel('\phi(r)')
yyaxis right
plt_4 = area(rgrid, learning_setup.rho,'FaceAlpha',0.4, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');
legend([plt_1, plt_2,plt_3,plt_4], 'Location', 'east')
xlim([learning_setup.rho_dx, L])
title('Kernel estimation')
ax = gca;
ax.YAxis(2).Color = [.7 .7 .7];
ylabel('\rho')




%%%%%%%%%%%%%%%% kernel estiamtion zoom in $%%%%%%%%%%%%%%%%
nexttile(3, [2, 1])
hold on;grid on;
L = length(learning_setup.rho)*learning_setup.rho_dx;
rgrid = learning_setup.rho_dx:learning_setup.rho_dx:L;
plt_1 = plot(rgrid, dyn_sys.phi_kernel(rgrid),'-', 'linewidth',4,'DisplayName','True');
plt_2 = plot(rgrid, phi_kernel_ALS(rgrid),'-o','linewidth',1,'DisplayName','ALS');
plt_3 = plot(rgrid, phi_kernel_ORALS(rgrid),'-*','linewidth',1,'DisplayName','ORALS');
xlabel('r')
ylabel('\phi(r)')

yyaxis right
plt_4 = area(rgrid, learning_setup.rho,'FaceAlpha',0.4, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');
% legend([plt_1, plt_2,plt_3,plt_4], 'Location', 'northeast')
xlim([0.8, 1.37])
title('Zoom in the rectangle')
ax = gca;
ax.YAxis(2).Color = [.7 .7 .7];
ylabel('\rho')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







% %%%%%%%%%% estimated traj%%%%%%%%%%
nexttile(4)
id = [1,2,3,4,5,6];
color_order = colororder;
% nexttile(2);
hold on;grid on;
text_t = 1;
% for i = 1:dyn_sys.N
%     % if i == 3
%     %     text(dyn_sys.tgrid(text_t), true_path(i, 1, text_t)-0.06, ['X_', num2str(i)]);
%     % elseif i == 1
%     %     text(dyn_sys.tgrid(text_t), true_path(i, 1, text_t)+0.09, ['X_', num2str(i)]);
%     % else
%         text(dyn_sys.tgrid(text_t), true_path(i, 1, text_t)+0.06, ['X_', num2str(i)]);
%     % end
% end
p1 = plot(dyn_sys.tgrid, squeeze(true_path(id, 1, :))', 'linewidth', 4, 'color', color_order(1, :),'DisplayName','True');
p2 = plot(dyn_sys.tgrid, squeeze(pred_path_ALS(id, 1, :))','-o', 'linewidth', 2, 'color', color_order(2, :), 'MarkerSize', 2,'DisplayName','ALS');
p3 = plot(dyn_sys.tgrid, squeeze(pred_path_ORALS(id, 1, :))', '--', 'linewidth', 2, 'color', color_order(3, :), 'MarkerSize', 2,'DisplayName','ORALS');
% legend([p1(1), p2(1), p3(1)])
ylabel('Coordinate 1')
title('Estimated trajectory')
xlim([dyn_sys.tgrid(1), dyn_sys.tgrid(end)])



nexttile(8)
hold on;grid on;
% for i = 1:dyn_sys.N
%     if i == 5
%         text(dyn_sys.tgrid(text_t), true_path(i, 2, text_t)+0.136, ['X_', num2str(i)]);
%     elseif i == 4
%         text(dyn_sys.tgrid(text_t), true_path(i, 2, text_t)-0.09, ['X_', num2str(i)]);
%     else
%         text(dyn_sys.tgrid(text_t), true_path(i, 2, text_t)+0.1, ['X_', num2str(i)]);
%     end
% end
p1 = plot(dyn_sys.tgrid, squeeze(true_path(id, 2, :))', 'linewidth', 4, 'color', color_order(1, :));
p2 = plot(dyn_sys.tgrid, squeeze(pred_path_ALS(id, 2, :))', '-o', 'linewidth', 2, 'color', color_order(2, :), 'MarkerSize', 2);
p3 = plot(dyn_sys.tgrid, squeeze(pred_path_ORALS(id, 2, :))', '--', 'linewidth', 2, 'color', color_order(3, :), 'MarkerSize', 2);

xlabel('Time')
ylabel('Coordinate 2')
xlim([dyn_sys.tgrid(1), dyn_sys.tgrid(end)])

% legend([p1(1), p2(1), p3(1)], {'True', 'ORALS', 'ALS'}, 'location', 'best')
set(gcf,'Position',[100 100 1300 300])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% set_positionFontsAll;

figname = [dyn_sys.PAPER_FIG_DIR, '/traj_kernel_true_est_',options, '.pdf'];
saveas(gcf, figname);

% figname = [SLIDE_FIG_DIR, '/traj_kernel_true_est_spiral.pdf'];
% saveas(gcf, figname);





%% auxilary figure
figure;
tiledlayout(1, 4, 'TileSpacing', 'compact', 'padding', 'compact', 'InnerPosition', [0.04, 0.11, 0.91, 0.75])
% This command is very efficient in setting a tight figutre
% You can customize the margin in a tiledlayout.

nexttile;
G = digraph(dyn_sys.A,'omitselfloops');
plot(G,'Layout','force', 'EdgeLabel', round(G.Edges.Weight,2),'linewidth',1,'markersize',10)
% title('True Graph', 'Interpreter', 'latex')
title('True Graph')





nexttile;
imagesc(dyn_sys.A')
colormap(fliplr(gray(15)')')
% colormap default
% colorbar('Ticks',colorbarticks)
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
xticks(1:6)
yticks(1:6)
colorbar
title('True a')
% title('Adjacency Matrix $$\mathbf{a}$$', 'Interpreter', 'latex')


nexttile;
hold on; view(3);grid on
for i = 1:dyn_sys.N
    plot3(squeeze(true_path(i, 1, :)), squeeze(true_path(i, 2, :)), dyn_sys.tgrid, 'linewidth', 3)
    if i == 4
        text(true_path(i, 1, 1)-0.09, true_path(i, 2, 1)+0.1, 0, ['X_', num2str(i)],'HorizontalAlignment','left','FontSize',10);
    elseif i == 3
        text(true_path(i, 1, 1)-0.1, true_path(i, 2, 1), 0, ['X_', num2str(i)],'HorizontalAlignment','left','FontSize',10);
    else
        text(true_path(i, 1, 1)+0.05, true_path(i, 2, 1), 0, ['X_', num2str(i)],'HorizontalAlignment','left','FontSize',10);
    end
end
xlabel('Coordinate 1', 'Rotation',27)
ylabel('Coordinate 2', 'Rotation',-40)
zlabel('Time')
title('True trajectory')

phi_kernel_ALS      = estALS.phi_kernel;
phi_kernel_ORALS    = estORALS.phi_kernel;




nexttile;
hold on;grid on;
L = length(learning_setup.rho)*learning_setup.rho_dx;
rgrid = learning_setup.rho_dx:learning_setup.rho_dx:L;
plt_1 = plot(rgrid, dyn_sys.phi_kernel(rgrid),'-', 'linewidth',4,'DisplayName','True');
plt_2 = plot(rgrid, phi_kernel_ALS(rgrid),'-o','linewidth',1,'DisplayName','ALS');
plt_3 = plot(rgrid, phi_kernel_ORALS(rgrid),'-*','linewidth',1,'DisplayName','ORALS');
xlabel('r')
ylabel('\phi(r)')

yyaxis right
plt_4 = area(rgrid, learning_setup.rho,'FaceAlpha',0.4, 'FaceColor',[.7 .7 .7], 'LineStyle','none','DisplayName','\rho');
legend([plt_1, plt_2,plt_3,plt_4], 'Location', 'north')
xlim([0.8, 1.37])
title('Kernel estimation')
ax = gca;
ax.YAxis(2).Color = [.7 .7 .7];
ylabel('\rho')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'Position',[100 100 1300 300])
set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

figname = [dyn_sys.PAPER_FIG_DIR, '/traj_kernel_graph_short_',options, '.pdf'];
saveas(gcf, figname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% save figure to paper figure folder
% figname = [PAPER_FIG_DIR, '/traj_pred.pdf'];
% saveas(gcf, figname);