function plot_graph(E)
figure;
tiledlayout(1,2, 'TileSpacing','compact', 'padding','compact');
nexttile;
% subplot(121);
G = digraph(E,'omitselfloops');
plot(G,'Layout','force', 'EdgeLabel', round(G.Edges.Weight,2),'linewidth',1,'markersize',10)
% set(gca, 'Position',[0.1 0.29 0.4 0.45]);
% title('graph')
% subplot(122); 
nexttile;
spy(E);
imagesc(E');
colormap(fliplr(gray(15)')')
colorbar
xlabel('')
% set(gca, 'Position',[0.6 0.29 0.4 0.45]);
set(gcf, 'Position', [100, 100, 500, 240]);
end