% Plots of estimation errors
warning off;

bigFig;h=gcf;set(h,'NumberTitle', 'off', 'Name', 'Estimation errors');hold on;
obs_std = dyn_sys.obs_std;
plotHandles = [];
legendStrings = {};
plotxlength = length(estALS.ALS_error_graph);

plotHandles(end+1)=plot( log10( estALS.ALS_error_kernel ),'r-o','color',[0.937 0.10196 0.1765],'linewidth',2,'MarkerSize',10,'MarkerFaceColor',[0.937 0.10196 0.1765]);                              % ALS kernel error, relative
legendStrings{end+1} = 'ALS - relL2 error kernel';

plotHandles(end+1)=plot( log10( estALS.ALS_error_graph ),'-x','color',[0.8290 0.2940 0.0250],'linewidth',2,'MarkerSize',10,'MarkerFaceColor',[0.8290 0.2940 0.0250]);                                % ALS graph error, relative
legendStrings{end+1} = 'ALS - relFro error graph';

if isfield(estALS.stats,'err_relmeanL2')
    plotHandles(end+1)=plot( log10( estALS.stats.err_relmeanL2 ),'-*','color',[0.7290 0.1940 0.050],'linewidth',2,'MarkerSize',10,'MarkerFaceColor',[0.7290 0.1940 0.050]);                                % ALS graph error, relative
else
    plotHandles(end+1)=plot( log10( estALS.pathTestErr)*ones(1, plotxlength),'-*','color',[0.7290 0.1940 0.050],'linewidth',2,'MarkerSize',10,'MarkerFaceColor', 'black');   % ORALS graph error, relative
end
legendStrings{end+1} = 'ALS - relL2 error paths';

if runORALS
    plotHandles(end+1)=plot( log10( kernel_err(estORALS.c, learning_setup)*ones(1, plotxlength) ),'b-o','linewidth',2,'MarkerSize',10,'MarkerFaceColor', 'black');  % ORALS kernel error, relative
    legendStrings{end+1} = 'ORALS - relL2 error kernel';

    plotHandles(end+1)=plot( log10( graph_err(estORALS.A, dyn_sys)*ones(1, plotxlength) ),'-x','color',[0.3010 0.7450 0.9330],'linewidth',2,'MarkerSize',10,'MarkerFaceColor', 'black');   % ORALS graph error, relative
    legendStrings{end+1} = 'ORALS - relFro error graph';

    plotHandles(end+1)=plot( log10( estORALS.pathTestErr)*ones(1, plotxlength),'-*','color',[0.2010 0.6450 0.7330],'linewidth',2,'MarkerSize',10,'MarkerFaceColor', 'black');   % ORALS graph error, relative
    legendStrings{end+1} = 'ORALS - relL2 error paths';
end
plotHandles(end+1)=plot( log10( 10*dyn_sys.obs_std*sqrt(dyn_sys.N*dyn_sys.d*dyn_sys.L)/estALS.meanL2traj*ones(1, plotxlength) ),'k:','linewidth',2);                                                % std observation error
legendStrings{end+1} = '10\sigma_{obs} on traj.';

plotHandles(end+1)=plot( log10( dyn_sys.viscosity*sqrt(dyn_sys.N*dyn_sys.d*dyn_sys.L)/estALS.meanL2traj*ones(1, plotxlength) ),'k-.','linewidth',2);                                                % std observation error
legendStrings{end+1} = '\sigma_{stoch. forcing} on traj.';

legend( plotHandles,legendStrings,'fontsize',16);

ylabel('log_{10}(Error)','fontsize',16); xlabel('Iteration number','fontsize',16);
title('ALS and ORALS estimation errors','fontsize',22)
axis tight; grid on;

%% Plots of estimated interaction kernels
bigFig;h=gcf;set(h,'NumberTitle', 'off', 'Name', 'Estimation of the interaction kernel');
phi_kernel          = get_kernel_from_c(learning_setup.c, learning_setup.dict);

if runORALS,phi_kernel_ORALS    = get_kernel_from_c(estORALS.c, learning_setup.dict); end
if runALS,  phi_kernel_ALS      = get_kernel_from_c(estALS.c, learning_setup.dict);   end

r = linspace(0,learning_setup.rho_bin_edges(end),1000);

subplot(1,2,1);hold on;
y1 = plot( r,phi_kernel(r),'k','linewidth',2);
if runALS,      y2 = plot( r,phi_kernel_ALS(r),'--','linewidth',2,'color',[0.937 0.10196 0.1765]); else y2=[]; end
if runORALS,    y3 = plot( r,phi_kernel_ORALS(r),'b:','linewidth',3); else y3=[]; end
xlabel('r','FontSize',16)
ylabel('Values of the interaction kernels','FontSize',16)
yyaxis right
y4 = plot( 0.5*(learning_setup.rho_bin_edges(1:end-1)+learning_setup.rho_bin_edges(2:end)), learning_setup.rho,'-','color',[0.9,0.5,0] );
ylabel('\rho','FontSize',16)
axis tight
ylim([0,1]);
%yyaxis left;ylimcur = ylim;ylim([ylimcur(1)-abs(min(phi_kernel(r))),ylimcur(2)])
legend( [y1,y2,y3],{'\phi_{true}','\phi_{ALS}','\phi_{ORALS}'},'FontSize',16 )

subplot(1,2,2);hold on;
if runALS, y2 = plot( r,phi_kernel_ALS(r)-phi_kernel(r),'r-','linewidth',2); end
if runORALS, y3 = plot( r,phi_kernel_ORALS(r)-phi_kernel(r),'b-','linewidth',2); end
xlabel('r','FontSize',16)
ylabel('Values of the interaction kernels','FontSize',16)
yyaxis right
y4 = plot( 0.5*(learning_setup.rho_bin_edges(1:end-1)+learning_setup.rho_bin_edges(2:end)), learning_setup.rho,'-','color',[0.9,0.5,0] );
ylabel('\rho','FontSize',16)
axis tight
ylim([0,1]);
legend( [y2,y3],{'\phi_{ALS}-\phi_{true}','\phi_{ORALS}-\phi_{true}'},'FontSize',16 );

%% Plots of estimated graph adjacency matrix
bigFig;h=gcf;set(h,'NumberTitle', 'off', 'Name', 'Estimation of A');
subplot(2,3,1); imagesc(dyn_sys.A); colormap('copper');colorbar; axis equal; axis tight; c1 = caxis; title('A_{true}');
if runALS, subplot(2,3,2); imagesc(estALS.A); colorbar; axis equal; axis tight; c2 = caxis; title('A_{ALS}'); else c2 = 0; end
if runORALS, subplot(2,3,3); imagesc(estORALS.A); colorbar; axis equal; axis tight; c3 = caxis; title('A_{ORALS}'); else c3 = 0; end
c3 = [min([c1 c2 c3]), max([c1 c2 c3])]; caxis(c3); colorbar off; cb=colorbar(gca,'EastOutside');cb.Position(1)=cb.Position(1)+0.08;
subplot(2,3,1); caxis(c3); colorbar off; subplot(2,3,2); caxis(c3); colorbar off;
subplot(2,3,5); imagesc(dyn_sys.A-estALS.A); colorbar; axis equal; axis tight; title('A_{true}-A_{ALS}')
c2 = caxis;
if runORALS, subplot(2,3,6); imagesc(dyn_sys.A-estORALS.A); colorbar; axis equal; axis tight; title('A_{true}-A_{ORALS}'); end
c1 = caxis; c3 = [min([c1 c2]), max([c1 c2])]; caxis(c3); colorbar off; cb=colorbar(gca,'EastOutside');cb.Position(1)=cb.Position(1)+0.08;
subplot(2,3,5); colorbar off;

subplot(2,3,4); colorbar off; hold on;
plot(svd(full(dyn_sys.A)),'k','linewidth',2);
if runALS, plot(svd(estALS.A),'r--','linewidth',2); end
if runORALS, plot(svd(estORALS.A),'b:','linewidth',3); end
axis tight; title('SVD of adjacency matrices')
legend({'\Sigma(A)','\Sigma(A_{ALS})','\Sigma(A_{ORALS})'});

%% Plots of some trajectoris
bigFig;
plot(squeeze(testingPathsObj.paths{4}(2,1,:)),'linewidth',2); hold on;
if runALS, plot(squeeze(estALS.estTestPathObj.paths{4}(2,1,:)),'--','linewidth',2); end
if runORALS, plot(squeeze(estORALS.estTestPathObj.paths{4}(2,1,:)),':','linewidth',3); end
axis tight;grid on;
legend({'original','predicted (ALS)','predicted (ORALS)'})


warning on;
