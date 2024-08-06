function K = Kuramoto_main(dyn_sys, learning_setup, M, reg_ORALS)
fprintf('\nGenerating trajectories, M = %i...',M);tic
trainingPathsObj = get_paths( dyn_sys, M, 'ParforProgressON', 1,'saveON', 0, 'loadON', 0 );
fprintf('done (%.2f sec).\n',toc);


learning_setup                = update_dict_mat( learning_setup, trainingPathsObj.paths);

rho = learning_setup.rho;
rho_dx = learning_setup.rho_dx;
xgrid = learning_setup.rho_bin_edges(2:end);
N = length(rho);
n = learning_setup.n;
%%
% N = length(I.rho);
% xgrid = I.rho_dx:I.rho_dx:I.rho_dx*N;
% N = length(xgrid);
b = sin(xgrid);
A = zeros(N, n);
for i = 1:n
    A(:, i) = learning_setup.dict{i}(xgrid);
end
W = diag(rho);

learning_setup.c                = (A'*W*A)\(A'*W*b');
learning_setup.kernel_norm      = sqrt(sum(sin(xgrid).^2.*rho*rho_dx));  
%% Joint learning using ORALS
estALS = learn_kernel_graph_ALS( trainingPathsObj.paths, dyn_sys, learning_setup, 'niter', 10, 'normalizeON', 1, ...
    'stop_thres_testpaths', dyn_sys.viscosity*10, 'stop_thres_relDelta_A', 1e-3, 'stop_thres_relDelta_c', 1e-3, 'plotON', 1, ...
    'reg_methodK', 'None', 'reg_methodA', 'lsqnonneg', ...
    'A_sparsity_thres', 0,'A_sparsity_prior',1, 'test_paths', [] );% 1e-2*norm(I.A,'fro') );
graph_err ( estALS.A, dyn_sys ),
kernel_err( estALS.c, learning_setup ),


matFactorize  = 'ALS';    % %  'SVD', 'ALS', 'SVD+ORALS'
estORALS = learn_kernel_graph_ORALS_B(trainingPathsObj.paths, dyn_sys, matFactorize, ...
    learning_setup, 'plotON', 0, 'reg_method', reg_ORALS);

%%
K.M = M;
K.estALS = estALS;
K.estORALS = estORALS;
K.learning_setup = learning_setup;

% err.graph_ALS = graph_err(E_ALS, I);
% err.graph_ORALS = graph_err(E_ORALS, I);
% err.kernel_ALS = kernel_err_non_para(c_ALS, I);
% err.kernel_ORALS = kernel_err_non_para(c_ORALS, I);
% 
% %% Check kernel estimation
% K.phi_kernel_ALS = get_kernel_from_c(c_ALS, I.dict);
% K.phi_kernel_ORALS = get_kernel_from_c(c_ORALS, I.dict);
% K.phi_kernel_c_true = get_kernel_from_c(I.c_true, I.dict);
% K.I = I;
% K.M = M;
% K.err = err;
% 
% figure;hold on;grid on;
% 
% rgrid = I.rgrid;
% plot(rgrid, K.phi_kernel_ALS(rgrid),'linewidth',1);
% plot(rgrid, K.phi_kernel_ORALS(rgrid),'linewidth',1);
% plot(rgrid, I.phi_kernel(rgrid),':','linewidth',2);
% % plot(rgrid, phi_kernel_c_true(rgrid),'linewidth',2);
% 
% ylim([-1.2, 1.2])
% 
% % yyaxis right
% % plot(rgrid, I.rho)
% 
% fill_x = [rgrid, fliplr(rgrid)];
% fill_y = [I.rho, 0*rgrid];
% fill(fill_x, fill_y, [0.9, 0.9, 0.9],'LineStyle','none');
% 
% xlim([I.rho_dx, L])
% ylim([-1.2, 1.2])
% legend( 'ALS', 'ORALS', 'True','rho', 'Location', 'southwest')
% 
% 
% set(gca, 'SortMethod', 'depth');
% set_positionFontsAll;
% title(['M = ', num2str(M)])
% %%
% figname = [I.PAPER_FIG_DIR, '/Kuramoto_M_', num2str(M), '.pdf'];
% saveas(gcf, figname);
% close all

end