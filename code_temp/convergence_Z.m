% Convergence of graph and kernel error
%% Test Convergence in M for both ALS and ORALS
% generate data of multple trajectories
Mseq        = 10*2.^(1:1:4);
MM = Mseq(end)*10;  % Generate more data and sample from these data when testing convergence

fprintf('The sequence of M is');disp(Mseq);
fprintf('Generating trajectories, M = %i...', MM);tic
all_xpath = get_M_path(I, MM, 'ParforProgressON', 1,'saveON', 1,'loadON', 0);
fprintf('done (%.2f sec).\n',toc);

Mseq_length = length(Mseq);
batch_num = 5;
error_Z_orals  = zeros(batch_num, Mseq_length);
error_Z_orsvd  = zeros(batch_num, Mseq_length);
error_Z_als = zeros(batch_num, Mseq_length);
Z_true = get_Z_from_E_c(I.E, I.c_true);

%%
for i = 1:Mseq_length
    i
    M = Mseq(i);
    for b = 1:batch_num
        path_id = randperm(MM, M);
        test_path = all_xpath(path_id); % sample from a large collection of path
        
        %         if b == 1
        %             plotON = 1;
        %         else
        %             plotON = 0;
        %         end
        
        [E_ORSVD, c_ORSVD, E_ORALS, c_ORALS, Z_OR] = learn_kernel_graph_ORALS(test_path, I);
        Z_ORALS = get_Z_from_E_c(E_ORALS, c_ORALS);
        error_Z_orals(b, i) = get_Z_error(Z_ORALS, Z_true);
        Z_ORSVD = get_Z_from_E_c(E_ORSVD, c_ORSVD);
        error_Z_orsvd(b, i) = get_Z_error(Z_ORSVD, Z_true);
        
        error_k_orals(i, b) = kernel_err(c_ORALS, I);
        error_k_orsvd(i, b) = kernel_err(c_ORSVD, I);
%         error_k_als(i, b) = kernel_err(c_ALS, I);
        
        error_g_orals(i, b) = graph_err(E_ORALS, I);
        error_g_orsvd(i, b) = graph_err(E_ORSVD, I);
        
        
        %         [E_ALS, c_ALS]  = learn_kerne_graph_ALS(test_path, I, 'niter', 20, 'normalizeON', 1, 'plotON', plotON);
        %         Z_ALS = get_Z_from_E_c(E_ALS, c_ALS);
        %         error_Z_als(b, i) = get_Z_error(Z_ALS, Z_true);
    end
end

%%
figure;hold on;
% plot(log10(Mseq), log10(error_Z_als), 'b.', 'LineWidth', 2, 'MarkerSize', 12);
% f=get(gca,'Children');f_als = f(1);
plot(log10(Mseq), log10(error_Z_orals), 'r.', 'LineWidth', 2, 'MarkerSize', 12);
f=get(gca,'Children');f_orals = f(1);
plot(log10(Mseq), log10(error_Z_orsvd), 'y.', 'LineWidth', 2, 'MarkerSize', 12);
f=get(gca,'Children');f_orsvd = f(1);
% legend([f_als,f_orals, f_orsvd], 'ALS', 'ORALS', 'ORSVD')
xlabel('log_{10} M')
ylabel('log_{10} Z error')
title('Compare the estimation of Z')


%%
% save('data.mat', 'error_Z_als', 'error_Z_or')