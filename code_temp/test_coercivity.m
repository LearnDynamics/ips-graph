% Test the coercivity
close all;  clear all;   clc;    addpaths;   rng(10);

%% System settings
% Specify the settings here
% In order to fix the figures for paper

I = system_settings(); % setting of the IPS and graph and its integrator
I.viscosity = 0;    % viscosity (noise)
I.N = 6;               % number of agents
I.d = 2;               % dim of state vectors
I.t0       = 0;
I.dt       = 1e-3;     % time steps
I.steps    = 1;      % time steps
I.obs_std  = 0;     % observation noise
I.E = set_graph(I.N, 'sparsity', 0.4, 'plotON', 0);
I.initial = 'N_0_0.1';
I.basis_case = 3;



I = update_system_settings(I);
I.n = length(I.dict);       % Dimension of the hypothesis space
I.phi_kernel = get_kernel_from_c(I.c_true, I.dict);
I.TN = I.steps + 1;                                         % Total time steps
I.tgrid = I.t0:I.dt:(I.steps)*I.dt;
I.X0 = set_particle_initial_all_dim(I.N, I.d, I.initial);   % Initial condition
I.Z_true = get_Z_from_E_c(I.E, I.c_true);     % Z is the product of E and c
I.graph_norm  = norm(I.E,'fro');


%% Coercivity condition for the original IPS
M = 100;
fprintf('\nGenerating trajectories, M = %i...',M);tic
[test_path, I] = get_M_path(I, M, 'ParforProgressON', 1,'saveON', 0, 'loadON', 0);
fprintf('done (%.2f sec).\n',toc);

% % Joint learning using ORALS
% [E_ORSVD, c_ORSVD, E_ORALS, c_ORALS, Z_OR, Aarray, condA_orals] = learn_kernel_graph_ORALS(test_path, I, ...
%     'plotON', 0, 'reg_method', 'RKHS');  % reg_methods: ID, RKHS, None, lsqminnorm
%
% % Joint learning using alternating least square:
% [E_ALS, c_ALS]  = learn_kerne_graph_ALS(test_path, I, 'niter', 10, 'normalizeON', ...
%     1, 'plotON', 0, 'reg_method', 'RKHS');
% [RIP, COND, time, RIP2] = get_RIP_COND(test_path, I);

%%

B = 5000;
M = 5000;

n = I.n;
dict = I.dict;

all_cH = zeros(B, 1);
parfor b = 1:B
    c = randn(n, 1);
    test_kernel = get_kernel_from_c(c, dict);
    X0 = set_particle_initial_all_dim(2*M, d, initial);
    dif = X0(1:M, :) - X0(M+1:2*M, :);
    dis = sqrt(sum(dif.^2, 2));
    samples_KX = test_kernel(dis).*dif./dis;
    E_KX_2 = sum(samples_KX.^2, 'all')/M;
    E_KX = sum(samples_KX, 1)/M;
    var_KX = E_KX_2 - sum(E_KX.^2);
    norm_kernel_2 = c'*dict_mat*c;
    all_cH(b) = var_KX/norm_kernel_2;
    b
end
figure;
histogram(all_cH)

%%
B = 100;
M = 1000;

n = I.n;
dict = I.dict;

all_cH = zeros(B, 1);
for b = 1:B
    c = randn(n, 1);
    test_kernel = get_kernel_from_c(c, dict);
    samples_KX = zeros(m, d);
    for m = 1:M
        X0 = set_particle_initial_all_dim(2, d, initial);
        dif = X0(1, :) - X0;
        dif(1, :) = [];
        dis = sqrt(sum(dif.^2, 2));
        samples_KX(m, :) = test_kernel(dis).*dif./dis;
    end
    E_KX_2 = sum(samples_KX.^2, 'all')/M;
    E_KX = sum(samples_KX, 1)/M;
    var_KX = E_KX_2 - sum(E_KX.^2);
    norm_kernel_2 = c'*dict_mat*c;
    all_cH(b) = var_KX/norm_kernel_2;
end




%%
% B = 100;
% kernel_num = 1000;
% samples_ratio = zeros(kernel_num, 1);
% n = I.n;
% dict_mat = I.dict_mat;
% dict = I.dict;
% d = I.d;
% initial = I.initial;
% parfor k = 1:kernel_num
%     c = randn(n, 1);
%     test_kernel_norm = sqrt(c'*dict_mat*c);
%     c = c/test_kernel_norm;
%     test_kernel = get_kernel_from_c(c, dict);
%
%     samples_KX = zeros(B, 1);
%     for i = 1:B
%         X0 = set_particle_initial_all_dim(B, d, initial);
%         dif = X0(1, :) - X0;
%         dis = sum(dif.^2, 2);
%         dis = dis(2:end, :);
%         samples_KX(i) = var(test_kernel(dis));
%     end
%     samples_ratio(k) = mean(samples_KX);
%     k
% end
%
% %
% figure;histogram(log10(samples_ratio));
% title('distribution of coercivity constant')
% xlabel('$\log_{10}(c_\mathcal{H})$','fontweight', 'bold', 'Interpreter', 'latex');


%% test
% M = 1000;
% N = I.N;
% d = I.d;
% initial = I.initial;
% obs_std = I.obs_std;
% parfor i = 1:M
%     X0 = set_particle_initial_all_dim(N, d, initial);
%     xpath        = graph_forward_model(I, X0, 0, 0);      % using RK4, the result will not converge.
%     all_path{i} = xpath + obs_std* randn(size(xpath));
% end
% [N, d, steps_] = size(test_path{1});
% steps = steps_-1;
%
%
% B = 1;
% all_ratio = zeros(B, 1);
% for b = 1:B
%     c = randn(I.n, 1);
%     test_kernel = get_kernel_from_c(c, dict);
%     v = zeros(N-1, d, steps, M);
%     for m = 1:M
%         for t = 1:steps
%             X = all_path{m}(:, :, t);
%             dif = X(:, :) - X(1, :);
%             dif(1, :) = [];
%             dis = sum(dif.^2, 2);
%             s = test_kernel(dis).*dif./dis;
%             v(:, :, t, m) = s;
%         end
%     end
%
%     mat = zeros(N-1, N-1);
%     for i = 1:N-1
%         for j = 1:N-1
%             mat(i, j) = sum(v(i, :, :, :).*v(j, :, :, :), 'all')/(steps*M);
%         end
%     end
%
%     rho_norm = sum(v.^2, 'all')/((N-1)*steps*M);
%     mat = mat/rho_norm;
%     s = svd(mat);
%     all_ratio(b) = max(s)/min(s);
%
%     b
% end

%% test
% clc
% B = 100000;
%
% S = zeros(B, 1);
% T = zeros(B, 1);
% X = randn(4, B);
% S = test_kernel(X(3, :) - X(1, :));
% T = test_kernel(X(2, :) - X(1, :));
% cov(S, T)