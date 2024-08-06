function [E, c, error_graph, error_kernel] = learn_kerne_graph_ALS_optimized(all_xpath, I, Gamma, dX, varargin)

%% Input parser
p = inputParser;
addRequired(p, 'all_xpath');
addRequired(p, 'I');
addRequired(p, 'Gamma');
addRequired(p, 'dX');
addOptional(p, 'niter', 50);   % sparsity is the percentage of nonzero entries connected to a node
addOptional(p, 'plotON', 0);
addOptional(p, 'normalizeON', 1);
addOptional(p, 'auto_stop_at', 1e-12); % If detect convergence, stop

parse(p,all_xpath, I, Gamma, dX, varargin{:});
niter = p.Results.niter;
normalizeON = p.Results.normalizeON;
plotON = p.Results.plotON;

n = I.n;
N = I.N;
d = I.d;
steps = I.steps;
M = length(all_xpath);




%% Joint learning using alternating least square:
fprintf('Joint learning using ALS ...');tic



%% Estimate the kernel
tic
AA = zeros(n, n);
bb = zeros(n, 1);
for m = 1:M
    for i = 1:N
        for t = 1:steps
            A = zeros(d, n);
            for k = 1:n
                temp = squeeze(Gamma(i, :, m, k, t, :)).*I.E(:, i);
                temp(isnan(temp)) = 0;
                A(:, k) = sum(temp, 1)';
            end
            
            b = squeeze(dX(i, m, t, :));
            AA = AA + A'*A;
            bb = bb + A'*b;
        end
    end
end
AA = AA/(N*steps*M);
bb = bb/(N*steps*M);
            

toc




% Initialize the kernel and the graph for iteration
I = I_true;
I.c_true = randn(I.n, 1);
I.E = set_graph(I.N, 'plotON', 0);

error_kernel = zeros(niter, 1);
error_graph = zeros(niter, 1);

% iteration
for i = 1:niter
    % Estimate the kernel
    [I.c_true,~,~]      = estimate_kernel(I, all_xpath);
    I                   = combine_dict_coef(I);
    error_kernel(i)     = sqrt((I.c_true - I_true.c_true)'*I.dict_mat*(I.c_true - I_true.c_true));
    % eigAmat = eig(Amat);        --- returns error because the trajectory is divergent
    
    % Estimate the graph
    I.E                 = estimate_graph(I, all_xpath);
    error_graph(i)      = norm(I_true.E - I.E, 'fro');
    
    % normalize the graph
    if normalizeON
        column_norm = sqrt(sum(I.E.^2, 1));
        I.E = I.E./column_norm;
    end
    
    %     if auto_s
    
end

% It there is no normalization for each step,
% perform a final normalization step at the end.
%%
if ~ normalizeON
    column_norm = sqrt(sum(I.E.^2, 1));
    disp('If no normalization in each step, the final column norms are')
    disp(column_norm)
    mean_column_norm = mean(column_norm);
    I.E = I.E./mean_column_norm;
    I.c_true = I.c_true * mean_column_norm;
    I = combine_dict_coef(I);
    
    
    error_kernel_final = sqrt((I.c_true - I_true.c_true)'*I.dict_mat*(I.c_true - I_true.c_true));
    error_graph_final = norm((I_true.E - I.E).^1, 'fro');
    
    disp('The kernel error after normalization is ')
    disp(error_kernel_final);
    
    disp('The graph error after normalization is ')
    disp(error_graph_final);
    
end

fprintf('done (%.2f sec).\n',toc);

c = I.c_true;
E = I.E;
%% plot the decay of the errors
% figure;hold on;
% plot(log10(error_kernel),'k-','linewidth',1);
% plot(log10(error_graph),'b-.','linewidth',1);
% legend('Kernel estimator error','Graph estimator error');
% ylabel('log_{10} error'); xlabel('n iterations')


%%
if plotON == 1
    figure;hold on;
    obs_std = I.obs_std;
    y1 = error_kernel/sqrt((I_true.c_true'*I.dict_mat*I_true.c_true));
    y2 = error_graph/(norm(I_true.E,'fro'));
    y3 = obs_std^2*ones(1,niter);
    plot(log10(y1),'k-','linewidth',1);
    plot(log10(y2),'b-.','linewidth',1);
    plot(log10(y3),'r--','linewidth',1);
    legend('relMSE Kernel estimator error','relFroSq Graph estimator error','observation noise level');
    ylabel('log_{10} error'); xlabel('n iterations');
    title('ALS')
    
    %% Want to implement auto stop
    % detect convergence automatically
    % test
    plot(diff(log10(error_graph)))
end

%%
% figure;hold on;
% fplot(I_true.phi_kernel, [0, 5]);
% fplot(I.phi_kernel, [0, 5]);
% legend('True','est')
% title('Interaction kernel')
% disp('The difference of estimated graph')
% disp(I_true.E - I.E);

% function get_c_from_E(Gamma, dX, n, M, N, )
%% Estimate the kernel
AA = zeros(n, n);
bb = zeros(n, 1);
for m = 1:M
    for i = 1:N
        for t = 1:steps
            A = zeros(d, n);
            for k = 1:n
                temp = squeeze(Gamma(i, :, m, k, t, :)).*I.E(:, i);
                temp(isnan(temp)) = 0;
                A(:, k) = sum(temp, 1)';
            end
            
            b = squeeze(dX(i, m, t, :));
            AA = AA + A'*A;
            bb = bb + A'*b;
        end
    end
end
AA = AA/(N*steps*M);
bb = bb/(N*steps*M);
end