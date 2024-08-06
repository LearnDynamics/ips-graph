function [E, c, error_graph, error_kernel, time] = learn_kerne_graph_ALS_social_range(all_xpath, I_true, varargin)

%% Input parser
p = inputParser;

% These are all regularization methods available
expected_reg_methods = {'None','pinv','lsqminnorm', 'ID', 'RKHS', 'RKHS_plain'};

addRequired(p, 'all_xpath');
addRequired(p, 'I_true');
addOptional(p, 'niter', 50);   % sparsity is the percentage of nonzero entries connected to a node
addOptional(p, 'plotON', 0);
addOptional(p, 'normalizeON', 1);
addOptional(p, 'auto_stop_at', 10); % If detect convergence, stop. Threshold  = 1e(auto_stop_at)
addOptional(p, 'reg_method', 'None', @(x) any(validatestring(x, expected_reg_methods)));

parse(p,all_xpath, I_true,varargin{:});
niter = p.Results.niter;
normalizeON = p.Results.normalizeON;
plotON = p.Results.plotON;
auto_stop_at = p.Results.auto_stop_at;
reg_method = p.Results.reg_method;

%% Joint learning using alternating least square:
fprintf('\nALS - ALS\t M = %i ... regu method=%s ,', length(all_xpath),reg_method);tic
% Initialize the kernel and the graph for iteration
% I_true = update_dict_mat(I_true, all_xpath);

I = I_true;

I.UU = randn(I.N, 2);
I.VV = randn(I.n, 2);
I.coef_mat = I.UU*I.VV';
I.E = set_graph(I.N, 'plotON', 0);

error_kernel = zeros(niter+1, 1);
error_graph = zeros(niter+1, 1);

error_kernel(1) = norm(I.coef_mat - I_true.coef_mat, 'fro');
error_graph(1)  = norm(I.E - I_true.E, 'fro');

% iteration

% if 1
%     %% DEBUG
%     estimate_UU_social_range(I_true, all_xpath, 'reg_method', reg_method);
% end

% if 1
%     %% DEBUG
%     estimate_VV_social_range(I_true, all_xpath, 'reg_method', reg_method);
% end
% I1 = I;
%%
% all_I = cell(niter, 1);
for i = 1:niter
%     i
%     all_I{i} = I;
    % Estimate the kernel
    I.E = estimate_graph_social_range(I, all_xpath, 'reg_method', reg_method);
    I.E = I.E ./ sqrt(sum(I.E.^2));
%     error_kernel(i+1)     = sqrt((I.c_true - I_true.c_true)'*I.dict_mat*(I.c_true - I_true.c_true));
    % eigAmat = eig(Amat);        --- returns error because the trajectory is divergent
    
    % Estimate the graph
    I.UU = estimate_UU_social_range(I, all_xpath, 'reg_method', reg_method);
    I.UU = I.UU./sqrt(sum(I.UU.^2));
    
    I.VV = estimate_VV_social_range(I, all_xpath, 'reg_method', reg_method);
    I.coef_mat = I.UU*I.VV';
    
    % normalize the graph
    %     if normalizeON
    %         column_norm = sqrt(sum(I.E.^2, 1));
    %         I.E = I.E./column_norm;
    %     end
    
    switch normalizeON
        case 1
            column_norm = sqrt(sum(I.E.^2, 1));
            I.E = I.E./column_norm;
        case 'mat'
            mat_norm = norm(I.E, 'fro');
            I.E = I.E ./ mat_norm;
    end
    
%     if normalizeON
%         I.UU = 
    %     if auto_s
    
    error_graph(i+1) = norm(I_true.E - I.E, 'fro')./norm(I_true.E, 'fro');
    error_kernel(i+1) = norm(I_true.coef_mat - I.coef_mat, 'fro')./norm(I.coef_mat, 'fro');
    if auto_stop_at && i > 2
        dif = log10(abs(error_graph(i+1) - error_graph(i)));
        if -dif > auto_stop_at
            fprintf('converge at step %i with threshold 1e-%i...', i, auto_stop_at);
            break;
        end
    end
end

% If there is no normalization for each step,
% perform a final normalization step at the end.
%%
% pause;

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
time = toc;
c = I.coef_mat;
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
    y1 = error_kernel;
    y2 = error_graph;
    y3 = obs_std*ones(1, i+1);
    plot(log10(y1),'k-','linewidth',1);
    plot(log10(y2),'b-.','linewidth',1);
    plot(log10(y3),'r--','linewidth',1);
    legend('relMSE Kernel estimator error','relFroSq Graph estimator error','observation noise level');
    ylabel('log_{10} error'); xlabel('n iterations');
    title('ALS')
    ylim([log10(obs_std)-1, 2])
    grid on
    
    %% Want to implement auto stop
    % detect convergence automatically
    % test
    %     plot(diff(log10(error_graph)))
end



end
%%
% figure;hold on;
% fplot(I_true.phi_kernel, [0, 5]);
% fplot(I.phi_kernel, [0, 5]);
% legend('True','est')
% title('Interaction kernel')
% disp('The difference of estimated graph')
% disp(I_true.E - I.E);


