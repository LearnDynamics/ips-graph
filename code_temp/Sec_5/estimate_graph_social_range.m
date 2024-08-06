function E = estimate_graph_social_range(I, X, varargin)
% This code estimate the graph given a kernel
% The input can be single or multiple trajectories

p = inputParser;
addRequired(p, 'I');
addRequired(p, 'X');
addOptional(p, 'reg_method', 0);

parse(p, I, X, varargin{:});
reg_method = p.Results.reg_method;

%% load parameters
N = I.N;
%% If the input is a single trajectory
if ~isa(X, 'cell')
    X = {X};        % For single trajectory case
end
B = length(X);


%%
batch_A = zeros(N-1, N-1, N);
batch_b = zeros(N-1, 1, N);
for i = 1:B
    [all_A, all_b] = get_all_A_b_graph(I, X{i});
    batch_A = batch_A + all_A;
    batch_b = batch_b + all_b;
end

%% get graph
E_temp = zeros(N-1, N);
for i = 1:N
    A = batch_A(:, :, i);
    b = batch_b(:, :, i);
    
%     %     E_temp(:, i) = A\b;
%     E_temp(:, i) = pinv(A)*b;
%     
    switch reg_method
        case 'None'
            E_temp(:, i) = A\b;
        case 'pinv'
            E_temp(:, i) = pinv(A)*b;
        case 'lsqminnorm'
            E_temp(:, i) = lsqminnorm(A, b);
        case 'ID'
            [~, E_temp(:, i)] = L_curve_standard_form(A, b, 0);  % to supply the Basis matrix
        case {'RKHS', 'RKHS_plain'}
            [~, E_temp(:, i)] = L_curve(A, b, 'RKHS', 0);
    end
    
    
    
end

E = zeros(N, N);
for i = 1:N
%     Z_i = all_Z(:,:,i);
%     [U, S, V] = svd(Z_i);
    graph = E_temp(:, i);
%     coef = V(:, 1)*S(1,1);
    
    % make the graph to be positive
    [~, ind] = max(abs(graph));
    sgn = sign(graph(ind));
    
    graph = graph/sgn;
%     coef = coef/sgn;
    
    % store the output
    E(:, i) = [graph(1:(i-1));0;graph(i:end)];
%     all_coef(:, i) = coef;
end

% E = E';
% E = reform_E(E_temp);



end







function [all_A, all_b] = get_all_A_b_graph(I, X)
% apply the least square for each particle.
% There are N-1 edges start from the i-th particle
% Hence we need to apply optimization N times
% Each time we optimize N-1 parameters.

%% load parameters
N = I.N;
dt = I.dt;
%
coef_mat = I.coef_mat;
%% loop
% Store N linear inverse problems
% with A of the size N-1 * N-1;
% and b of the size 1 * N-1;
all_A = zeros(N-1, N-1, N);
all_b = zeros(N-1, 1, N);


for k = 1:N
    dX = (X(k,:, 2:end)- X(k,:,1:end-1))/dt;
    Xk = X(k, :, :);
    dif = Xk - X;
    dis = sqrt(sum(dif.^2, 2));
    
    phi_kernel = get_kernel_from_c(coef_mat(k, :), I.dict);% The kernel for the k-th particle
    K = phi_kernel(dis)./dis.*dif;
    K(k, :, :) = [];
    A = zeros(N-1, N-1);
    for i = 1:N-1
        for j = 1:i
            A(i, j) = sum(K(i, :, 1:end-1).*K(j, :, 1:end-1), 'all');
            A(j, i) = A(i, j);
        end
    end
    
    b = zeros(N-1, 1);
    for i = 1:N-1
        b(i) = sum(K(i, :, 1:end-1).*dX, 'all');
    end
    
    
    all_A(:, :, k) = A;
    all_b(:, :, k) = b;
end
end