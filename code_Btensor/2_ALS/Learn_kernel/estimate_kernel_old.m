function [c, all_A, all_b] = estimate_kernel(I, X, varargin)
%

p = inputParser;
addRequired(p, 'I');
addRequired(p, 'X');
addOptional(p, 'reg_method', 'None');

parse(p, I, X, varargin{:});
reg_method = p.Results.reg_method;

%% Generate A and b for least square estimation of kernel
if ~isa(X, 'cell')
    X = {X};        % For single trajectory case
end
lengthX = length(X);

all_A = zeros(I.n, I.n);
all_b = zeros(I.n, 1);

for i = 1:lengthX
    [A, b] = get_A_b_single_traj(I, X{i});
    all_A = all_A + A;
    all_b = all_b + b;
end

switch reg_method
    case 'None'
        c = all_A\all_b;
    case 'pinv'
        c = pinv(all_A)*all_b;
    case 'lsqminnorm'
        c = lsqminnorm(all_A, all_b);
    case 'ID'
        [~, c] = L_curve_standard_form(all_A, all_b, 0);  % to supply the Basis matrix
    case 'RKHS'
        [~, c] = L_curve(all_A, all_b, 'RKHS', 0, I.dict_mat);
    case 'RKHS_plain'
        [~, c] = L_curve(all_A, all_b, 'RKHS', 0);
end




end


function [A, b] = get_A_b_single_traj(I, X)

%% Generate A and b for least square estimation of kernel
% Change the shape of [X*d, tN] to [N, d, tn]
d = I.d;
N = I.N;
TN = I.TN;
dt = I.dt;
E = I.E;
% X = reshape(xpath, [d, N, TN]);
% X = permute(X,[2 1 3]);

%% load parameters
basis = I.dict;
n = length(basis);
K_basis_conv_U_at_traj = cell(n, 1);
dX   = (X(:,:, 2:end)- X(:,:,1:end-1))/dt;
for i = 1:n
    f = basis{i};
    K_conv_U_at_traj = zeros(N, d, TN);
    for k = 1:N
        Xk = X(k, :, :);
        dif = Xk - X;
        dis = sqrt(sum(dif.^2, 2));
        % nonlinear function
        K_phi = (f(dis)./dis).*dif;
        K_phi(isnan(K_phi)) = 0;
        K_conv_U_at_traj(k, :, :) = squeeze(sum(K_phi.*E(:, k)));
    end
    K_basis_conv_U_at_traj{i} = K_conv_U_at_traj;
end

%% debug
% Xk = X(1, :, :);
% dif = Xk - X;
% dis = sqrt(sum(dif.^2, 2));
% temp = I.phi_kernel(dis)./dis.*dif;
% temp = temp.*E(:, k);
% temp(isnan(temp)) = 0;
% temp = squeeze(sum(temp, 1));
%
% test = 0;
% for i = 1:n
%     test = test + I.c_true(i) * K_basis_conv_U_at_traj{i};
% end
%%
A = zeros(n, n);
b = zeros(n, 1);

for i = 1:n
    for j = 1:i
        A(i,j) = sum(K_basis_conv_U_at_traj{i}(:,:, 1:end-1) .* K_basis_conv_U_at_traj{j}(:,:, 1:end-1), 'all');
        A(j,i) = A(i,j);
    end
end

% dX = gradient(obs_path);
% NO, use Euler approximation: obs_path(t+1)-obs_path(t): here we have stochastic force, the gradient uses central difference

for i = 1:n
    b(i)  = sum(K_basis_conv_U_at_traj{i}(:,:, 1:end-1).*dX, 'all');
end

% Const = sum(dX.^2, 'all');
% c = A\b;

end










%
% function K_conv_U_at_traj = get_function_conv_U_at_traj_Monte_Carlo(sysInfo, xpath, obs_path, f, dist_info)
% % This function approximate the term K_phi * U(Xk, t) by Monte Carlo
% % Integration
% [a, TN] = size(obs_path);
% d = sysInfo.d;
% K = floor(a/d);
% % N = sysInfo.N;
% [s1, ~] = size(xpath);
% N = s1/d;
%
%
% if exist('dist_info','var')
%     obs_all_diff = dist_info.obs_all_diff;
%     obs_all_dist = dist_info.obs_all_dist;
%
%
%     temp_val = f(obs_all_dist)./obs_all_dist;
%     temp_vec = obs_all_diff;
%     for k = 1:d
%         temp_vec(k:d:end, :, :) = temp_vec(k:d:end, :, :).*temp_val;
%     end
%     temp_vec(isnan(temp_vec)) = 0;
%     K_conv_U_at_traj = reshape(sum(temp_vec, 2), a, TN)/N;
%
%
% else
%
%     K_conv_U_at_traj = zeros(d*K, TN);
%     for k = 1:K
%         Xk = obs_path((k-1)*d+1:k*d, :);
%
%         Xk_Xn = repmat(Xk, N, 1) - xpath;
%
%         % distance between the given particle and other particles
%         Xk_Xn_distance = zeros(N, TN);
%         for j = 1:d
%             Xk_Xn_distance = Xk_Xn_distance + Xk_Xn(j:d:end,:).^2;
%         end
%         Xk_Xn_distance = sqrt(Xk_Xn_distance);
%         % nonlinear function
%         phi_distance_divide_distance = f(Xk_Xn_distance)./Xk_Xn_distance;
%
%         temp = Xk_Xn;
%         for j = 1:d
%             temp(j:d:end,:) = temp(j:d:end,:) .* phi_distance_divide_distance;
%         end
%         s = sum(isnan(temp))/d;
%         temp(isnan(temp)) = 0;
%
%
%
%         for j = 1:d
%             K_conv_U_at_traj((k-1)*d+j, :) = sum(temp(j:d:end, :), 1)./(N - s);
%         end
%     end
% end
% end
%
