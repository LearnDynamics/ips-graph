function [c, all_A, all_b] = estimate_VV_social_range(I, X, varargin)
%% Inpute paser
p = inputParser;
addRequired(p, 'I');
addRequired(p, 'X');
addOptional(p, 'reg_method', 'None');

parse(p, I, X, varargin{:});
reg_method = p.Results.reg_method;

%% Generate A and b for least square estimation of VV
if ~isa(X, 'cell')
    X = {X};        % For single trajectory case
end
lengthX = length(X);

n = I.n;

all_A = zeros(2*n, 2*n);
all_b = zeros(2*n, 1);

for i = 1:lengthX
    [A_N, b_N] = get_A_b_single_traj(I, X{i});
    all_A = all_A + A_N;
    all_b = all_b + b_N;
end

N = I.N;
% U_est = zeros(N, 2);
% VV = I.VV;
% VV_line = reshape(VV, [], 1);

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

c = [c(1:I.n), c(I.n+1:end)];
% norm(all_A*VV_line - all_b)
% norm(c - VV_line)
end



function [AA, bb] = get_A_b_single_traj(I, X)

%% Generate A and b for least square estimation of kernel
d = I.d;
N = I.N;
TN = I.TN;
dt = I.dt;
E = I.E;
U = I.UU;
% VV = I.VV;
% VV_line = reshape(VV, [], 1);
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
% id = 4;
% U = I.UU;
% temp = zeros(d, TN-1);
% ZZ = cell(2*n, 1);
% for i = 1:n
%     for r = 1:2
%         Z = squeeze(U(id, r)*K_basis_conv_U_at_traj{i}(id, :, 1:end-1));
%         ZZ{i+(r-1)*n} = Z;
%         temp = temp + I.VV(i, r)*Z;
%     end
% end
%
%
%
% temp
% squeeze(dX(id, :, :))
%
% % ZZZ = {ZZ{:, 1},ZZ{:, 2}};
%
% A = zeros(24, 24);
% b = zeros(24, 1);
% for i = 1:24
%     for j = 1:24
%         A(i, j) = sum(ZZ{i}.*ZZ{j}, 'all');
%     end
% end
%
% for i = 1:24
%     b(i) = sum(ZZ{i}.*squeeze(dX(id, :, :)), 'all');
% end
%
%
% norm(A*VV_line - b)
% norm(A\b - VV_line)
% a = 1;
%% debug

% %
% id = 10;
% U = I.UU(id, :);
% temp = zeros(d, TN-1);
% for i = 1:2
%     temp = temp + squeeze(U(i)*VV_K_basis_conv_U_at_traj{i}(id, :, 1:end-1));
% end
% temp
% squeeze(dX(id, :, :))

%%
%%
A_all = zeros(2*n, 2*n, N);
b_all = zeros(2*n, 1, N);




% id = 4;
% U = I.UU;
% temp = zeros(d, TN-1);
% ZZ = cell(n, 2);
% for i = 1:n
%     for r = 1:2
%         Z = squeeze(U(id, r)*K_basis_conv_U_at_traj{i}(id, :, 1:end-1));
%         ZZ{i, r} = Z;
%         temp = temp + I.VV(i, r)*Z;
%     end
% end
% temp
% squeeze(dX(id, :, :))
%
% ZZZ = {ZZ{:, 1},ZZ{:, 2}};
%
% A = zeros(24, 24);
% b = zeros(24, 1);
% for i = 1:24
%     for j = 1:24
%         A(i, j) = sum(ZZZ{i}.*ZZZ{j}, 'all');
%     end
% end
%
% for i = 1:24
%     b(i) = sum(ZZZ{i}.*squeeze(dX(id, :, :)), 'all');
% end
%
%
%
% a = 1;


% k is the number of particles
for k = 1:N
    
    % compute the part U times R
    ZZ = cell(2*n, 1);
    for i = 1:n
        for r = 1:2
            Z = squeeze(U(k, r)*K_basis_conv_U_at_traj{i}(k, :, 1:end-1));
            ZZ{i+(r-1)*n} = Z;
            %             temp = temp + I.VV(i, r)*Z;
        end
    end
    
    % compute the matrix A and b
    
    for i = 1:2*n
        for j = 1:2*n
            A(i, j) = sum(ZZ{i}.*ZZ{j}, 'all');
        end
    end
    
    for i = 1:2*n
        b(i) = sum(ZZ{i}.*squeeze(dX(k, :, :)), 'all');
    end
    
    %
    %     A = zeros(2*n, 2*n);
    %     b = zeros(2*n, 1);
    %     for r = 1:2
    %         for s = 1:2b
    %             for i = 1:n
    %                 for j = 1:n
    %                     id_1 = i+(r-1)*n;
    %                     id_2 = j+(s-1)*n;
    %                     A(id_1, id_2) = sum(UU(k, r)*K_basis_conv_U_at_traj{i}(k,:, 1:end-1) .* UU(k, s).*K_basis_conv_U_at_traj{j}(k,:, 1:end-1), 'all');
    % %                     A(id_2, id_1) = A(id_1, id_2);
    %                 end
    %             end
    %         end
    %     end
    %
    %     for r = 1:2
    %         for i = 1:2
    %             id_1 = i+(r-1)*n;
    %             b(id_1) = sum(UU(k, r)*K_basis_conv_U_at_traj{i}(k,:, 1:end-1).*dX(k, :, :), 'all');
    %         end
    %     end
    
    A_all(:, :, k) = A;
    b_all(:, :, k) = b;
    %     k
    
    %     norm(A*VV_line - b)
    % norm(A\b - VV_line)
    
    
end

a = 1;
%% debug

id = 4;
AA = sum(A_all, 3);
bb = sum(b_all, 3);
% norm(AA*VV_line - bb)
% norm(AA\bb - VV_line)
end

