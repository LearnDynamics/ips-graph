function [v_hat, A_kernel, y,condA] = ALS_3_estimate_v( B, y, A_cur, u_hat, dict_mat, varargin)

%
% function [c, A_kernel, y] = estimate_kernel( B, y, A_cur, varargin)
%
% IN:
%   B : cell array of N cells, each of size N x (d*(L-1)*M) x n, as constructed by build_B_tensor
%   y : NdLM vector from observations (velocities, or approximations thereof), will be right-hand of A_kernel c = y
%
% OUT:
%   c           : n vector of coefficients for interaction kernel
%   A_kernel    : NdLM x n matrix for linear system A_kernel c = y
%

% Computational cost:   Solves a linear system of size (dLMN^2)*n, at a computational cost of (dLMN^2)*n^2, assuming n<dLMN. 
%                       Assembling of the systems takes time N_par*NdLMn.

% (c) Mauro Maggioni 9/2023

p = inputParser;
addRequired(p, 'B');
addRequired(p, 'y');
addRequired(p, 'A_cur');
addOptional(p, 'system_params',[]);
addOptional(p, 'reg_method', 'None');
addOptional(p, 'test_regON', '0');
parse(p, B, y, A_cur, varargin{:});
reg_method = p.Results.reg_method;
test_regON = p.Results.test_regON; 
system_params = p.Results.system_params;

N           = system_params.N;
d           = system_params.d;
L           = system_params.L;
M           = system_params.M;
[~, Q]      = size(u_hat);
b           = reshape(permute(reshape( y, [N,d,L-1, M] ),[2,3,4,1]),[d*(L-1)*M,N]);


% Generate A and b for linear system for estimation of interaction kernel
[N,dLM,n] = size(B{1});
A_adj   = A_cur';


% debug
% c_mat   = system_params.c_mat;
% u       = system_params.u;
% v       = system_params.v;
%% Compute the product of the pre-computed basis kernels on data with the adjacency matrix
v_hat = zeros(Q, N);
parfor i = 1:N                                                                                                                                                                                                                                    
    A_i    = squeeze( pagemtimes( A_adj(i,:), B{i} ) )*u_hat;
    b_i    = b(:, i);
    switch reg_method
        case 'None'
            v_hat(:, i) = A_i\b_i;
        case 'pinv'
            v_hat(:, i) = pinv(A_i)*b_i;
        case 'lsqminnorm'
            v_hat(:, i) = lsqminnorm(A_i, b_i);
        case 'lsqlin'
            % if DEBUG>0
            %     optimOpts = optimoptions(optimOpts,'Display','iter-detailed');
            % else
            %     optimOpts = optimoptions(optimOpts,'Display','none');
            % end
%            if ~isempty(A_0)
%                [v_i,~,~,~,stats_solver{i}] = lsqlin(A_i, b_i, constraintMat, constraintVec, [],[],[],[],A_0(i,:),optimOpts);
%            else
                v_hat(:, i) = lsqlin( double(A_i), double(b_i), constraintMat, constraintVec, [],[],[],[] );
%            end
            % if DEBUG>0
            %     fprintf('\nIters=%d',stats_solver{i}.iterations);
            % end
        case 'ID'
            [~, v_hat(:, i)] = L_curve_standard_form(A_i, b_i, 0);  % to supply the Basis matrix
        case {'RKHS', 'RKHS_plain'}
            [~, v_hat(:, i)] = L_curve(A_i, b_i, 'RKHS', 0);
    end
    
end
v_hat = v_hat';
end