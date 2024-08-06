function [u, A_kernel, y, condA] = ALS_3_estimate_u( B, y, A_cur, v, dict_mat, varargin)

%
% function [c, A_kernel, y] = ALS_3_estimate_u( B, y, A_cur, v, varargin)
% Three fold ALS in the learning of multiple type kernel



% IN:
%   B : cell array of N cells, each of size N x (d*(L-1)*M) x n, as constructed by build_B_tensor
%   y : NdLM vector from observations (velocities, or approximations thereof), will be right-hand of A_kernel c = y
%
% OUT:
%   u           : coefficient matrix of shape n x Q,
%   A_kernel    : NdLM x n matrix for linear system A_kernel c = y
%

% Computational cost:   Solves a linear system of size (dLMN^2)*n, at a computational cost of (dLMN^2)*n^2, assuming n<dLMN.
%                       Assembling of the systems takes time N_par*NdLMn.

% (c) Mauro Maggioni, Quanjun Lang 1/2024

p = inputParser;
addRequired(p, 'B');
addRequired(p, 'y');
addRequired(p, 'A_cur');
addRequired(p, 'v');
addOptional(p, 'system_params',[]);
addOptional(p, 'reg_method', 'None');
addOptional(p, 'test_regON', '0');
parse(p, B, y, A_cur, v, varargin{:});
reg_method = p.Results.reg_method;
test_regON = p.Results.test_regON;
system_params = p.Results.system_params;

N           = system_params.N;
d           = system_params.d;
L           = system_params.L;
M           = system_params.M;
[~, Q]      = size(v);
b           = reshape(permute(reshape( y, [N,d,L-1, M] ),[2,3,4,1]),[d*(L-1)*M,N]);

% Generate A and b for linear system for estimation of interaction kernel
[N,dLM,n] = size(B{1});
A_adj   = A_cur';

%%%%%%%%%%%%%%%%%%%%%%%%%%% Debug settings
% c_mat   = system_params.c_mat;
% u       = system_params.u;
% v       = system_params.v;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Debug settings

%% Compute the product of the pre-computed basis kernels on data with the adjacency matrix
for i = 1:N
    A_i    = squeeze( pagemtimes( A_adj(i,:), B{i} ) );
    A_reshape_i = reshape(A_i', [n, 1, dLM]);
    Av_i = pagemtimes(A_reshape_i, v(i, :));

    %sanity check:
    % squeeze(sum(Av{i}.*u, [1, 2]))
    % b(:, i)

    % because the above two equations equals, we reshape Av to be a matrix
    B{i} = reshape(Av_i, [n*Q, dLM])';
end

A_u     = vertcat(B{:});
b_u     = reshape(b, [N*dLM, 1]);


switch reg_method
    case 'None'
        u_temp = A_u\b_u;
    case 'pinv'
        u_temp = pinv(A_u)*b_u;
    case 'lsqminnorm'
        u_temp = lsqminnorm(A_u, b_u);
    case 'lsqnonneg'
        u_temp = lsqnonneg( double(A_u), double(b_u) );
    case 'lsqlin'

        u_temp = lsqlin( double(A_u), double(b_u), constraintMat, constraintVec, [],[],[],[] );

    case 'ID'
        [~, u_temp] = L_curve_standard_form(A_u, b_u, 0);  % to supply the Basis matrix
    case {'RKHS', 'RKHS_plain'}
        [~, u_temp] = L_curve(A_u, b_u, 'RKHS', 0);
end

u = reshape(u_temp, [n, Q]);

end