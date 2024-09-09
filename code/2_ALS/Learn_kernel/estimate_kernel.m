function [c, A_kernel, y,condA] = estimate_kernel( B, y, A_cur,dict_mat, varargin)

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

MIN_N_PARFOR = Inf;                                                                                                             % minimum value of N for which code is parallelized

p = inputParser;
addRequired(p, 'B');
addRequired(p, 'y');
addRequired(p, 'A_cur');
addOptional(p, 'system_params',[]);
addOptional(p, 'reg_method', 'None');
addOptional(p, 'test_regON', '0');
parse(p, B, y, A_cur, varargin{:});
reg_method      = p.Results.reg_method;
test_regON      = p.Results.test_regON; 
system_params   = p.Results.system_params;
% Generate A and b for linear system for estimation of interaction kernel
[N,dLM,n] = size(B{1});
A_adj   = A_cur';

%% Compute the product of the pre-computed basis kernels on data with the adjacency matrix
if N > MIN_N_PARFOR                                                                                                             % extra memory but in parallel
    parfor i = 1:N
        B{i} = squeeze( pagemtimes( A_adj(i,:), B{i} ) );
    end    
    A_kernel = reshape(horzcat(B{:})',[n,N*dLM])';                                                                              % reshape to form the matrix for the least squares problem
else                                                                                                                            % without extra memory, but harder to parallelize        
    A_kernel   = zeros(N*dLM,n);
    for i = 1:N
        A_kernel(i:N:N*dLM,:) = pagemtimes( A_adj(i,:), B{i} );
    end
end

% Solve the least squares problem for the (coefficients of the) new estimated interaction kernel
switch reg_method
    case 'None'
        c = A_kernel\y;
    case 'pinv'
        c = pinv(A_kernel)*y;
    case 'lsqminnorm'
        c = lsqminnorm(A_kernel, y);
    case 'ID'    % later, we should make it using A, b
        Abar = A_kernel'* A_kernel; bbar = A_kernel'*y; 
        [~, c] = L_curve_standard_form(Abar, bbar, 0);  % to supply the Basis matrix
%          [~, c] = L_curve_standard_form(A_kernel, y, 0);
    case 'RKHS'
        Abar = A_kernel'* A_kernel; bbar = A_kernel'*y; 
        [~, c] = L_curve(Abar, bbar, 'RKHS', 0, dict_mat);                                                                      % MM: TBD: dict_mat porbably not passed correctly, or can't handle rectangular least squares problems
    case 'RKHS_plain'
        Abar = A_kernel'* A_kernel; bbar = A_kernel'*y; 
        [~, c] = L_curve(Abar, bbar, 'RKHS', 0);
end

condA = 0; 
if test_regON ==1
    condA = cond(A_kernel); 
end
end