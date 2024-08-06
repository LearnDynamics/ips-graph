function [A_adj,stats_solver] = estimate_graph( B, y, c, varargin)
%
% function [A_adj,stats_solver] = estimate_graph( B, y, c, varargin)
%
% IN:
%   B is a cell array of N cells, each of size N x (d*(L-1)*M) x n
%   y has size N*d*(L-1)*M x 1
%   c has size n

% Computational complexity: constructs and solves N systems of size d(L-1)M x N for a solution of size N, taking total
%   time N_par*min(d(L-1)M x N^2,(d(L-1)M)^2 x N); the construction takes time N_par*Nd(L-1)Mn.
% Memory use: about d(L-1)MN if not parallelized, up to N*d(L-1)MN if fully parallelized (in N).

DEBUG = 0;

%% parse inputs
p = inputParser;
addRequired(p, 'B');
addRequired(p, 'y');
addRequired(p, 'c');
addOptional(p, 'system_params',[]);
addOptional(p, 'reg_method', 'None');
addOptional(p, 'optimOpts', []);
addOptional(p, 'warmStart', []);
parse(p, B, y, c, varargin{:});

reg_method  = p.Results.reg_method;
N           = p.Results.system_params.N;
d           = p.Results.system_params.d;
L           = p.Results.system_params.L;
M           = p.Results.system_params.M;   
optimOpts   = p.Results.optimOpts;
A_0         = p.Results.warmStart';

%% initialize in view of iterations
A_adj       = zeros(N,N);
stats_solver= [];
b           = reshape(permute(reshape( y, [N,d,L-1, M] ),[2,3,4,1]),[d*(L-1)*M,N]);
c           = c';

constraintMat   = -speye(N);
constraintVec   = zeros(N,1);

%% Given B and current guess for the interaction kernel, construct the matrix Bb of the linear system for the unknown A
parfor i = 1:N         
    A_i = squeeze( pagemtimes( c,permute( B{i}, [3,1,2] ) ) )';    
    b_i = b(:,i);

    switch reg_method
        case 'None'
            A_adj(i,:) = A_i\b_i;
        case 'pinv'
            A_adj(i,:) = pinv(A_i)*b_i;
        case 'lsqminnorm'
            A_adj(i,:) = lsqminnorm(A_i, b_i);
        case 'lsqnonneg'
            A_adj(i,:) = lsqnonneg( double(A_i), double(b_i) );
        case 'lsqlin'
            % if DEBUG>0
            %     optimOpts = optimoptions(optimOpts,'Display','iter-detailed');
            % else
            %     optimOpts = optimoptions(optimOpts,'Display','none');
            % end
%            if ~isempty(A_0)
%                [A_adj(i,:),~,~,~,stats_solver{i}] = lsqlin(A_i, b_i, constraintMat, constraintVec, [],[],[],[],A_0(i,:),optimOpts);
%            else
                A_adj(i,:) = lsqlin( double(A_i), double(b_i), constraintMat, constraintVec, [],[],[],[] );
%            end
            % if DEBUG>0
            %     fprintf('\nIters=%d',stats_solver{i}.iterations);
            % end
        case 'ID'
            [~, A_adj(i,:)] = L_curve_standard_form(A_i, b_i, 0);  % to supply the Basis matrix
        case {'RKHS', 'RKHS_plain'}
            [~, A_adj(i,:)] = L_curve(A_i, b_i, 'RKHS', 0);
    end
end

A_adj = A_adj';

end
