function est =  learn_kernel_graph_ORALS_B( path_data,dyn_sys,matFactorize, learning_set, varargin)
% [E, c, all_Z, condA_orals, condL, est] = ... %
% Estimate graph and coefs by first operator regression, then matrix factorization
% OUT
%   est : a structure with the following fields:
%       A       : estimator for the graph
%       c       : estimator for the coefficients of the interaction kernel
%       time    : a structure with various timing information

% Quanjun Lang, Mauro Maggioni, Fei Lu (c), 2023

reg_testON = 0;    % if 1, compute the condition number of A and (A,B);
%% Input parser
p      = inputParser;
% These are all regularization methods available
expected_reg_methods = {'None', 'pinv', 'pinvreg','lsqminnorm', 'ID', 'RKHS', 'RKHS_plain'};

addRequired(p,'path_data');
addRequired(p,'dyn_sys');
addOptional(p, 'plotON', 0);
addOptional(p, 'normalizeON', 1);
addOptional(p, 'A_sparsity_thres', 0);
addOptional(p, 'A_sparsity_prior', 1);                                                                                          % sparsity is the percentage of nonzero entries connected to a node
addOptional(p, 'reg_method', 'None', @(x) any(validatestring(x, expected_reg_methods)));
addOptional(p, 'test_paths',[])

parse(p,path_data, dyn_sys,varargin{:});
normalizeON     = p.Results.normalizeON;
reg_method      = p.Results.reg_method;
A_sparsity_thres= p.Results.A_sparsity_thres;
A_sparsity_prior= p.Results.A_sparsity_prior;

%% Build B tensor from path data
time_orals        = tic;

% Build B, a cell array of N cells, each of size N x (d*(L-1)*M) x n <<< B{i} = ( phi_k(X_t^i-X_t^j) ):
time_OR      = tic;
[B,dX]                        = build_B_tensor_orals( dyn_sys, learning_set, path_data );
est.time_build_B_tensor       = toc(time_OR);
%%%   B{i}      = B_tmp(j,dtm,n)
%%%   a(:,i)*B{i}*c  = a_j*B{i}(j,dtm,n)*c_n   >> z_jn = a_j*c_n,  z is (N-1) x p 
%%%   BB*z      = sum_jn BB(dtm,jn)*z_jn  
%%%  dX has size N*d*(L-1)*M x 1

%% Operator regression: get all_Z 
time_OR      = tic;
[all_Z,N,n,condA_orals,condL] = operator_regression(B,dX,learning_set,dyn_sys,reg_method,reg_testON); 
est.time_solveOR      = toc(time_OR);

%% Matrix factorization to get A, c from Z_i's
time_factorizeZ  = tic; 
switch matFactorize 
    case 'SVD'
        [Esvd,csvd]    = factorizeZs_svd(N,n,all_Z,normalizeON);
    case 'ALS'; niter  = 10;
        [Esvd, csvd]    = factorizeZs_svd(N,n,all_Z,normalizeON);
        [E,c,Eseq,cseq] = ALS_inORALS(all_Z, csvd, niter, normalizeON);
end

        % normalize the adjacency matrix estimator  ------ enhance sparcity, reduce the errors 
        [Esvd,~]   = normalizeAdj( Esvd, normalizeON, A_sparsity_thres, A_sparsity_prior );
        [E,~]      = normalizeAdj( E, normalizeON, A_sparsity_thres, A_sparsity_prior );
        
% If there is no normalization for each step, perform a final normalization step now
if ~normalizeON
    [E,~]            = normalizeAdj( E, normalizationType, A_sparsity_thres, A_sparsity_prior );
    [E,c,Eseq,cseq]  = ALS_inORALS_initial_A(all_Z, E, 1, normalizeON);
end
est.time_factorizeZ = toc(time_factorizeZ);

%% Output
est.time_ORALS    = toc( time_orals);
est.A             = E;
est.c             = c;
est.phi_kernel    = combine_dict_coef ( learning_set, c );
est.condA_orals   = condA_orals; est.condL  = condL; 
est.matFactorize  = matFactorize;
est.all_Z         = all_Z; 
est.Esvd  = Esvd';      est.csvd = csvd; 
est.Eseq  = Eseq;      est.cseq = cseq; 
end




