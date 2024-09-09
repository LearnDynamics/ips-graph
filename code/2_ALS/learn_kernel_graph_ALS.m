function est = learn_kernel_graph_ALS( training_paths , dyn_sys, learning_set, varargin)

%
% OUT
%   est : a structure with the following fields:
%       A       : estimator for the graph
%       c       : estimator for the coefficients of the interaction kernel
%       time    : a structure with various timing information
%

% Quanjun Lang, Mauro Maggioni (c), 2023
%

clear estimate_kernel2

%% Input parser
p = inputParser;

% These are all regularization methods available
expected_reg_methods = {'None','pinv','lsqminnorm', 'lsqnonneg', 'lsqlin','ID', 'RKHS', 'RKHS_plain'};

addRequired(p, 'training_paths');
addRequired(p, 'I_true');
addOptional(p, 'niter', 50);   % sparsity is the percentage of nonzero entries connected to a node
addOptional(p, 'plotON', 0);
addOptional(p, 'normalizeON', 1);
addOptional(p, 'A_sparsity_thres', 0);
addOptional(p, 'A_sparsity_prior', 1);                                                                                          % sparsity is the percentage of nonzero entries connected to a node
addOptional(p, 'stop_thres_relDelta_c', Inf);
addOptional(p, 'stop_thres_relDelta_A', Inf);
addOptional(p, 'stop_thres_testpaths', Inf);
addOptional(p, 'reg_methodK', 'None', @(x) any(validatestring(x, expected_reg_methods)));
addOptional(p, 'reg_methodA', 'None', @(x) any(validatestring(x, expected_reg_methods)));
addOptional(p, 'test_paths',[])

parse(p,training_paths, dyn_sys,varargin{:});
niter           = p.Results.niter;
normalizeON     = p.Results.normalizeON;
plotON          = p.Results.plotON;
reg_methodK     = p.Results.reg_methodK;
reg_methodA     = p.Results.reg_methodA;
A_sparsity_thres= p.Results.A_sparsity_thres;
A_sparsity_prior= p.Results.A_sparsity_prior;
test_paths      = p.Results.test_paths;

stop_thres_testpaths    = p.Results.stop_thres_testpaths;
stop_thres_relDelta_c   = p.Results.stop_thres_relDelta_c;
stop_thres_relDelta_A   = p.Results.stop_thres_relDelta_A;


%% Initialization for the iterations
est.time.ALS        = tic;
system_params       = struct( 'N',dyn_sys.N,'d',dyn_sys.d,'M',length(training_paths),'L',size(training_paths{1},3) );
c_hat               = randn(learning_set.n, 1);                                                                                 % not needed as first update is going to be on c_hat, given A_hat
A_hat               = set_graph(dyn_sys.N, 'plotON', 0);

%% Construct data tensor
est.time.build_B_tensor       = tic;
[B,dX]                        = build_B_tensor( dyn_sys, learning_set, training_paths );
est.time.build_B_tensor       = toc(est.time.build_B_tensor);

est.stats.c_hat(:,1)          = c_hat;                                                                                          % save info about this step
est.stats.A_hat(:,:,1)        = A_hat;

est.dyn_sys                   = dyn_sys;
if ~isempty(test_paths)
    est.dyn_sys.initial       = zeros(length(test_paths),dyn_sys.N,dyn_sys.d);
    for m = 1:length(test_paths)
        est.dyn_sys.initial(m,:,:) = test_paths{m}(:,:,1);
    end
end
if stop_thres_testpaths>0 && ~isempty(test_paths)
    est.time.testpaths = tic;
    [est.stats.pathTestErr(1),~,est.stats.meanL2traj(1)] = getPathTestError( dyn_sys, learning_set, A_hat, c_hat, test_paths ); % if test paths are provided, evaluate current estimator on them
    est.time.testpaths = toc(est.time.testpaths);
end

switch(reg_methodA)                                                                                                             % optimization options
    case 'lsqnonneg'
        optimOptsA  = optimset( 'TolX',10*max(size(c_hat))*1e-4 );
    case 'lsqlin'
        optimOptsA  = optimoptions( 'lsqlin','TolX',10*max(size(c_hat))*1e-4 );           
    otherwise
        optimOptsA  = [];
end

%% ALS iterations
dict_mat = learning_set.dict_mat;
for q = 1:niter
    [c_hat,~,~,condA]   = estimate_kernel   ( B, dX, A_hat, dict_mat,'reg_method', reg_methodK, 'system_params', ...
                                                            system_params,'test_regON',1);                                      % Estimate the interaction kernel
    A_hat               = estimate_graph    ( B, dX, c_hat, 'reg_method', reg_methodA, 'system_params', system_params,...
                                                            'optimOpts', optimOptsA, 'warmStart', [] );                         % estimate the graph
    A_hat               = normalizeAdj      ( A_hat, normalizeON, A_sparsity_thres, A_sparsity_prior );                         % normalize the graph    

    est.stats.c_hat(:,q+1)        = c_hat;                                                                                      % save info about this step
    est.stats.A_hat(:,:,q+1)      = A_hat;
    est.condA_kernel(q+1)  = condA; 

    % check stopping criterion
    if q >= 2
        relDelta_c = norm(est.stats.c_hat(:,q+1)-est.stats.c_hat(:,q) )/norm(est.stats.c_hat(:,q));
        relDelta_A = norm(est.stats.A_hat(:,:,q+1)-est.stats.A_hat(:,:,q),'fro' )/norm(est.stats.A_hat(:,:,q),'fro');
        stop_satisfied = relDelta_c<stop_thres_relDelta_c & relDelta_A<stop_thres_relDelta_A;

        if stop_thres_testpaths<Inf && ~isempty(test_paths)                                                                     % using test paths error for stopping  
            [est.stats.pathTestErr(q+1),~,est.stats.meanL2traj(q+1)] = ...
                                            getPathTestError( dyn_sys, learning_set, A_hat, c_hat, test_paths );                % if test paths are provided, evaluate current estimator on them
            stop_satisfied = stop_satisfied & abs(est.stats.pathTestErr(q+1))*est.stats.meanL2traj(q+1) < stop_thres_testpaths;
        end

        if stop_satisfied, break; end
    end
end

%% It there is no normalization at each step, perform a final normalization step now
if ~normalizeON
    A_hat           = normalizeAdj      ( A, 1 );
    [c_hat,~,~]     = estimate_kernel   ( B, dX, A_hat, 'reg_method', reg_methodK, 'system_params', system_params);             % update the interaction kernel after the normalization of the graph
end

est.time.ALS      = toc( est.time.ALS );
est.ALS_n_iter    = q; 
est.A             = A_hat;
est.c             = c_hat;
est.phi_kernel    = combine_dict_coef ( learning_set, c_hat );
if stop_thres_testpaths>0 && ~isempty(test_paths)
    est.pathTestErr   = est.stats.pathTestErr(end);
    est.meanL2traj    = est.stats.meanL2traj(end);
end

if plotON==1
    fprintf('ALS number of iterations: %d  ', q);
end
end
