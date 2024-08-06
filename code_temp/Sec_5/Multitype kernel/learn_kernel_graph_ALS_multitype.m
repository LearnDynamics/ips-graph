function est = learn_kernel_graph_ALS_multitype( training_paths , dyn_sys, learning_set, varargin)

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
addOptional(p, 'stop_thres_relDelta_u', Inf);
addOptional(p, 'stop_thres_relDelta_v', Inf);
addOptional(p, 'stop_thres_relDelta_A', Inf);
addOptional(p, 'stop_thres_testpaths', Inf);
addOptional(p, 'reg_methodK', 'None', @(x) any(validatestring(x, expected_reg_methods)));
addOptional(p, 'reg_methodA', 'None', @(x) any(validatestring(x, expected_reg_methods)));
addOptional(p, 'test_paths',[])
addOptional(p, 'K_means_on_v', false)

parse(p,training_paths, dyn_sys,varargin{:});
niter           = p.Results.niter;
normalizeON     = p.Results.normalizeON;
plotON          = p.Results.plotON;
reg_methodK     = p.Results.reg_methodK;
reg_methodA     = p.Results.reg_methodA;
A_sparsity_thres= p.Results.A_sparsity_thres;
A_sparsity_prior= p.Results.A_sparsity_prior;
test_paths      = p.Results.test_paths;
K_means_on_v    = p.Results.K_means_on_v;


stop_thres_testpaths    = p.Results.stop_thres_testpaths;
stop_thres_relDelta_u   = p.Results.stop_thres_relDelta_u;
stop_thres_relDelta_v   = p.Results.stop_thres_relDelta_v;
stop_thres_relDelta_A   = p.Results.stop_thres_relDelta_A;


% initialize for the iterations
est.time.ALS        = tic;
system_params       = struct( 'N',dyn_sys.N,'d',dyn_sys.d,'M',length(training_paths),'L',size(training_paths{1},3) );

Q                   = learning_set.num_kernel_choices;
u_hat               = randn(learning_set.n, Q);   
v_hat               = randn(dyn_sys.N, Q);        
A_hat               = set_graph(dyn_sys.N, 'plotON', 0);

est.time.build_B_tensor       = tic;
[B,dX]                        = build_B_tensor( dyn_sys, learning_set, training_paths );
est.time.build_B_tensor       = toc(est.time.build_B_tensor);




est.stats.u_hat(:,:,1)          = u_hat;        
est.stats.v_hat(:,:,1)          = v_hat;    
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
        optimOptsA  = optimset( 'TolX',10*max(size(u_hat))*1e-4 );
    case 'lsqlin'
        optimOptsA  = optimoptions( 'lsqlin','TolX',10*max(size(u_hat))*1e-4 );           
    otherwise
        optimOptsA  = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% Debug settings
% u_hat = learning_set.UU;
% v_hat = learning_set.VV;
% A_hat = dyn_sys.A;
% system_params.u         = learning_set.UU;
% system_params.v         = learning_set.VV;
% system_params.c_mat     = learning_set.coef_mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% ALS iteration
dict_mat = learning_set.dict_mat;
for q = 1:niter
    u_hat               = ALS_3_estimate_u   ( B, dX, A_hat, v_hat, dict_mat,'reg_method', reg_methodK, 'system_params', system_params, 'test_regON',1);         % Estimate the interaction kernel% Normalize the graph    
    v_hat               = ALS_3_estimate_v   ( B, dX, A_hat, u_hat, dict_mat,'reg_method', reg_methodK, 'system_params', system_params,'test_regON',1);         % Estimate the interaction kernel
    v_hat               = normalize_v        ( v_hat , K_means_on_v);
    A_hat               = ALS_3_estimate_graph    ( B, dX, u_hat, v_hat, 'reg_method', reg_methodA, 'system_params', system_params,...
                                                            'optimOpts', optimOptsA, 'warmStart', [] );                         % Estimate the graph
    A_hat               = normalizeAdj      ( A_hat, normalizeON, A_sparsity_thres, A_sparsity_prior );                         % Normalize the graph    



    est.stats.v_hat(:,:,q+1)        = v_hat;
    est.stats.u_hat(:,:,q+1)        = u_hat;
    est.stats.coef_mat(:,:,q+1)     = u_hat*v_hat';
    est.stats.A_hat(:,:,q+1)        = A_hat;
    % est.condA_kernel(q+1)  = condA; 
    % stopping criterion
    if q >= 2
        relDelta_u = norm(est.stats.u_hat(:,q+1)-est.stats.u_hat(:,q) )/norm(est.stats.u_hat(:,q));
        relDelta_v = norm(est.stats.v_hat(:,q+1)-est.stats.v_hat(:,q) )/norm(est.stats.v_hat(:,q));
        relDelta_A = norm(est.stats.A_hat(:,:,q+1)-est.stats.A_hat(:,:,q),'fro' )/norm(est.stats.A_hat(:,:,q),'fro');
        stop_satisfied = relDelta_A<stop_thres_relDelta_A & relDelta_u<stop_thres_relDelta_u & relDelta_v<stop_thres_relDelta_v;

        if stop_thres_testpaths<Inf && ~isempty(test_paths)                % using test paths error for stopping  
            [est.stats.pathTestErr(q+1),~,est.stats.meanL2traj(q+1)]    = getPathTestError( dyn_sys, learning_set, A_hat, c_hat, test_paths );  % if test paths are provided, evaluate current estimator on them
            stop_satisfied = stop_satisfied & abs(est.stats.pathTestErr(q+1))*est.stats.meanL2traj(q+1) < stop_thres_testpaths;
        end

        if stop_satisfied, break; end
    end
end


%%%%%%%%%%%%%%%%%%%%%%% final Kmeans for v, since we know there are Q types of kernels
[N, Q] = size(v_hat);[ind, x] = kmeans(v_hat, Q);
v_hat = zeros(N, Q);
for i = 1:N;v_hat(i, :) = x(ind(i), :);end



est.time.ALS      = toc( est.time.ALS );
est.ALS_n_iter    = q+1; 
est.A             = A_hat;
est.u             = u_hat;
est.v             = v_hat;
est.coef_mat      = u_hat*v_hat';
est.kernel_idx    = ind;



if stop_thres_testpaths>0 && ~isempty(test_paths)
    est.pathTestErr   = est.stats.pathTestErr(end);
    est.meanL2traj    = est.stats.meanL2traj(end);
end

fprintf('\n3 fold ALS number of iterations: %d  ', q); 
end
