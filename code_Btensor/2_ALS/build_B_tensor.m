function [B,dX,dX_L] = build_B_tensor( dyn_sys, learning_setup, paths )

%
% Construct internally a data-driven tensor B_tmp of size N x d x L x M x N x n, and l.h.s. vector y of size N*d*(L-1)*M x 1
%   B is a cell array of N cells, each of size N x (d*(L-1)*M) x n
%   dX has size N*d*(L-1)*M x 1
% Computational cost: M*N^2*d*L*n with a constant likely dominated by the cost of evaluating a basis function at a point,
%   and all the way down to N^2*d*L*n with ideal parallelization over M, or, which is the current implementation,
%   down to M*N*d*L*n with ideal parallelization over N.
%

% (c) Mauro Maggioni

M       = length(paths);
d       = dyn_sys.d;
N       = dyn_sys.N;
L       = dyn_sys.L;
basis   = learning_setup.dict;
n       = length(basis);

%% Create tensor B_tmp and cell array B
B   = cell( N,1 );

parfor i = 1:N                                                                                                                  % Run through agents
    Btmp = zeros( N,d,L-1,M,n, 'single' );
    for m = 1:M                                                                                                                 % Run through the trajectories
        X   = paths{m};                                                                                                         %#ok<PFBNS>
        Xi  = X(i, :, :);
        dif = Xi - X;
        dif = dif(:,:,1:L-1);
        dis = sqrt(sum(dif.^2, 2));
        divdistimesdif = dif./dis;
        %z=cellfun( @(f) f(dis).*divdistimesdif,basis,'UniformOutput',false );                                                  % evaluate all basis functions at dis, obtain terms in the sum on the r.h.s. of i-th equation in the ODE system, and return a d x L x N tensor of results
        for k = 1:n                                                                                                             % Run through basis functions
            f                   = basis{k};
            K_phi               = f(dis).*divdistimesdif;
            K_phi(isnan(K_phi)) = 0;
            Btmp(:,:,:,m,k)     = single(K_phi);                                                                                % note: K_phi has size N x d x L
            %z{k}(isnan(z{k})) = 0;
            %B_tmp(:,:,:,i,m,k)  = single(z{k} );
        end
    end
    Btmp = permute( Btmp, [1,4,2,3,5] );
    B{i} = double(reshape(Btmp,[N,M*d*(L-1),n]));
end

%% Create the dX vector
dX          = zeros(N,M,d,L-1);
dt          = dyn_sys.dt;
parfor m    = 1:M
    dX(:,m,:,:) = (paths{m}(:,:,2:end)-paths{m}(:,:,1:end-1))/dt;
end

% temp = dX;

dX = reshape(dX,[dyn_sys.N*M*dyn_sys.d*(L-1),1]);

%% Create the dX_L vector
if nargout>=3
    dX_L        = zeros(N,M,d,L-1);
    dt          = dyn_sys.T/dyn_sys.L;
    parfor m    = 1:M
        dX_L(:,m,:,:) = (paths{m}(:,:,2:end)-paths{m}(:,:,1:end-1))/dt;
    end
    dX_L = reshape(dX_L,[dyn_sys.N*M*dyn_sys.d*(L-1),1]);
end


% dX_L = temp;
end
