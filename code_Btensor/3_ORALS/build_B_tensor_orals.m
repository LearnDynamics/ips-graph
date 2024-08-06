function [B,dX,dX_L] = build_B_tensor_orals( dyn_sys, learning_setup, paths )
%
% Construct internally a data-driven tensors
%   B is a cell array of N cells, each B{i} is of size N x (d*(L-1)*M) x n
%          B_tmp(j,d,t,m,n)[i]  = (basis{k}(X^i_t-X^j_t))  >>> the i-th Btmp(j,d,t,m,n)
%          B{i}                = B_tmp(j,dtm,n)
%          a*B{i}*c            = sum_j sum_n a(j,i)*B{i}(j,dtm,n)*c_n
%   dX has size N x d*(L-1)*M 
% Computational cost: M*N^2*d*L*n with a constant likely dominated by the cost of evaluating a basis function at a point,
%       - parallel: can be either in M or N 
% (c) Mauro Maggioni, Fei Lu

%  Difference from build_B_tensor.m in ALS: parallel in M; dX as array; %%  MM's permuted reshape of B 5D tensor: not used here

M       = length(paths);
d       = dyn_sys.d;
N       = dyn_sys.N;
L       = dyn_sys.L;
basis   = learning_setup.dict;
n       = length(basis);

%% Create tensor B_tmp and cell array B
B   = cell( N,1 );
dLM = d*(L-1)*M; 
for i = 1:N                                                                            % Run through agents   -- parallel 
    Btmp  = zeros( N,d,L-1,M,n); % Btmp = zeros( N,d,L-1,M,n, 'single' );
    parfor m = 1:M    % can use parfor here                                            % Run through the M-trajectories
        for k = 1:n                                                                    % Run through basis functions
            X   = paths{m};
            Xi  = X(i, :, :);
            dif = Xi - X;
            dif = dif(:,:,1:L-1);
            dis = sqrt(sum(dif.^2, 2));
            normal_dif= dif./dis; 
            f                   = basis{k};
            K_phi               = f(dis).*normal_dif;                          
            K_phi(isnan(K_phi)) = 0;                                          % note: K_phi has size N x d x L
            Btmp(:,:,:,m,k)    = K_phi;
        end
    end

    B{i} = double(reshape(Btmp,[N,dLM,n]));                              % N, M, d, L, p >> [N,M*d*(L-1),n])
end

%% Create the dX array
dX          = zeros(N,d,L-1,M); 
dt          = dyn_sys.dt;
parfor m    = 1:M
    dX_temp     = (paths{m}(:,:,2:end)-paths{m}(:,:,1:end-1));
    dX(:,:,:,m) = dX_temp; 
end
dX = dX/dt; 

dX = reshape(dX,[N,dLM]);                                  % N, M, d, L >> [N,M*d*(L-1)]  

% verify if dX(i,:)= A(:,i)'* B{i} *c
testON = 0; 
if testON ==1 && dLM < 100    % verified: dX(i,:)= A(:,i)'* B{i} *c 
    a_test_Btensor(dyn_sys,learning_setup,Btmp,dX,n,N,B,2)
end

%% Create the dX_L array
if nargout>=3
    dX_L        = zeros(N,d,L-1,M);
    dt          = dyn_sys.T/dyn_sys.L;   % This step is not used --- debatable, return to the setting of the observation
    parfor m    = 1:M
        dX_L(:,:,:,m) = (paths{m}(:,:,2:end)-paths{m}(:,:,1:end-1))/dt;
    end
    dX_L = reshape(dX_L,[N,dLM]);
end

end
