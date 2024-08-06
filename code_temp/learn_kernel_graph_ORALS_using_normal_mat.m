function [E, c, E2, c2, all_Z, all_A, condA_orals, condL, time] = learn_kernel_graph_ORALS(path_data, dyn_sys, learning_setup, varargin)
% joint learning of kernel and graph from data trajectories
p = inputParser;

% These are all regularization methods available
expected_reg_methods = {'None', 'pinv', 'pinvreg','lsqminnorm', 'ID', 'RKHS', 'RKHS_plain'};

testON = 0; % if 1, compute the condition number of A and (A,B);

addRequired(p,'path_data');
addRequired(p,'dyn_sys');
addOptional(p, 'plotON', 0);
addOptional(p, 'normalizeON', 1);
addOptional(p, 'reg_method', 'None', @(x) any(validatestring(x, expected_reg_methods)));
addOptional(p, 'pinvregeps',eps)

parse(p, path_data, dyn_sys, varargin{:});
plotON      = p.Results.plotON;
normalizeON = p.Results.normalizeON;
reg_method  = p.Results.reg_method;
pinvregeps  = p.Results.pinvregeps;

M       = length(path_data);  % path_data is a cell array Mx1
n       = learning_setup.n;
N       = dyn_sys.N;
d       = dyn_sys.d;
L       = dyn_sys.L;
dt      = dyn_sys.dt;
basis   = learning_setup.dict;

% AA = zeros(N*n*(N-1), n*(N-1));     % Store the matrix A for N particles
% Ab = zeros(N*n*(N-1), 1);           % Store the vector b for N particles
% % This is because the graph has N columns, and each time we only estimate
% % one column of it. But in the parfor loop, we must extract information for
% % all columns.
%fprintf('\nORALS \t M = %i ...', M);
time = tic;


all_A = zeros(n*(N-1), n*(N-1), N);
all_b = zeros(n*(N-1), 1, N);

parfor m = 1:M
    % For the m-th trajectory
    X = path_data{m};
    A_M = zeros(n*(N-1), n*(N-1), N);
    b_M = zeros(n*(N-1), 1, N);
    
    for i = 1:N
        % consider the i-th particle
        Xi = X(i, :, :);
        % Using the trajectory of the i-th particle to estimate the i-th
        % column of the graph and the kernel
        
        dXi = (X(i,:, 2:end)- X(i,:,1:end-1))/dt;
        dif = Xi(:, :, 1:end-1) - X(:, :, 1:end-1);
        dis = sqrt(sum(dif.^2, 2));
        
        all_K = zeros(d*(L-1), n*(N-1));
        for k = 1:n
            f = basis{k};
            K_phi = (f(dis)./dis).*dif;
            K_phi(isnan(K_phi)) = 0;
            %             all_K(k, :, :, :) = K_phi;
            K_phi_reshape = reshape(K_phi, [N, d*(L-1)]);
            % Original shape : N*d*L
            % changed into shape: N*(d*L)
            % The first coordinates stays the same.
            % Combine the 2nd, 3rd coordinates, in the order of
            % X(d1, t1), X(d2, t1), X(d1, t2), X(d2, t2), X(d1, t3), X(d2, t3);
            K_phi_reshape(i, :) = [];
            K_phi_reshape = K_phi_reshape';
            all_K(:, (k-1)*(N-1)+1:k*(N-1)) = K_phi_reshape;
        end
        
        % all_K store the regression matrix for ORALS
        % where the columns are of the order
        % X(d1, t1), X(d2, t1), X(d1, t2), X(d2, t2), X(d1, t3), X(d2, t3)...
        % The rows are of the order
        % particle first, then basis
        dXi_reshape = reshape(dXi, [], 1);
        
        A_i = all_K'*all_K;
        b_i = all_K'*dXi_reshape;
        
        % A_i and b_i are the regression matrix/vector for estimating the
        % i-th column of the graph and the kernel
        % derived from the m-th trajectory data
        A_M(:, :, i) = A_i;
        b_M(:, :, i) = b_i;
    end
    all_A = all_A + A_M;
    all_b = all_b + b_M;
end

all_A = all_A/M; all_b = all_b/M;

%% Solve linear system and get graph and coefficients



% get all_Z
% Z is the product of graph and coef
all_Z = zeros(N-1, n, N);

graph_mat   = eye(dyn_sys.N-1);
dyn_sys.ORALS_mat = kron( learning_setup.dict_mat, graph_mat );                                                                                      %MM: if this is ORALS-specific, it should not be computed here, nor be put in I

condA_orals = zeros(1, N);
condL = zeros(1, N);

for i = 1:N
    %extract the matrix
    A_i = all_A(:, :, i);
    b_i = all_b(:, :, i);
    % solve linear system
    
    switch reg_method
        case 'None'
            Z_i = A_i\b_i;
        case 'pinv'
            Z_i = pinv(A_i)*b_i;
        case 'pinvreg'
            Z_i = pinv(A_i,pinvregeps)*b_i;            
        case 'lsqminnorm'
            Z_i = lsqminnorm(A_i, b_i);
        case 'ID'
            [~, Z_i] = L_curve_standard_form(A_i, b_i, 0);  % to supply the Basis matrix
        case 'RKHS'
            [~, Z_i] = L_curve(A_i, b_i, 'RKHS', 0, dyn_sys.ORALS_mat);
        case 'RKHS_plain'
            [~, Z_i] = L_curve(A_i, b_i, 'RKHS', 0);
    end
    
    if testON==1
        condA_orals(i) = cond(A_i);
        warning off
        condL(i) = cond(dyn_sys.ORALS_mat\A_i);
        warning on
    end
    
    % reshape
    all_Z(:, :, i) = reshape(Z_i, [N-1, n]);
end
if testON==1
    fprintf('regu Method = %s, Average condA = 10^%2.0f, condL = 10^%2.0f   ...',reg_method, log10(mean(condA_orals)),log10(mean(condL)));
else
    %fprintf('regu Method = %s,   ...',reg_method);
end

%% SVD for Z

E = zeros(N, N);
all_coef = zeros(n, N);

for i = 1:N
    Z_i = all_Z(:,:,i);
    [U, S, V] = svd(Z_i);
    graph = U(:, 1);
    coef = V(:, 1)*S(1,1);
    
    % make the graph to be positive
    [~, ind] = max(abs(graph));
    sgn = sign(graph(ind));
    
    graph = graph/sgn;
    coef = coef/sgn;
    
    % store the output
    E(i, :) = [graph(1:(i-1));0;graph(i:end)];
    all_coef(:, i) = coef;
end
%%
% figure;hold on
% ylim([-2, 2])
% fplot(get_kernel_from_c(I.c_true, I.dict), [0, 5])
% for i = 1:10
%     fplot(get_kernel_from_c(all_coef(:, i), I.dict), [0, 5])
% end






%%

switch normalizeON
    case {1, 0}
        E = E';
        c = mean(all_coef, 2);
    case 'mat'
        c_norm = sqrt(sum(all_coef.^2, 1));
        all_coef = all_coef./c_norm;
        c = mean(all_coef, 2);
        E = E';
        E = E.*c_norm;
        
        E_norm = norm(E, 'fro');
        E = E/E_norm;
        c = c*E_norm;
end

% How to use the property that E is positive?

%% ALS starting from SVD result
niter = 5;
[graph2,  coef2] = ALS_inORALS( all_Z, c, niter, learning_setup, plotON, normalizeON );

E2 = 0*E;
for i =1:N
    % make the graph to be positive
    graph0= graph2(:,i);
    [~, ind] = max(abs(graph0));
    sgn = sign(graph0(ind));
    graph0 = graph0/sgn;
    
    E2(i, :) = [graph0(1:(i-1));0;graph0(i:end)];
end

E2 = E2';            % transpose here: Because svd gives column vectors, and the graph uses row representation
c2 = coef2;

time = toc( time );

%fprintf('done (%.2f sec).\n',time);

end


function [E, c] = ALS_inORALS(Z_arrays, c0, niter, learning_setup, plotON, normalizeON)
% the alternating least squares part of ORALS:
% estimate array a and vector c, with c shared between rows of a;


% I tried this code with obs_std = 0. It does not coverge to the true
% value. I will check the details.

% Numerical results shows that, ORALS is always worse than ORSVD.

[n, ~, N] = size(Z_arrays);

% err_g = zeros(niter, 1);
err_k = zeros(niter+1, 1);
err_k(1) = kernel_err(c0, learning_setup);
for t = 1:niter
    E = zeros(n, N);
    c = zeros(size(c0));
    switch normalizeON
        case {1, 0}
            for i = 1:N
                Z = Z_arrays(:, :, i);
                u_half = Z*c0;
                if normalizeON
                    u1 = u_half/norm(u_half);
                else
                    u1 = u_half;
                end
                E(:,i) = u1;
                c = c + Z'*u1;
            end
            c = c/N;
            
        case 'mat'
            for i = 1:N
                Z = Z_arrays(:, :, i);
                E(:,i) = Z*c0;
            end
            
            E_norm = norm(E, 'fro');
            E = E/E_norm;
            
            all_c = zeros(length(c0), N);
            for i = 1:N
                Z = Z_arrays(:, :, i);
                u1 = E(:,i);
                all_c(:, i) = Z'*u1./(norm(u1).^2);
            end
            
            %             c_norm = sqrt(sum(all_c.^2, 1));
            %             all_c = all_c./c_norm;
            c = mean(all_c, 2);
            
    end
    err_k(t+1) = kernel_err( c, learning_setup );
    c0 = c;
end

%%
if plotON
    figure;
    plot(err_k);
    title('The change of the kernel error with the iteration of ORALS')
end
end


%
%             all_c = zeros(length(c0), N);
%             for i = 1:N
%                 Z = Z_arrays(:, :, i);
%                 u1 = Z*c0;
%                 E(:,i) = u1;
%                 %                 c = c + Z'*u1;
%                 all_c(:, i) = c;
%             end
%             c_norm = sqrt(sum(all_c.^2, 1));
%             all_c = all_c./c_norm;
%             c = mean(all_c, 2);
%
%             E = E.*c_norm;
%             E_norm = norm(E, 'fro');
%             E = E/E_norm;
%             c = c*E_norm;
%
%             all_c = zeros(length(c0), N);
%             for i = 1:N
%                 Z = Z_arrays(:, :, i);
%                 u1 = Z*c0;
%                 E(:,i) = u1;
%                 %                 c = c + Z'*u1;
%                 all_c(:, i) = c;
%             end
