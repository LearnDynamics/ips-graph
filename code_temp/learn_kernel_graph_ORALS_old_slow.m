function [E, c, E2, c2, Z_arrays, A_arrays, condA_orals,time] = learn_kernel_graph_ORALS_old_slow(path_data, I, varargin)
% joint learning of kernel and graph from data trajectories
p = inputParser;

% These are all regularization methods available
expected_reg_methods = {'None','pinv','lsqminnorm', 'ID', 'RKHS'};

testON = 1; % if 1, compute the condition number of A and (A,B); 

addRequired(p,'path_data');
addRequired(p,'I');
addOptional(p, 'plotON', 0);
addOptional(p, 'reg_method', 'None', @(x) any(validatestring(x, expected_reg_methods)));

parse(p, path_data, I, varargin{:});
plotON = p.Results.plotON;
reg_method = p.Results.reg_method;

M = length(path_data);  % path_data is a cell array Mx1
n = I.n;
N = I.N;
d = I.d;
steps = I.steps;
dt = I.dt;

AA = zeros(N*n*(N-1), n*(N-1));     % Store the matrix A for N particles
Ab = zeros(N*n*(N-1), 1);           % Store the vector b for N particles
% This is because the graph has N columns, and each time we only estimate
% one column of it. But in the parfor loop, we must extract information for
% all columns.


fprintf('\nORALS M = %i ...', M); 
tic
parfor m = 1:M            % to QL: can we reduce the for-loops? 
    path = path_data{m};
    for t = 1:steps
        AA_all_particle = [];
        Ab_all_particle = [];
        for i = 1:N
            A = zeros(d, n*(N-1));
            s = 1;
            b = (squeeze(path(i, :, t+1) - path(i, :, t))/dt)';
            for k = 1:n
                for j = 1:N
                    if j ~= i
                        dif = squeeze(path(i, :, t) - path(j, :, t))';
                        dist = norm(dif);
                        A(:, s) = I.dict{k}(dist).*(dif/dist);
                        s = s+1;
                    end
                end
            end
            AA_all_particle = [AA_all_particle;A'*A];
            Ab_all_particle = [Ab_all_particle;A'*b];
        end
        
        AA = AA + AA_all_particle;
        Ab = Ab + Ab_all_particle;
    end
end

AA = AA/M; Ab = Ab/M; 

%% Solve linear system and get graph and coefficients

E = zeros(N, N);
all_coef = zeros(n, N);

% get the Z arrays
% Z is the product of graph and coef
Z_arrays = zeros(N-1, n, N);
A_arrays = zeros((N-1)*n, (N-1)*n, N);

condA_orals = zeros(1,N);  condL = condA_orals;
for i = 1:N
    %extract the matrix
    AA_i = AA((i-1)*(N-1)*n+1:i*(N-1)*n, :);
    Ab_i = Ab((i-1)*(N-1)*n+1:i*(N-1)*n, :);
    % solve linear system

    switch reg_method
        case 'None'
            Z = AA_i\Ab_i;
        case 'pinv'
            Z = pinv(AA_i)*all_b;
        case 'lsqminnorm'
            Z = lsqminnorm(AA_i, Ab_i);
        case 'ID'
            [~, Z] = L_curve_standard_form(AA_i, Ab_i, 0);  % to supply the Basis matrix
        case 'RKHS'
            [~, Z] = L_curve(AA_i, Ab_i, 'RKHS', 0, I.ORALS_mat);
    end

    if testON==1
        condA_orals(i) = cond(AA_i);
        condL(i) = cond(I.ORALS_mat\AA_i);
    end
    % reshape
    A_arrays(:,:,i) = AA_i;
    Z_arrays(:,:,i) = reshape(Z, [N-1, n]);
end
if testON==1
    fprintf('regu Method = %s, Average condA = 10^%2.0f, condL = 10^%2.0f   ...',reg_method, log10(mean(condA_orals)),log10(mean(condL)));
else
    fprintf('regu Method = %s,   ...',reg_method);
end

%% SVD for Z
for i = 1:N
    Z = Z_arrays(:,:,i);
    [U, S, V] = svd(Z);
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
E = E';
c = mean(all_coef, 2);

% How to use the property that E is positive?

%% ALS starting from SVD result
niter = 5;
[graph2,  coef2] = ALS_inORALS(Z_arrays, c, niter, I, plotON);

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



fprintf('done (%.2f sec).\n',toc);
time = toc; 
end


function [E, c] = ALS_inORALS(Z_arrays, c0, niter, I, plotON)
% the alternating least squares part of ORALS:
% estimate array a and vector c, with c shared between rows of a;


% I tried this code with obs_std = 0. It does not coverge to the true
% value. I will check the details.

% Numerical results shows that, ORALS is always worse than ORSVD.

[n, ~, N] = size(Z_arrays);

% err_g = zeros(niter, 1);
err_k = zeros(niter+1, 1);
err_k(1) = kernel_err(c0, I);
for t = 1:niter
    E = zeros(n, N);
    c = zeros(size(c0));
    for i = 1:N
        Z = Z_arrays(:, :, i);
        u_half = Z*c0;
        u1 = u_half/norm(u_half);
        E(:,i) = u1;
        %        v1 = Z'*u0;  coef0 = v1; % update coef0 in each row---
        c = c + Z'*u1;
    end
    c = c/N;
    err_k(t+1) = kernel_err(c, I);
    c0 = c;
end

%%
if plotON
    figure;
    plot(err_k);
    title('The change of the kernel error with the iteration of ORALS')
end
end

