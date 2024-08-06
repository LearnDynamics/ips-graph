function [E, c] = learn_kernel_graph_ORALS(M, I)
% joint learning of kernel and graph with M trajectories  

progressON = 0;
rk4ON = 0;
n = I.n;
N = I.N;
d = I.d;
steps = I.steps;
initial = I.initial;
dt = I.dt;
% Gamma = zeros(N, N, n, steps, d);

AA = zeros(N*n*(N-1), n*(N-1));     % Store the matrix A for N particles
Ab = zeros(N*n*(N-1), 1);           % Store the vector b for N particles
% This is because the graph has N columns, and each time we only estimate
% one column of it. But in the parfor loop, we must extract information for
% all columns.


fprintf('\nJoint learning using ORALS M= %i ...',M);tic
for m = 1:M
    X0 = set_particle_initial_all_dim(N, d, initial);
    xpath        = graph_forward_model(I, X0, progressON, rk4ON);      % using RK4, the result will not converge.
    path = xpath + I.obs_std* randn(size(xpath));
    
    % for the i-th particle
    %     for i = 1:N
    for t = 1:steps
        AA_all_particle = [];
        Ab_all_particle = [];
        for i = 1:N
            
            A = zeros(d, n*(N-1));
            s = 1;
            b = (squeeze(path(i, :, t+1) - path(i, :, t))/dt)';
            for k = 1:I.n
                for j = 1:I.N
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


fprintf('done (%.2f sec).\n',toc);

%% Solve linear system and get graph and coefficients

E = zeros(N, N);
all_coef = zeros(n, N);
for i = 1:N
    %extract the matrix
    AA_i = AA((i-1)*(N-1)*n+1:i*(N-1)*n, :);
    Ab_i = Ab((i-1)*(N-1)*n+1:i*(N-1)*n, :);
    
    % solve linear system
    Z = AA_i\Ab_i;
    
    % reshape and svd
    Z = reshape(Z, [N-1, n]);
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


end


