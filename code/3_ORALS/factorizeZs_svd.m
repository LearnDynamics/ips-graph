
function [E,c] = factorizeZs_svd(N,n,all_Z,normalizeON)
%% Matrix factorization to get A, c from Z_i's by SVD
E         = zeros(N, N);
all_coef  = zeros(n, N);

for i = 1:N
    Z_i       = all_Z(:,:,i);
    [U, S, V] = svds(Z_i,1);   % only the first singular value 
    graph     = U(:, 1);
    coef      = V(:, 1)*S(1,1);
    % [U, S, V] = svd(Z_i);  graph = U(:, 1);  coef = V(:, 1)*S(1,1);
    
    % make the graph to be positive
    [~, ind] = max(abs(graph));     sgn = sign(graph(ind));
    graph    = graph/sgn; 
    coef     = coef/sgn;

    E(i, :)  = [graph(1:(i-1));0;graph(i:end)];      % store the output;  E= E' below
    all_coef(:, i) = coef;
end

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