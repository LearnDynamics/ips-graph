function [E, social_degree] = set_graph_social_range(N, varargin)
% randomly set a graph
% sparsity of the graph, [0, 1].
% connectivity of each node. Min = 1, Max = N-1
% connectivity ~~ (N-1)*sparsity.
% extrem values are considered so no error message will pop.
% e.g. sparsity = -0.5. That means connectivity = 1.

%% Input parser
p = inputParser;
addRequired(p, 'N');
addOptional(p, 'social_ratio', 0.2);   % sparsity is the percentage of nonzero entries connected to a node
addOptional(p, 'social_affect_ratio', 0.3);
addOptional(p, 'associal_affect_ratio', 0.1);
addOptional(p, 'plotON', 1);
addOptional(p, 'normalizeON', 1);

parse(p, N, varargin{:});
social_ratio = p.Results.social_ratio;
social_affect_ratio = p.Results.social_affect_ratio;
associal_affect_ratio = p.Results.associal_affect_ratio;

normalizeON = p.Results.normalizeON;
plotON = p.Results.plotON;
%%
E = zeros(N, N);
% randomly select the index of leaders
social_num = floor(N*social_ratio);
social_id = randsample(N, social_num);

for i = 1:length(social_id)
    id = social_id(i);
    affect_num = max(min(N-1 - floor(social_affect_ratio*(N-1)), N-2), 0);
    r = rand(N-1, 1);
    r(randperm(N-1, affect_num)) = 0;
    E(:, id) = [r(1:id-1);0;r(id:end)];
    
    r = rand(N-1, 1);
    r(randperm(N-1, affect_num)) = 0;
    E(id, :) = [r(1:id-1);0;r(id:end)]';
end

E_associal = set_graph(N, 'sparsity', associal_affect_ratio, 'plotON', 0, 'normalized', 0);

E = E + E_associal';

% deal with the case that E has a zero column
col_norm = sqrt(sum(E.^2, 1));
zero_id = (col_norm == 0);
for i = 1:N
    if zero_id(i)
        if i == N
            E(1, i) = rand();
        else
            E(i+1, i) = rand();
        end
    end
end

% normalization
switch normalizeON
    case 1
        column_norm = sqrt(sum(E.^2, 1));
        E = E./column_norm;
    case 'mat'
        mat_norm = norm(E, 'fro');
        E = E ./ mat_norm;
end


if plotON
    plot_graph(E)
    %     ttl = ['sparsity = ', num2str(sparsity), ', N = ', num2str(N)];
    %     title(ttl)
    tightfig;
end

social_degree = zeros(N, 1);

for i = 1:N
    social_degree(i) = sum(abs(E(i, :))+abs(E(:, i))');
end


end
