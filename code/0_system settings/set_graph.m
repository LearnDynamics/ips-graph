function A = set_graph(N, varargin)
% randomly set a graph

% sparsity of the graph, [0, 1].
% connectivity of each node. Min = 1, Max = N-1
% connectivity ~~ (N-1)*sparsity.
% extrem values are considered so no error message will pop.
% e.g. sparsity = -0.5. That means connectivity = 1.

%% Input parser
p = inputParser;
addRequired(p,'N');
addOptional(p, 'sparsity',1);                                                                                                   % sparsity is the percentage of nonzero entries connected to a node
addOptional(p, 'plotON', 1);
addOptional(p, 'normalized', 1);

parse(p, N, varargin{:});
sparsity = p.Results.sparsity;
normalized = p.Results.normalized;
plotON = p.Results.plotON;

%%
A = zeros(N, N);
for i = 1:N
    r = rand(N-1, 1); 
    % randomly set (1-sparsity)*(N-1) entries to zero
    num = max(min(N-1 - floor(sparsity*(N-1)), N-2), 0);
    r(randperm(N-1, num)) = 0;
    A(:, i) = [r(1:i-1);0;r(i:end)];
end

if normalized
    column_norm = sqrt(sum(A.^2, 1));
    A = A./column_norm;
end

if plotON
    plot_graph(A)
    ttl = ['sparsity = ', num2str(sparsity), ', N = ', num2str(N)];
    title(ttl)
    tightfig;
end

end