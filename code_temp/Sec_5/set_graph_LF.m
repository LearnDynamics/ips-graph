function E = set_graph_LF(N, varargin)
% randomly set a graph
% sparsity of the graph, [0, 1].
% connectivity of each node. Min = 1, Max = N-1
% connectivity ~~ (N-1)*sparsity.
% extrem values are considered so no error message will pop.
% e.g. sparsity = -0.5. That means connectivity = 1.

%% Input parser
p = inputParser;
addRequired(p, 'N');
addOptional(p, 'leader_ratio', 0.2);   % sparsity is the percentage of nonzero entries connected to a node
addOptional(p, 'leader_affect_ratio', 0.3);
addOptional(p, 'follower_affect_ratio', 0.1);
addOptional(p, 'plotON', 1);
addOptional(p, 'normalizeON', 1);

parse(p, N, varargin{:});
leader_ratio = p.Results.leader_ratio;
leader_affect_ratio = p.Results.leader_affect_ratio;
follower_affect_ratio = p.Results.follower_affect_ratio;

normalizeON = p.Results.normalizeON;
plotON = p.Results.plotON;
%%
E = zeros(N, N);
% randomly select the index of leaders
leader_num = floor(N*leader_ratio);
leader_id = randsample(N, leader_num);

for i = 1:N
    % randomly set (1-ratio)*(N-1) entries to zero
    if ismember(i, leader_id)
        num = max(min(N-1 - floor(leader_affect_ratio*(N-1)), N-2), 0);
    else
        num = max(min(N-1 - floor(follower_affect_ratio*(N-1)), N-2), 0);
    end
    r = rand(N-1, 1);
    r(randperm(N-1, num)) = 0;
    E(:, i) = [r(1:i-1);0;r(i:end)];
end

E = E';

% if normalized
%     mat_norm = norm(E, 'fro');
%     E = E ./ mat_norm;
% end

% if normalized
%     column_norm = sqrt(sum(E.^2, 1));
%     E = E./column_norm;
% end

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

end
