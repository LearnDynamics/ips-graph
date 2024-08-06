function G = createGraph( graphType, opts )

%
% function G = createGraph( graphType, opts )
%
% IN:
%   graphType  : type of graph to be built. Any of the following:
%                   'Tree'      : Tree graph
%                             opts must contain:
%                               K       : number of branches at each node
%                               Levels  : number of levels
%                             The suggested embedding vX has the root at
%                             (0,0).
%                   'Complete'  : complete graph
%                               N       : number of vertices
%                   'Circle'    : circle graph
%                               N       : number of vertices
%                   'random_sparse' : random sparse graph
%                               N       : number of vertices
%                             opts must contain:
%                               sparsity    : only sparsity*(N-1) non-zero edges per vertex
%                               normalized  : true/false for l^2 normalization of the columns of adjacency matrix
%
%
%   opts       : structure of options, some of which specific (see cGraphType), some generic:
%                   [Undirected]    : 0,1. Default: 1
%                   [Embed]         : 0,1. Default: 0
%                   [selfLoops]     : 0,1. representing 'no self loops','self loops as in the model'. Default: 1
%
%
% OUT:
%   G.W          : adjacency matrix
%   G.X          : D by N matrix of D suggested coordinates for embedding
%                   the graph.
%   G.D          : N column vector of outdegrees
%
% EXAMPLE
% G = CreateGraph('Tree',struct('K',3,'Levels',4,'Embed',true) );
% DrawSimpleGraph(G);
% BGobj=biograph(G.W);BGobj.view;
%

% (c) Mauro Maggioni

G.W = [];
G.X = [];

if ~isfield(opts,'Undirected'), opts.Undirected = 1;    end
if ~isfield(opts,'Embed'),      opts.Embed = 0;         end
if ~isfield(opts,'selfLoops'),  opts.selfLoops = 1;     end

switch lower(graphType)
    case {'tree'}
        idx = 1;
        % Allocate memory
        I((opts.K^(opts.Levels+1)-1)/(opts.K-1)-1) = 0;
        J((opts.K^(opts.Levels+1)-1)/(opts.K-1)-1) = 0;
        if opts.Embed
            G.X(2,(opts.K^(opts.Levels+1)-1)/(opts.K-1)-1) = 0;
        end
        if opts.Embed
            G.X(:,1) = [0;0];
        end
        % Run through the levels
        for j = 1:opts.Levels
            % Connect the vertices at level J-1...
            I(idx:(idx+opts.K^j-1)) = floor((idx-1+(0:opts.K^j-1))/opts.K)+1;
            % ...with the new leaves at level J
            J(idx:(idx+opts.K^j-1)) = idx+(0:(opts.K^j-1))+1;
            % Construct embedding of the graph if requested
            if opts.Embed
                G.X(1,(idx:(idx+opts.K^j-1))+1) = idx+(0:(opts.K^j-1))+1;
                G.X(1,(idx:(idx+opts.K^j-1))+1) = G.X(1,idx:(idx+opts.K^j-1))-mean(G.X(1,idx:(idx+opts.K^j-1)));
                G.X(2,(idx:(idx+opts.K^j-1))+1) = -j;
            end
            idx = idx + opts.K^j;
        end
        % Assemble the sparse adjacency matrix
        G.W = sparse(I,J,1,idx,idx);
    case {'complete'}
        G.W = ones(opts.N);
        if opts.Embed
            % For lack of better ideas, let's do a 2D layout with points around
            % a circle
            lthetas = linspace(0,2*pi*(1-1/opts.N),opts.N);
            G.X = [cos(lthetas);sin(lthetas)];
        end
    case {'circle'}
        e               = 0.25*ones(opts.N,1);
        G.W             = spdiags([e 2*e e],-1:1,opts.N,opts.N);
        G.W(1,opts.N)   = 0.25;
        G.W(opts.N,1)   = 0.25;
        if opts.Embed            
            lthetas = linspace(0,2*pi*(1-1/opts.N),opts.N);
            G.X = [cos(lthetas);sin(lthetas)];
        end
    case {'random_sparse'}
        G.W = zeros(opts.N, opts.N);
        for i = 1:opts.N
            r   = rand(opts.N-1, 1);            
            num = max(min(opts.N-1 - floor(opts.sparsity*(opts.N-1)), opts.N-2), 0);                                            % randomly set (1-sparsity)*(opts.N-1) entries to zero
            r(randperm(opts.N-1, num)) = 0;
            G.W(:, i) = [r(1:i-1);0;r(i:end)];
        end
        
    otherwise
        fprintf('\n CreateGraph: Unknown graph type.\n');
        return;
end

G.D = sum(G.W,2);

% Rescale all the coordinates so that the graph embedding sits in the unit square
if opts.Embed
    for k = 1:size(G.X,1)
        G.X(k,:) = G.X(k,:)/max(abs(G.X(k,:)));
    end
end

% Undirect the edges if requested to do so
if opts.Undirected
    G.W = G.W+G.W';
end

% Add self-loops
if ~opts.selfLoops
    G.W = G.W-spdiags(G.D,0,size(G.W,1),size(G.W,2));
end

% Normalizes columns if requested
if opts.normalized
    column_norm = sqrt(sum(G.W.^2, 1));
    G.W = G.W./column_norm;
end

return