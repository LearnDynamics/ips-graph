function [RIP, COND, time, RIP2, RIP3] = get_RIP_COND(path_data, I,learning_setup)
%
% % joint learning of kernel and graph from data trajectories
% p = inputParser;
% addRequired(p,'path_data');
% addRequired(p,'I');
% addOptional(p, 'plotON', 0);
% parse(p, path_data, I, varargin{:});
% plotON = p.Results.plotON;


M = length(path_data);  % paths_data is a cell array Mx1
n = I.n;
N = I.N;
d = I.d;
steps = I.L;

AA = zeros(N*n*(N-1), n*(N-1));     % Store the matrix A for N particles
% This is because the graph has N columns, and each time we only estimate
% one column of it. But in the parfor loop, we must extract information for
% all columns.
dictmat_largeM = learning_setup.dict_mat;      %  try dict_mat with large M 
learning_setup   = update_dict_mat( learning_setup, path_data);                                          % compute the exploration measure rho and the basis matrix using all paths

dict_matM = learning_setup.dict_mat; 
dict      = learning_setup.dict;

fprintf('\n RIP and COND, M = %i ...', M);tic
parfor m = 1:M
    path = path_data{m};
    for t = 1:steps
        AA_all_particle = [];
        for i = 1:N           
            A = zeros(d, n*(N-1));
            s = 1;
            for k = 1:n
                for j = 1:N
                    if j ~= i
                        dif = squeeze(path(i, :, t) - path(j, :, t))';
                        dist = norm(dif);
                        A(:, s) = dict{k}(dist).*(dif/dist);
                        s = s+1;
                    end
                end
            end
            AA_all_particle = [AA_all_particle;A'*A];
        end
        
        AA = AA + AA_all_particle;
    end
end


%% Get condition number and RIP number

COND_seq = zeros(N, 1);


for i = 1:N
    % This is the matrix for ORALS linear inversion
    AA_i = AA((i-1)*(N-1)*n+1:i*(N-1)*n, :);
    COND_seq(i) = cond(AA_i);
end
COND = mean(COND_seq);
%%


% For the ith particle, 
% we need to combine the M, steps and d dimension
all_AA = cell(N, 1);
for i = 1:N
    AA = zeros(N-1, n, 1);
    s = 1;
    for m = 1:M
        path = path_data{m};
        for t = 1:steps
            AA_temp = zeros(N, n, d);
            for k = 1:n
                for j = 1:N
                    if j ~= i
                        dif = squeeze(path(i, :, t) - path(j, :, t))';
                        dist = norm(dif);
                        AA_temp(j, k, :) = dict{k}(dist).*(dif/dist);
                    end
                end
            end
            AA_temp(i, :, :) = [];
            
            AA(:, :, d*s+2-d:d*s+1) = AA_temp;
            s = s + 1;
        end
    end
    AA(:, :, 1) = [];
    % AA has dimension (N-1)*n*(M*d*steps)
    % AA store the many matrices in the sense of Matrix sensing.
    all_AA{i} = AA;
end

Sample_num = 1000;
[~, ~, MM] = size(all_AA{1});
all_ratio = zeros(N*Sample_num,1);
all_ratio2 = all_ratio;

for i = 1:N
    AA = all_AA{i};
    ratio = zeros(Sample_num, 1);
    ratio2 = ratio;
    parfor s = 1:Sample_num
        g = randn(N-1, 1);
        k = randn(n, 1);
        X = g*k';
        AA_X = squeeze(sum(AA.*X, [1, 2]))/sqrt(MM);
        AA_X_norm = norm(AA_X);
%         X_norm = norm(X, 'fro');
        normg  = norm(g); 
        X_norm = normg * sqrt((k'*dict_matM*k));
        X_norm2 = normg * sqrt((k'*dictmat_largeM*k));
        ratio(s) = (AA_X_norm/X_norm)^2;
        ratio2(s) = (AA_X_norm/X_norm2)^2;
        
        
        
    end
    r = max(ratio)/min(ratio);
    delta = (r-1)/(r+1);
    all_ratio(1+Sample_num*(i-1):Sample_num*i) = ratio;
    all_ratio2(1+Sample_num*(i-1):Sample_num*i) = ratio2;
end



fprintf('done (%.2f sec).\n',toc);
time = toc; 

% confidence interval with 95%
RIP.CI95 = prctile(all_ratio, 97.5)./prctile(all_ratio, 2.5);

RIP.min = min(all_ratio);
RIP.max = max(all_ratio);
r = RIP.max/RIP.min;
RIP.delta = (r-1)/(r+1);
ratio_grid = linspace(floor(log10(RIP.min)), ceil(log10(RIP.max)), 50);
RIP.ratio_val = histcounts(log10(all_ratio), ratio_grid)/(Sample_num*N);
RIP.ratio_grid = ratio_grid(2:end);
% figure;plot(ratio_grid(1:end-1), ratio_val);

RIP2.CI95 = prctile(all_ratio2, 97.5)./prctile(all_ratio2, 2.5);
RIP2.min = min(all_ratio2);
RIP2.max = max(all_ratio2);
r = RIP2.max/RIP2.min;
RIP2.delta = (r-1)/(r+1);
ratio_grid = linspace(floor(log10(RIP2.min)), ceil(log10(RIP2.max)), 50);
RIP2.ratio_val = histcounts(log10(all_ratio2), ratio_grid)/(Sample_num*N);
RIP2.ratio_grid = ratio_grid(2:end);
end