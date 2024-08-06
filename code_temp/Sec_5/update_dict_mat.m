function learning_setup = update_dict_mat( learning_setup, paths, option )

if ~exist('option', 'var')
    option = 'rho';
end

%% Estimate rho from data
M           = length(paths);
[N, ~, L]   = size(paths{1});

rho_sample  = cell(M,1);%zeros(M*N*(N-1)/2*(L-1),1);
parfor m = 1:M
    X       = paths{m};
    idx     = 1;
    for i = 1:N
        Xi  = X(i, :, :);
        dif = Xi(:, :, 1:end-1) - X(1:i-1, :, 1:end-1);
        dist = sqrt(sum(dif.^2, 2));
        dist = dist(dist~=0);
        rho_sample{m}(idx:idx+length(dist)-1) = dist;
        idx = idx+length(dist);
    end
end

rho_sample = [rho_sample{:}];

%% Compute histogram of rho and update dictionary matrix
switch option
    case 'rho'
        f = figure; set(f,'visible','off');
        [rho,bin_edges] = histcounts(rho_sample,'Normalization','pdf');
        bin_mids = 0.5*(bin_edges(1:end-1)+bin_edges(2:end));

        dx = bin_edges(2)-bin_edges(1);

        n = learning_setup.n;
        dict_mat = zeros(n, n);
        for i = 1:n
            for j = i:n
                dict_mat(i,j) = sum(learning_setup.dict{i}(bin_mids).*learning_setup.dict{j}(bin_mids).*rho)*dx;
                if j>i, dict_mat(j,i) = dict_mat(i,j); end
            end
        end
        learning_setup.rho           = rho;
        learning_setup.rho_dx        = dx;
        learning_setup.rho_bin_edges = bin_edges;
    case 'rho_support'
        max_supp_rho = max(rho_sample);

        n = learning_setup.n;
        dict_mat = zeros(n, n);
        for i = 1:n
            for j = i:n
                dict_mat(i, j) = integral( @(x) learning_setup.dict{i}(x).*learning_setup.dict{j}(x), 0, max_supp_rho );
                if j>i, dict_mat(j,i) = dict_mat(i,j); end
            end
        end
end

learning_setup.dict_mat  = dict_mat;

% fprintf('Cond(dict-mat) = 10^%2.0f, cond(ORALS dict-mat) = 10^%2.0f   ... \n ', log10(cond(dict_mat)),log10(cond(I.ORALS_mat)));

end