function learning_setup = update_dict_mat( learning_setup, paths, option, new_basis_flag)

if ~exist('option', 'var') || isempty(option)
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

if ~exist('new_basis_flag')
    change_basis_fn = 0;
else
    change_basis_fn = 1;  % not used yet, as there will be many changes later: the coeficient c, and Z matrix
    % It is good for non-parametric regression with ill-conditioned/deficient ranked matrices.
end
learning_setup.basis_changed    = change_basis_fn;
if change_basis_fn == 1
    [dictFn,dict_mat,n_new,transfrom_U,flag] = update_basisFn(learning_setup.dict,dict_mat);
    if flag ==1   %otherwise, no need of any changes
        learning_setup.dict_orig      = learning_setup.dict;
        learning_setup.dict_mat_orig  = learning_setup.dict_mat;


        learning_setup.n                = n_new;
        learning_setup.dict             = dictFn;
        learning_setup.transfromU_for_c = transfrom_U;  % % c_original = U*c_new;
        learning_setup.dict_mat         = dict_mat;

        if isfield(learning_setup, 'c')
            learning_setup.c_orig         = learning_setup.c;
            learning_setup.c                = transfrom_U*learning_setup.c_orig;  % % c_original = U*c_new;
        end
    end
    learning_setup.basis_changed    = flag;
end

% fprintf('Cond(dict-mat) = 10^%2.0f, cond(ORALS dict-mat) = 10^%2.0f   ... \n ', log10(cond(dict_mat)),log10(cond(I.ORALS_mat)));

end


function [dictFn,B,n_new, U,flag] = update_basisFn(dictFn,B)
% update the basis functions if the basis matrix in L2(rho) is near singular 
% c_original = U*c_new

if rank(B) == length(dictFn)
    n_new = [];% length(dictFn); 
    U = [];% eye(n_new); 
    flag = 0;
    return;               % if dict_mat is full rank, return
end
flag = 1; 
fprintf('Basis matrix not full rank! Change basis functions keeping eig >1e-12. \n'); 

% [U,S,V] = svd(Bmat);  % 

[U,S]   =eig(B); % since Bmat is real symmetric, eig returns U being orthonormal, B*U = U*S; 
[eigB_sort ,indx]   = sort(diag(S),'descend');




% n_new   = find(eigB_sort>1e-12);      % QL: This is problematic. function  find provides an array, 
n_new   = find(eigB_sort>1e-12, 1, 'last');
n = length(B);



indx    = indx(1:n_new); 
U_sort  = U(:,indx);

dictFn_new = cell(n_new,1);
for i=1:n_new
    dictFn_new{i} = @(x) 0*x;
    for j=1:n
        dictFn_new{i} = @(x) dictFn_new{i}(x) + dictFn{j}(x)*U_sort(j,i);
    end
end
B = diag(eigB_sort(1:n_new));  % B is a diagonal matrix. We can avoid its inversion in DARTR or iDARR

end