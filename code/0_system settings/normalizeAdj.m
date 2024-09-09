function [A,cols_norm] = normalizeAdj ( A, normalizationType, A_sparsity_thres, A_sparsity_prior )

% Normalization of columns of adjacency matrix

N = size(A,1);
if nargin>3 && A_sparsity_prior<1
    for i = 1:N                                                                                                             % set (1-sparsity)*(N-1) entries to zero in each row; pick the smallest entries
        [~,entries_idx] = sort(A(i,:),'ascend');
        A(i,entries_idx(1:floor(N-1-A_sparsity_prior*(N-1)))) = 0;
    end
end
if nargin>2 && A_sparsity_thres>0
    A = A.*( abs(A)>A_sparsity_thres );
end
switch normalizationType
    case 1
        cols_norm = sqrt(sum(A.^2, 1));
        cols_norm(cols_norm<eps) = 1;
        A = A./cols_norm;
    case 'mat'
        mat_norm = norm(A, 'fro');
        A = A ./ mat_norm;
end

end