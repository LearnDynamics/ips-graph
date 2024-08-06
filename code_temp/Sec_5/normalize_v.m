function v_hat = normalize_v (v, K_means_on_v)


if K_means_on_v
    [N, Q] = size(v);
    [ind, x] = kmeans(v, Q);
    v = zeros(N, Q);
    for i = 1:N
        v(i, :) = x(ind(i), :);
    end
    % a = 1;
end

% Normalization of type matrix v
[U, S, V] = svd(v);
v_hat = U*(S>0)*V';



end