function err = kernel_multi_err(coef_mat_est, coef_mat_true,  learning_set)

[Q, N] = size(coef_mat_est);
err = 0;
for i = 1:N
    c_est   = coef_mat_est(:, i);
    c_true  = coef_mat_true(:, i);
    c_diff  = c_true - c_est;

    L2_err      = sqrt(c_diff'*learning_set.dict_mat*c_diff);
    L2_err_rel  = L2_err/sqrt(c_true'*learning_set.dict_mat*c_true);
    
    err = err + L2_err_rel;
end

err = err/N;
end




