function err = get_Z_error( Z, dyn_sys, learning_setup )

%% Not using the Z matrix
% err = norm(reshape(Z - I.Z_true, [], 1))./norm(reshape(I.Z_true, [], 1));

learning_setup.ORALS_mat = kron(learning_setup.dict_mat, eye(dyn_sys.N-1));                                                                                     

%% Using the Z matrix get_Z_error(Z_OR, I)
[~, ~, N] = size(Z);

err = 0;
for i = 1:N
    Zi = Z(:, :, i);
    Zi_vec = reshape(Zi, [], 1);
    
    Zi_true = learning_setup.Z_true(:, :, i);
    Zi_true_vec = reshape(Zi_true, [], 1);
    
    Zi_true_norm = sqrt(Zi_true_vec'*learning_setup.ORALS_mat*Zi_true_vec);

    Zi_loss = sqrt((Zi_vec - Zi_true_vec)'*learning_setup.ORALS_mat*(Zi_vec - Zi_true_vec));
    err_i = Zi_loss/Zi_true_norm;
    
    err = err + err_i;
end

err = err/N;

end
