function [err_relL2,err_relLinfy] = kernel_err( c, learning_set )

% Returns the relative L^2 error between the kernel with coefficients c and the kernel with the coefficients learning_setup.c
% (over the dictionary in learning_setup)

% compute relative L^2 error
err_relL2       = sqrt((c - learning_set.c)'*learning_set.dict_mat*(c - learning_set.c));
err_relL2       = err_relL2./sqrt(learning_set.c'*learning_set.dict_mat*learning_set.c+eps);

% compute (coarse approximation to) relative L^\infty error
if nargout>1
    phi_true        = get_kernel_from_c(learning_set.c, learning_set.dict);
    phi_hat         = get_kernel_from_c(learning_set.c, learning_set.dict);
    r               = linspace(0,learning_set.rho_bin_edges(end),1000);
    err_relLinfy    = max( abs(phi_true(r)-phi_hat(r)) );
end

end