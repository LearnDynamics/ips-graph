function err = kernel_err_non_para(c, I)
    f = get_kernel_from_c(c, I.dict);
    phi_est = f(I.rgrid);
    phi_true = I.phi_kernel(I.rgrid);
    
    phi_err = sum((phi_est - phi_true).^2.*I.rho);
    phi_norm = sum(phi_true.^2.*I.rho);
    
    err = phi_err/phi_norm;
end