function r = get_ratio(c, M, d, initial, dict, dict_mat)
test_kernel = get_kernel_from_c(c, dict);
X0 = set_particle_initial_all_dim(2*M, d, initial);
dif = X0(1:M, :) - X0(M+1:2*M, :);
dis = sqrt(sum(dif.^2, 2));
samples_KX = test_kernel(dis).*dif./dis;
E_KX_2 = sum(samples_KX.^2, 'all')/M;
E_KX = sum(samples_KX, 1)/M;
var_KX = E_KX_2 - sum(E_KX.^2);
norm_kernel_2 = c'*dict_mat*c;
r = var_KX/norm_kernel_2;
end