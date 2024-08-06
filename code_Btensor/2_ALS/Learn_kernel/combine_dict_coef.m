function phi_kernel = combine_dict_coef( learning_setup,c )

if nargin<2,    c = learning_setup.c_true;   end

phi_kernel  = @(x) 0*x;
for i = 1:learning_setup.n
    phi_kernel = @(x) phi_kernel(x) + c(i)*learning_setup.dict{i}(x);
end

end