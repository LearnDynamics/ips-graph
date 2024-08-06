function multi_kernel = get_multi_type_kernel_from_id(coef_mat, id, Q, dict)

multi_kernel = {};

for i = 1:Q
    temp = find(id == i, 1);
    if length(temp) == 1
        multi_kernel{i} = get_kernel_from_c(coef_mat(:, temp), dict);
    end
end


end