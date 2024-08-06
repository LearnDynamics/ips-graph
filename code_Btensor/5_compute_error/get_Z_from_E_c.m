function Z = get_Z_from_E_c(E, c)

% reconstruct the tensor Z from E and c to compare ORSVD and ORALS

N = length(E);
n = length(c);
Z = zeros(N-1, n, N);
for j = 1:N
    temp        = E(:, j);
    temp(j)     = [];
    Z(:, :, j)  = temp*c';
end

end