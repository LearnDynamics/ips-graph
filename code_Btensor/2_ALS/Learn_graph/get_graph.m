function E = get_graph(all_A, all_b)
% Use constrained optimization to estimate the adjacent matrix E
[~, ~, N] = size(all_b);


%% loop
E_temp = zeros(N-1, N);
for i=1:N
    A = all_A(:, :, i);
    b = all_b(:, :, i);
    E_temp(:, i) = get_e(A, b, 'none');
end
E = reform_E(E_temp);
end

function e = get_e(A, b, opt)
N = length(b) +1;
switch opt
    case 'none'
        e = A\b;
    case 'normalized'
        D = [eye(N-2, N-2);-ones(1, N-2)];
        r = [zeros(N-2, 1);1];
        
        AA = D'*A*D;
        bb = D'*A*r + D'*b;
        x = AA\(-bb);
        
        e = r + D*x;
    case 'lsqlin'
        C = chol(A);
        d = -C'\b;
        
        Aeq = ones(1, N-1);
        beq = 1;
        lb = zeros(N-1, 1);
        ub = ones(N-1, 1);
        
        x0 = randn(size(d));
        options = optimoptions('lsqlin','Algorithm','interior-point','Display','off');
        e = lsqlin(C,d,[],[],Aeq,beq,lb, ub, x0, options);
        %         e = lsqlin(C,d,[],[],Aeq,beq,lb, ub);
end
end




