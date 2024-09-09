function [basis, knots] = spline_basis_MATALB(lb, rb, knot_num, deg, free_bdry)
% Using the default MATLAB function to generate thet B-spline basis
% The result is a function handle, can be used with numerical integrator
% Input: [lb, rb] is the domain for generating the spline functions
%        knot_num: number of knots
%        deg: degree of the polynomials +1 .
%        free_bdry: If true, the boundary value has freedom. Otherwise
%                   boundary value is set to be 0.


%%
if free_bdry
    knots = [lb*ones(1, deg-1), linspace(lb, rb, knot_num+1), rb*ones(1, deg-1)];
    N = knot_num+deg-1;
else
    knots = linspace(lb, rb, knot_num+1);
    N = knot_num-deg+1;
end
basis = cell(N, 1);
for i = 1:N
    coef = zeros(N, 1)';
    coef(i) = 1;
    sp = spmak(knots, coef);
%     fnbrk(sp)
    basis{i} = @(x) fnval(sp, x);
end

end



%% debug
% lb = -4;
% rb = 2;
% knot_num = 10;
% deg = 4;
%
% coef(1) = 1;
% sp = spmak(knots, coef);
% fnbrk(sp)
% g1 = @(x) fnval(sp, x);
%
% coef = zeros(knot_num+deg-1, 1)';
% coef(2) = 1;
% sp = spmak(knots, coef);
% fnbrk(sp)
% g2 = @(x) fnval(sp, x);
%
% fplot(g1, [lb, rb]);hold on;fplot(g2, [lb, rb]);



