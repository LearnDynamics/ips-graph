function [lambda_opt, c] = L_curve(A, b, method, plotON, B, x)
%{
In this code, we wish to minimize

        c'Ac - 2b'c + b'A^{-1}b + lambda*(c-x)'B^{-1}(c-x)

Where the term (c-x)'B^{-1}(c-x) could be thought as the Gaussian prior
N(B, x). We define the following standard form
        c'Ac - 2b'c + b'A^{-1}b + lambda c'c
This problem could be solved using the code
        L_curve_standard_form(A, b, plotON)
We could then transfer the current problem to the standard form, simply
        AA = sqrt(B)*A*sqrt(B)
        bb = sqrt(B)*b - sqrt(B)*A*x
        cc = L_curve_standard_form(AA, bb, plotON)
        c = sqrt(B)*cc + x
        
Notice that the input matrix B is the Gaussian prior covariance

***************************************************************************
B SHOULD BE THE INVERSE OF THE REGULARIZATION MATRIX OF THE LOSS FUNCTIONAL
***************************************************************************

Different choices of methods corresponding to the following settings
    'ID'        B = Id,         x = 0 (No need to specify B and x)
    'RKHS'      B = A,          x = 0 (No need to specify B and x)
                *If B is given, B is taken to be the measure rho matrix Rho == should it be the basis matrix in L2rho? or in RKHS? 
                *Rho matrix must have full rank
    'custom'    B and x are given by the user
                ** If B is not given, B is set to be Id
                ** If x is not given, x is set to be zeros

[lambda_opt, c] = L_curve(A, b, 'ID,        1)
[lambda_opt, c] = L_curve(A, b, 'RKHS',     1)
[lambda_opt, c] = L_curve(A, b, 'RKHS',     1, B)
[lambda_opt, c] = L_curve(A, b, 'custom',   1)
[lambda_opt, c] = L_curve(A, b, 'custom', 	1, B)
[lambda_opt, c] = L_curve(A, b, 'custom', 	1, B, x)

2022-10-03 Quanjun Lang
%}
switch method
    case 'ID' 
        [lambda_opt, c] = L_curve_standard_form(A, b, plotON);
    case 'RKHS'
        if nargin == 4
            [U, S, ~] = svd(A);
            sqrt_A = U*sqrt(S)*U';
            bb = sqrt_A*b;
            AA = A*A;
            [lambda_opt, cc] = L_curve_standard_form(AA, bb, plotON);
            c = sqrt_A*cc;
        elseif nargin >= 5
                assert(rank(B) == length(B), 'B must have full rank');      % TODO: treat this case by chol, similar to below, but avoid \ 
                % Find the corresponding regularization matrix
                % Using svd to compute the generalized eigenvalue problem
                % Increase stability
                % Get BB = Binv*A*Binv by svd instead pinv, to increase stability 
                L0 = chol(B);          % B = L'*L, L = sqrtB 
                AA = (L0'\A)*inv(L0);   % AA = sqrtBinv * A * sqrtBinv
                [U, S, ~] = svd(AA);  % AA = U*S*U'
                V = L0\U;              % V = sqrtBinv*U
                BB = V*S*V';          % BB = sqrtBinv*U*S*U'*sqrtBinv = Binv*A*Binv

                % Use the following customized method
                [U, S, ~] = svd(BB);
                s = diag(S);
                L = U*diag(sqrt(s))*U';
                AA = L'*A*L;
                bb = L'*b;
                [lambda_opt, cc] = L_curve_standard_form(AA, bb, plotON);
                
                c = L*cc;                            
        end
    case 'custom'
        if nargin == 4
            B = eye(n);
            x = zeros(n, 1);
        elseif nargin == 5
            x = zeros(n, 1);
        end
        [U, S, ~] = svd(B);
        s = diag(S);
        L = U*diag(sqrt(s))*U';
        AA = L'*A*L;
        bb = L'*b - L'*A*x;
        [lambda_opt, cc] = L_curve_standard_form(AA, bb, plotON);
        
        c = L*cc + x;
end
if plotON
    sgtitle(method)
end
end
