function [c, lambda_opt] = L_curve_paper(A, b, method, varargin)
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
                *If B is given, B is taken to be the measure rho matrix Rho
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

%%
n = length(A);
expectedMethods = {'ID','RKHS','custom','true_prior','plain'};

p = inputParser;
addRequired(p,  'A',         @(x) ismatrix(x));
addRequired(p,  'b',         @(x) isvector(x));
addRequired(p,  'method',    @(x) any(validatestring(x,expectedMethods)));
addParameter(p, 'Bmat',     eye(n),      @(x) isnumeric(x));
addParameter(p, 'prior_m',  zeros(n, 1),      @(x) isnumeric(x));
addParameter(p, 'prior_Q',  eye(n),      @(x) isnumeric(x));
addOptional(p, 'sigma',    inf,      @(x) isnumeric(x));
addParameter(p, 'plotON',   false,  @(x) islogical(x));

parse(p, A, b, method,varargin{:});

Bmat = p.Results.Bmat;
plotON = p.Results.plotON;
Q = p.Results.prior_Q;
m = p.Results.prior_m;
sigma = p.Results.sigma;
%%

n = length(A);

switch method
    case 'ID'
        assert(length(A) == length(b), 'The size of A and b do not match');
        [lambda_opt, c] = L_curve_standard_form(A, b, plotON);
        if plotON
            sgtitle(method)
        end
    case 'RKHS'
        use_Binv = 2;    % use Binv as in paper if 1, otherwise, the one without Binv, which is better in general. 
        if use_Binv ==1       
            bb = A*pinv(Bmat)*b;
            AA = (A*pinv(Bmat))^2;
            [lambda_opt, cc] = L_curve_standard_form(AA, bb, plotON);
            c = Bmat\cc;
        else  
            assert(rank(Bmat) == length(Bmat), 'B must have full rank');
            % Find the corresponding regularization matrix
            % Using svd to compute the generalized eigenvalue problem
            % Increase stability
            L = chol(Bmat);
            AA = (L'\A)*inv(L);
            [U, S, ~] = svd(AA);
            V = L\U;
            BB = V*S*V';
            
            % Use the following customized method
            [U, S, ~] = svd(BB);
            s = diag(S);
            L = U*diag(sqrt(s))*U';
            AA = L'*A*L;
            bb = L'*b;
            [lambda_opt, cc] = L_curve_standard_form(AA, bb, plotON);           
            c = L*cc;
            
            %         [U, S, ~] = svd(A);
            %         sqrt_A = U*sqrt(S)*U';
            %         bb = sqrt_A*b;
            %         AA = A*A;
            %         [lambda_opt, cc] = L_curve_standard_form(AA, bb, plotON);
            %         c = sqrt_A*cc;
        end
        if plotON;sgtitle(method);end
    case 'custom'
        [U, S, ~] = svd(Q);
        s = diag(S);
        L = U*diag(sqrt(s))*U';     %L is square root of B
        AA = L'*A*L;
        bb = L'*b - L'*A*m;
        [lambda_opt, cc] = L_curve_standard_form(AA, bb, plotON);
        if plotON
            sgtitle(method)
        end
        c = L*cc + m;
    case 'true_prior'
        assert(isfinite(sigma), 'Must specify the noise level if using the true prior');
        c = (A + sigma^2*Q)\b;
        lambda_opt = sigma^2;
    case 'plain'
        c = A\b;
        lambda_opt = sigma^2;
end

% if isfinite(sigma)
%     lambda_opt = lambda_opt/(sigma^2);
% end


end
