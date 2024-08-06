
function [all_Z,N,n,condA_orals,condL]= operator_regression(B,dX,learning_set,dyn_sys,reg_method,testON)
%  operator regression to estimate all_Z: (N-1, n, N) 
%                Z(:,:,i) = a(:,i)*c,         a(:,i) removes a(i,i)=0  --- not necessary, but may help 
% Regression to estimate z_i= vec(Z(:,:,i))
%  Input
%      B    cell(N,1);   each B{i} is an N x dLM x n array.   dX(i,:) = a(:,i)'*B{i}*c
%     dX    array N x dLM
% Output: 

[N,dLM,n]  = size(B{1}); 
all_Z        = zeros(N-1, n, N); 
% compute the condition number of OR regression
    condA_orals = zeros(1, N);        condL       = zeros(1, N);
    graph_mat   = eye(dyn_sys.N-1);

if strcmp(reg_method,'RKHS') || testON == 1
    ORALS_mat   = kron( learning_set.dict_mat, graph_mat );
end

for i=1:N 
    indx_i    = 1:N;     indx_i(i) = [] ;            % remove index i, since a_ii= 0;   ---- is it necessary? Not
    Btemp     = B{i}(indx_i,:,:);                    % B{i} size = N x dLM x n 
    Btemp     = permute(Btemp,[2,1,3]);              %% XXXX
    A_i       = reshape(Btemp,[ dLM, (N-1)*n]); 
    b_i       = dX(i,:)'; 

    Z_i = 0; 
    % solve linear system
    switch reg_method
        case 'None';          Z_i = A_i\b_i;
        case 'pinv';          Z_i = pinv(A_i)*b_i;
        case 'lsqminnorm';    Z_i = lsqminnorm(A_i, b_i);
        case 'pinvreg'                                       % not a solid method
            pinvregeps = 1e-6; 
            Z_i = pinv(A_i,pinvregeps)*b_i;
        case 'ID' 
            Abar = A_i'*A_i/dLM; bbar =  A_i'*b_i/dLM;        
            [~,Z_i] = L_curve_standard_form(Abar, bbar, 0);  % to supply the Basis matrix
        case 'RKHS'   
            Abar = A_i'*A_i/dLM; bbar =  A_i'*b_i/dLM;   
            [~,Z_i] = L_curve(Abar, bbar, 'RKHS', 0, ORALS_mat);
        case 'RKHS_plain'
            Abar = A_i'*A_i/dLM; bbar =  A_i'*b_i/dLM;
            [~,Z_i] = L_curve(Abar, bbar, 'RKHS', 0);
    end

    if testON==1    % compute the condition number
        Abar = A_i'*A_i/dLM;
        condA_orals(i) = cond(Abar);               warning off
        condL(i)       = cond(ORALS_mat\Abar);     warning on
    end
    all_Z(:,:,i) = reshape(Z_i,[N-1, n]);
end
if testON == 1
    fprintf('regu Method = %s, Average condA = 10^%2.0f, condL = 10^%2.0f   ...\n',reg_method, log10(mean(condA_orals)),log10(mean(condL)));
else
    fprintf('ORALS, regu Method = %s,   ...\n ',reg_method);
end