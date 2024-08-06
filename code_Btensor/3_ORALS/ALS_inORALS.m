
function [Emat, c,E1_seq,cseq] = ALS_inORALS(all_Z, c0, niter, normalizeON)
% The alternating least squares part of ORALS: estimate array a and vector c
% Difference from SVD: c shared between rows of A_adj;
%       all_Z size: N-1 x n x  N   

% QL: I tried this code with obs_std = 0. It does not coverge to the true value. I will check the details.  
% Numerical results shows that, ORALS is always worse than ORSVD.           

[N1,n , N] = size(all_Z);

Emat   = zeros(N,N);
E1_seq = zeros(N,N,niter);   % the sequence of Adj mat in iterations
cseq   = zeros(n,niter);

for t = 1:niter
   E1 = zeros(N1, N); 
   c = zeros(size(c0));  
    switch normalizeON
        case {1, 0}        % normalize by L2  
            Etemp = zeros(N,N);
            for i = 1:N
                Z      = all_Z(:, :, i);
                u_half = Z*c0;
                u_half = lsqnonneg(eye(N1), u_half);          % making the entries non-negative   ******it is helpful.                
                if normalizeON;  u1 = u_half/norm(u_half);    % L2 normalization
                else;            u1 = u_half;
                end
                E1(:,i) = u1;    Etemp(:,i)= [u1(1:(i-1));0;u1(i:end)];
                c = c + Z'*u1;
            end
             if mean(E1,'all')< 0;  c = -c; E1 = -E1;   end

            c = c/N;  cseq(:,t) = c;     E1_seq(:,:,t) =  Etemp';   
        case 'mat'        % normalize by Frobenious norm of the matrix --- problematic. 
            all_c = zeros(length(c0), N); Etemp = zeros(N,N);
            for i = 1:N
                Z       = all_Z(:, :, i); 
                u1      = Z*c0;     
                E1(:,i) = u1;    Etemp(:,i) = [u1(1:(i-1));0;u1(i:end)];
                all_c(:, i) = Z'*u1./(norm(u1).^2);      
            end
            c     = mean(all_c, 2); 
            Enorm = norm(E1, 'fro'); 
            E1    = E1/Enorm;   
            if mean(E1,'all')< 0;  c = -c; E1 = -E1; Etemp= - Etemp;   end
            E1_seq(:,:,t) = Etemp'/Enorm; 
            cseq(:,t) = c;         
    end
    c0 = c;
end

%% fill the diagonal entries by zeros
for i =1:N   
    E_1row     = E1(:,i);    % TBD: it changes from row to column
    Emat(:,i) = [E_1row(1:(i-1));0;E_1row(i:end)];
end
% Emat = Emat';            % transpose here: Because svd gives column vectors, and the graph uses row representation
end


%{
% removed: error using the unknown true kernel/graph 
err_k = zeros(niter+1, 1);
err_k(1) = kernel_err(c0, learning_set);
err_k(t+1) = kernel_err( c, learning_set );
if plotON
    figure;     plot(err_k);
    title('The change of the kernel error with the iteration of ORALS')
end
%}

%{
the case of normalize by Frobenius norm seem problematic --- it uses L2 normalization below.  
      case 'mat'        % normalize by Frobenious norm of the matrix
            for i = 1:N
                Z = all_Z(:, :, i); 
                E1(:,i) = Z*c0;              
            end
            E1 = E1/norm(E1, 'fro');                            

            all_c = zeros(length(c0), N);
            for i = 1:N
                Z  = all_Z(:, :, i); 
                u1 = E1(:,i);
                all_c(:, i) = Z'*u1./(norm(u1).^2);      % --- not really just Fro: the c's are computed with norm(u1)
            end
            c = mean(all_c, 2);    cseq(:,t) = c;   
%}