function [Emat, c,E1_seq,cseq] = ALS_inORALS_initial_A(all_Z, E0, niter, normalizeON)
% The alternating least squares part of ORALS: estimate array a and vector c
% Difference from SVD: c shared between rows of A_adj;
%       all_Z size: N-1 x n x  N          

[~,n , N] = size(all_Z);

% estimate c0 from E0
E0 = E0';   % NOTE: E0 is the matrix a in a*B*c
c0 = zeros(n,1); 
for i=1:N
      Z  = all_Z(:, :, i); 
      ind1N_i = 1:N; ind1N_i(i) =[];
      ui = E0(ind1N_i,i); 
      c0 =c0+ Z'*ui;    
end
c0 = c0/N; 

% iteration starting from c0
[Emat, c,E1_seq,cseq] = ALS_inORALS(all_Z, c0, niter, normalizeON);
end

