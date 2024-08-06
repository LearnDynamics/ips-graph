%% test
%{
Both graph and kernel can be accurately estimated. 
- Interesting observations: 
  + graph error decays as n increases;     
%}
N=6; n=20; std_noise = 0.01; 
a_true  = set_graph(N, 'sparsity', 0.6, 'plotON', 0);                                                   % create influence graph
% a_true = rand(N,N); a_true= a_true - diag(diag(a_true)); 

c_true = randn(n,1);  % [1,2,3]';%

[all_Z,a_true] = Ac_toZ_normalizeA(a_true,c_true); 

all_Z = all_Z+ std_noise*randn(size(all_Z));
normalizeON = 1; 

%% 
rng(12);
c0 =1+ randn(n,1); 
niter = 10; 
 [Emat, c,E1_seq,cseq] = ALS_inORALS(all_Z, c0, niter, normalizeON);
 error_Emat = zeros(1,niter);
 error_c = zeros(1,niter);
 for t=1:niter
    error_Emat(t) = norm(E1_seq(:,:,t)-a_true,'fro');
    error_c(t)    = norm(cseq(:,t)-c_true);
 end
 figure; plot(1:niter,error_Emat,1:niter,error_c,'linewidth',2);
legend('Graph error','Coef error');
xlabel('Iteration number')
c 
% E0 = randn(N,N);
% [Emat, c,E1_seq,cseq] = ALS_inORALS_initial_A(all_Z, E0, niter, normalizeON);
% norm(cseq-c_true)


%% the SVD as start point does not affect the result much; 
[Esvd,csvd] = factorizeZs_svd(N,n,all_Z,normalizeON); 
[Emat1, c1,E1_seq,cseq] = ALS_inORALS(all_Z, csvd, niter, normalizeON);  
[Emat2, c2,E1_seq,cseq] = ALS_inORALS_initial_A(all_Z, Esvd', niter, normalizeON);

error_Emat =[ norm(Esvd'-a_true,'fro'), norm(Emat1-a_true,'fro'), norm(Emat2-a_true,'fro')]
error_c=[ norm(csvd-c_true), norm(c1-c_true), norm(c2-c_true)]


function [all_Z,a_true] = Ac_toZ_normalizeA(a_true,c_true)
% get Z from a and c 
% can also use function get_Z_from_E_c.m; but here we normalize A rows
N = length(a_true(:,1));
n = length(c_true);
all_Z  = zeros(N-1,n,N);
for i=1:N
    a_temp        = a_true(i,:);
    a_temp        = a_temp/norm(a_temp);
    a_true(i,:)   = a_temp;
    ind1N_i       = 1:N; ind1N_i(i) =[];
    a_remove_i    = a_temp(ind1N_i);
    all_Z(:,:,i)  = a_remove_i'*c_true';
end
end