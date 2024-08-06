%% test computational complexity of least squares: 
% order in n for size(A)= [m,n] 
%    lsqmininorm vs   normal matrix + backslash  
%   Result:  both O(n^2), but the later is faster  


close all; clear all; 

m     = 100; 

%%% test in n,m: 
nsimu = 10;
mseq  = 100*4.^(0:4);
nseq  = 10*2.^(1:7); 
time_seq = zeros(nsimu,2,length(nseq),length(mseq)); 
for k = 1:length(nseq)
    for j=1:length(mseq)
    time_seq(:,:,k,j) = timeLS(mseq(j),nseq(k),nsimu);
    end
end

tseq_mean = squeeze( mean(time_seq,1)); 
% tseq_std = std(time_seq,2,1);
%% 
figure;  
subplot(121)% order in n 
tseq_mean_n = tseq_mean(:,:,end);  
loglog(nseq,tseq_mean_n(1,:),'-x');hold on; 
loglog(nseq,tseq_mean_n(2,:),'-o');hold on; 
% loglog(nseq,(nseq/nseq(1)).^3*tseq_mean1(1,1) ,'--'); 
loglog(nseq,(nseq/nseq(1)).^2*tseq_mean_n(1,1) ,'-.'); 
legend('Lsqmininorm','normal Mat','O(n^2)'); 
xlabel('n'); ylabel('Time (seconds)')
% figname = 'lse_comput_time_n';
% tightfig; 
% print([figname,num2str(mseq(end)),'.pdf'],'-dpdf');  


%% 
subplot(122)% order in m
tseq_mean_m = tseq_mean(:,end,:);  
loglog(mseq,tseq_mean_m(1,:),'-x');hold on; 
loglog(mseq,tseq_mean_m(2,:),'-o');hold on; 
loglog(mseq,(mseq/mseq(1))*tseq_mean_m(1,1) ,'-.'); 
legend('Lsqmininorm','normal Mat','O(m)'); 
xlabel('m'); ylabel('Time (seconds)');
% figname = 'lse_comput_time_m';
% tightfig; 
% print([figname,num2str(nseq(end)),'.pdf'],'-dpdf');  
figname = 'lse_comput_time';
tightfig; 
print([figname,'.pdf'],'-dpdf');  


function t_n = timeLS(m,n,nsimu)
t_n   = zeros(nsimu,2); % lsqmininorm and LSE-using normal matrix
for i=1:nsimu
    A = randn(m,n);
    b = rand(m,n);

    % solve using lsqmininorm 
    tic
    x = lsqminnorm(A,b);
    t_n(i,1) = toc;

    % solve using normal matrix and A\b
    tic
    Abar = A'*A; 
    bbar = A'*b; 
    x = Abar\bbar; 
    t_n(i,2) = toc;
end
end


function t_n =time_Abackslash(n,nsimu)
t_n   = zeros(nsimu,1); % lsqmininorm and LSE-using normal matrix
for i=1:nsimu
    A = randn(n,n);
    b = rand(n);

    % solve using normal matrix and A\b
    tic
    x = A\b; 
    t_n(i,2) = toc;
end
end

