%%draw_graph_kernel_M_ALS
function draw_graph_kernel_M_ALS(Mseq, ind, I, error, regustr, ttl)
[~,serr] = size(error);
num0 = (I.N-1+I.n)/I.d/I.steps;
num1 = (I.N-1)*I.n/I.d/I.steps;

figure; % plot error of ORALS in M



for ii = 1:serr
    %     ii=1;
    loglog(Mseq(ind),error{ii}.g_als(ind,1),'-.x','linewidth',1.5); hold on;
    loglog(Mseq(ind),error{ii}.k_als(ind,1),'-.o','linewidth',1.5); hold on;
    grid on
    %     label1=num2str(regustr{serr});
    legendInfo{2*ii-1}=['Graph: ALS+' num2str(regustr{ii})];
    legendInfo{2*ii}=['Kernel: ALS+' num2str(regustr{ii})];
    %       legend('Graph: ALS+','Kernel: ALS+');
    %     legend('location','best');
end


xline(num1);
text(num1,I.obs_std/1.3,'n_1')
xline(num0);
xlim([20, inf]);
text(num0,I.obs_std/1.3,'n_0')
yline(I.obs_std);
text(min(Mseq(ind)),I.obs_std*1.3,'obs std')
ylim([I.obs_std/50, 0.1])
xlabel('Sample size M'); ylabel('Relative L2 Error');
legend(legendInfo);
title(ttl)
end


%% Previous code
% loglog(Mseq(ind),error.g_orals(ind,1),'-x','linewidth',1); hold on;
% loglog(Mseq(ind),error.k_orals(ind,1),'-o','linewidth',1); hold on;
% loglog(Mseq(ind),error.g_orsvd(ind,1),'--x','linewidth',1.5); hold on;
% loglog(Mseq(ind),error.k_orsvd(ind,1),'--o','linewidth',1.5); hold on;