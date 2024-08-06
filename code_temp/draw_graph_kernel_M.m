function draw_graph_kernel_M(Mseq, ind, I, error, ttl)
num0 = (I.N-1+I.n); % /I.d/I.steps;
num1 = (I.N-1)*I.n; % /I.d/I.steps;

figure; % plot error of ORALS in M
loglog(Mseq(ind),error.g_orals(ind,1),'-x','linewidth',1); hold on;
loglog(Mseq(ind),error.k_orals(ind,1),'-o','linewidth',1); hold on;
loglog(Mseq(ind),error.g_orsvd(ind,1),'--x','linewidth',1.5); hold on;
loglog(Mseq(ind),error.k_orsvd(ind,1),'--o','linewidth',1.5); hold on;
loglog(Mseq(ind),error.g_als(ind,1),'-.x','linewidth',1.5); hold on;
loglog(Mseq(ind),error.k_als(ind,1),'-.o','linewidth',1.5); hold on;
grid on
xline(num1)
text(num1,I.obs_std/1.3,'n_1')
xline(num0)
text(num0,I.obs_std/1.3,'n_0')
yline(I.obs_std);
text(min(Mseq(ind)),I.obs_std*1.3,'obs std')
ylim([I.obs_std/2, 10])
xlabel('Sample size M'); ylabel('Relative L2 Error');
legend('Graph: ORALS','Kernel: ORALS','Graph: ORSVD','Kernel: ORSVD','Graph: ALS','Kernel: ALS');
legend('location','best');
title(ttl)

end