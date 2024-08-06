function plot_computTime(comput_time,Mseq,figname)
timeORALS = comput_time(1,:,1);
timeALS   = comput_time(2,:,1);
figure; 
loglog(Mseq,timeALS,'-', 'LineWidth', 2);  hold on; 
loglog(Mseq,timeORALS,'-.', 'LineWidth', 2); 
legend('ALS','ORALS'); 
xlabel('Sample size: M'); ylabel('Computational Time (seconds)')
figname = strrep(figname,'conv_in_M','time');
set_positionFontsAll;
end