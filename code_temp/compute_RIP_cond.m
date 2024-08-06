function [RIP_seq,COND_seq]= compute_RIP_cond(rip_filename,all_path,I,Mseq,num0,num1,plotRIP)
%% compute the RIP and condition number
if ~exist(rip_filename,'file')
    Mseq_length = length(Mseq);
    RIP_seq     = zeros(Mseq_length, 1);
    COND_seq    = zeros(Mseq_length, 1);
    comput_time = zeros(Mseq_length,1);
    for i = 1:Mseq_length
        fprintf('\n M-sequence:  %i out of %i : \n',i, Mseq_length);
        M = Mseq(i);

        [RIP, COND, time1] = get_RIP_COND(all_path(1:M), I);
        RIP_seq(i) = RIP.delta;
        COND_seq(i) = COND;
        comput_time(i) = time1;
    end
    save(rip_filename,'I','comput_time','RIP_seq','COND_seq');
else
    load(rip_filename,'I','comput_time','RIP_seq','COND_seq');
end

%% RIP and COND
% Recall the RIP condition
% Using d stands for delta for simplicity
%
% (1-d)|x|^2 < |Ax|^2 < (1+d) |x|^2
% We sample x and compute |Ax|/|x|
% Note that up to a scalar multiplication, 
% we really just care about the range of |Ax|/|x|
% Take the ratio R of the max(|Ax|/|x|) and min(|Ax|/|x|)
% we should have R = (1+d)/(1-d)
% In order to make d > 0.5, we need to have R < 3. 
% the line R==3 is drawn in the following figures

if plotRIP ==1
    figure;
    subplot(121);
    loglog(Mseq, RIP_seq, '-', 'LineWidth', 3)
    title('Change of RIP ratio in ALS')
    yline(0.5); text(10,0.5,'\delta=0.5')
    ylim([min(RIP_seq)*0.9, max(RIP_seq)+0.1])
    xline(num0); text(num0,0.6,'n_0')
    xline(num1); text(num1,0.6,'n_1')
    xlabel('Sample size M'); ylabel('RIP')
    subplot(122);
    loglog(Mseq, COND_seq, '-', 'LineWidth', 3)
    title('Change of Conditional number in ORALS')
    xline(num1); text(num1,median(COND_seq),'n_1')
    xline(num0); text(num0, median(COND_seq),'n_0')
    str_RIP      = sprintf('RIP_cond_N%i_n%i_M%i_L%i_obsnr',I.N,I.n,Mseq(end),I.steps);
    str_viscosity= ['_visc',num2str(I.viscosity)]; 
    str_RIP      = [str_RIP, num2str(I.obs_std),str_viscosity,'_ic',I.initial];
    xlabel('Sample size M'); ylabel('Condition number')
    figname = [I.SAVE_DIR_fig,str_RIP];
    set_positionFontsAll;
end