% Convergence with respect to M and observation noise

close all;  clear all;   clc;    addpaths;   rng(10);

%% System settings
% Specify the settings here
% In order to fix the figures for paper

I = system_settings(); % setting of the IPS and graph and its integrator
I.viscosity = 0;    % viscosity (noise)
I.N = 7;               % number of agents
I.d = 2;               % dim of state vectors
I.t0       = 0;
I.dt       = 1e-3;     % time steps
I.steps    = 1;      % time steps
I.obs_std  = 0;     % observation noise
I.E = set_graph(I.N, 'sparsity', 0.4, 'plotON', 0);
I.initial = 'Unif_0_5';
I.basis_case = 8;

I = update_system_settings(I);

%% Test Convergence in stochastic forece for both ALS and ORALS
L_M = 25;
M_seq = floor(10.^linspace(1, 3.5, L_M));
L_o = 8;
obs_std_seq = 10.^(linspace(-2, -6, L_o));

MM = M_seq(end)*5;       % Generate more data and sample from them
test_num = 8;

% regu = 'RKHS'; % 'None'; %'lsqminnorm'; 'ID','RKHS_plain'
regu = 'RKHS';

loadON = 1;     % Load previous results or not
saveON = 1;     % Save

fprintf('\nThe sequence of M is');
for i = 1:L_M;fprintf('\t%d', M_seq(i));end;fprintf('\n');

str_name      = sprintf('conv_in_M_and_obs_std_N%i_basis_case%i_d%i_L%i_vis%i_', I.N, I.basis_case, I.d, I.steps, I.viscosity);
str_name      = [str_name,'_ic',I.initial, '_regu_', regu];

est_filename = [I.SAVE_DIR,str_name,'regu',regu,'.mat']; % saves estimator for each regu;
rip_filename = [I.SAVE_DIR,str_name,'_rip','.mat'];      %

fprintf(['File name is : ', str_name, '\n']);
%% compute the estimator; the different regularization method
if ~exist(est_filename,'file') || ~loadON
    
    all_c.ORALS_seq = cell(L_M, L_o, test_num);
    all_c.ORSVD_seq = cell(L_M, L_o, test_num);
    all_c.ALS_seq   = cell(L_M, L_o, test_num);
    
    all_E.ORALS_seq = cell(L_M, L_o, test_num);
    all_E.ORSVD_seq = cell(L_M, L_o, test_num);
    all_E.ALS_seq   = cell(L_M, L_o, test_num);
    
    
    % Generate data with largest M and 0 observation noise
    I.obs_std = 0;
    fprintf('Generating trajectories, M = %i, and update I ...', MM);tic
    [all_xpath,I] = get_M_path(I, MM, 'ParforProgressON', 1,'saveON', 1,'loadON', 0);
    fprintf('done (%.2f sec).\n',toc);
    
    for i = 1:L_M
        fprintf('\nM sequence:  %i out of %i : \n', i, L_M);
        M = M_seq(i);
        for o = 1:L_o
            obs_std = obs_std_seq(o);
            
            for b = 1:test_num
                path_id   = randperm(MM, M);
                test_path = all_xpath(path_id);         % sample from a large collection of path
                % Add observation noise to it
                for m = 1:M
                    test_path{m} = test_path{m} + obs_std*randn(size(test_path{m}));
                end
                
                
                [E_ORSVD, c_ORSVD, E_ORALS, c_ORALS, Z_OR, ~, ~, ~, ~] = learn_kernel_graph_ORALS(test_path, I, 'reg_method', regu);
                [E_ALS, c_ALS, ~, ~, ~]  = learn_kerne_graph_ALS(test_path, I, 'niter', 100, 'normalizeON', 1, 'plotON', 0,'reg_method', regu);
                
                all_c.ALS_seq{i, o, b} = c_ALS;
                all_c.ORALS_seq{i, o, b} = c_ORALS;
                all_c.ORSVD_seq{i, o, b} = c_ORSVD;
                
                all_E.ALS_seq{i, o, b} = E_ALS;
                all_E.ORALS_seq{i, o, b} = E_ORALS;
                all_E.ORSVD_seq{i, o, b} = E_ORSVD;
            end
        end
    end
    save(est_filename,'all_E','all_c', 'I');
else
    load(est_filename,'all_E','all_c', 'I');
end

%% Compute error
for i = 1:L_M
    for o = 1:L_o
        for b = 1:test_num
            E_ALS   = all_E.ALS_seq{i, o, b};
            E_ORSVD = all_E.ORSVD_seq{i, o, b};
            E_ORALS = all_E.ORALS_seq{i, o, b};
            
            c_ALS   = all_c.ALS_seq{i, o, b};
            c_ORALS = all_c.ORALS_seq{i, o, b};
            c_ORSVD = all_c.ORSVD_seq{i, o, b};
            
            error.k_orals(i, o, b) = kernel_err(c_ORALS, I);
            error.k_orsvd(i, o, b) = kernel_err(c_ORSVD, I);
            error.k_als(i, o, b) = kernel_err(c_ALS, I);
            
            error.g_orals(i, o, b) = graph_err(E_ORALS, I);
            error.g_orsvd(i, o, b) = graph_err(E_ORSVD, I);
            error.g_als(i, o, b) = graph_err(E_ALS, I);
        end
    end
end


%% plot M*
plot_M_star(M_seq, obs_std_seq, I, error.g_als, 'graph, ALS', 10)
plot_M_star(M_seq, obs_std_seq, I, error.k_als, 'kernrel, ALS', 10)
plot_M_star(M_seq, obs_std_seq, I, error.g_orals, 'graph, ORALS', 10)
plot_M_star(M_seq, obs_std_seq, I, error.k_orals, 'kernrel, ORALS', 10)
%%
function plot_M_star(M_seq, obs_std_seq, I, err, ttl, rate)
L_o = length(obs_std_seq);
set(gca,'colororder',parula(L_o))
rgb = get(gca,'colororder');
figure;


for o = 1:L_o
    obs_std = obs_std_seq(o);
    temp = squeeze(err(:, o, :));
    loglog(M_seq, temp, '.', 'color', rgb(o, :), 'MarkerSize', 15);hold on;
    YY = yline(rate*obs_std, 'color', rgb(o, :), 'Linewidth', 5, 'Label', ['10*\sigma']);
    val(o) = sum(mean(temp') > rate*obs_std);
    xx = xline(M_seq(val(o)), 'color', rgb(o, :), 'Linewidth', 5, 'Label', 'M_*', 'LabelOrientation', 'horizontal');
    xx.FontSize = 15;
end
ylabel('error')
title(['M_* such that error < ', num2str(rate), '*\sigma']);

xx = xline((I.N-1)*I.n, 'Label', 'N*n', 'LabelOrientation', 'horizontal');
xx.FontSize = 15;
xx = xline((I.N-1)+I.n, 'Label', 'N+n', 'LabelOrientation', 'horizontal');
xx.FontSize = 15;

xlim([((I.N-1)+I.n)/2, M_seq(end)*2])
Yticklabel = cell(L_o, 1);
for i = 1:L_o
    Yticklabel{i} = ['\sigma = ', sprintf('%.2e', obs_std_seq(i))];
end
cbh = colorbar; %Create Colorbar
cbh.Ticks = linspace(0, 1, L_o);
cbh.TickLabels = Yticklabel;


% title(ttl)
dim = [.7 .5 .2 .4];
annotation('textbox',dim,'String',ttl,'FitBoxToText','on');
end



%% plot 3D

%{
figure;
subplot(231)
plot_scatter_3(M_seq, obs_std_seq, error.g_als);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('graph error');title('ALS')

subplot(232);
plot_scatter_3(M_seq, obs_std_seq, error.g_orals);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('graph error');title('ORALS')

subplot(233)
plot_scatter_3(M_seq, obs_std_seq, error.g_orsvd);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('graph error');title('ORSVD')

subplot(234)
plot_scatter_3(M_seq, obs_std_seq, error.k_als);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('kernel error');title('ALS')

subplot(235)
plot_scatter_3(M_seq, obs_std_seq, error.k_orals);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('kernel error');title('ORALS')

subplot(236)
plot_scatter_3(M_seq, obs_std_seq, error.k_orsvd);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('kernel error');title('ORSVD')

%}

%% plot 2D

%{
figure;
subplot(231)
plot_scatter_2(M_seq, obs_std_seq, error.g_als);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('graph error');title('ALS')

subplot(232);
plot_scatter_2(M_seq, obs_std_seq, error.g_orals);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('graph error');title('ORALS')

subplot(233)
plot_scatter_2(M_seq, obs_std_seq, error.g_orsvd);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('graph error');title('ORSVD')

subplot(234)
plot_scatter_2(M_seq, obs_std_seq, error.k_als);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('kernel error');title('ALS')

subplot(235)
plot_scatter_2(M_seq, obs_std_seq, error.k_orals);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('kernel error');title('ORALS')

subplot(236)
plot_scatter_2(M_seq, obs_std_seq, error.k_orsvd);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('kernel error');title('ORSVD')

%}
%% plot 2D version 2

%{
figure;
subplot(231)
plot_scatter_2_v2(M_seq, obs_std_seq, error.g_als);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('graph error');title('ALS')

subplot(232);
plot_scatter_2_v2(M_seq, obs_std_seq, error.g_orals);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('graph error');title('ORALS')

subplot(233)
plot_scatter_2_v2(M_seq, obs_std_seq, error.g_orsvd);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('graph error');title('ORSVD')

subplot(234)
plot_scatter_2_v2(M_seq, obs_std_seq, error.k_als);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('kernel error');title('ALS')

subplot(235)
plot_scatter_2_v2(M_seq, obs_std_seq, error.k_orals);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('kernel error');title('ORALS')

subplot(236)
plot_scatter_2_v2(M_seq, obs_std_seq, error.k_orsvd);
xlabel('log_{10}(M)');ylabel('log_{10}(obs std)');zlabel('kernel error');title('ORSVD')

%}


%% plot
% err_top = max([error.k_als;error.k_orals;error.k_orsvd;error.g_als;error.g_orals;error.g_orsvd], [], 'all');
% err_bot = min([error.k_als;error.k_orals;error.k_orsvd;error.g_als;error.g_orals;error.g_orsvd], [], 'all');
%
% y_bot = min(err_bot, I.obs_std/2);
% y_top = err_top*2;
%
% figure;
% subplot(231);loglog(M_seq, error.k_als  , 'r.');ylim([y_bot, y_top]);grid on;xlabel('M');ylabel('kernel error');xline(I.N-1+I.n);yline(I.obs_std);xline((I.N-1)*I.n);title('ALS');
% subplot(232);loglog(M_seq, error.k_orals, 'r.');ylim([y_bot, y_top]);grid on;xlabel('M');ylabel('kernel error');xline(I.N-1+I.n);yline(I.obs_std);xline((I.N-1)*I.n);title('ORALS');
% subplot(233);loglog(M_seq, error.k_orsvd, 'r.');ylim([y_bot, y_top]);grid on;xlabel('M');ylabel('kernel error');xline(I.N-1+I.n);yline(I.obs_std);xline((I.N-1)*I.n);title('ORSVD');
%
% subplot(234);loglog(M_seq, error.g_als  , 'r.');ylim([y_bot, y_top]);grid on;xlabel('M');ylabel('graph error');xline(I.N-1+I.n);yline(I.obs_std);xline((I.N-1)*I.n);
% subplot(235);loglog(M_seq, error.g_orals, 'r.');ylim([y_bot, y_top]);grid on;xlabel('M');ylabel('graph error');xline(I.N-1+I.n);yline(I.obs_std);xline((I.N-1)*I.n);
% subplot(236);loglog(M_seq, error.g_orsvd, 'r.');ylim([y_bot, y_top]);grid on;xlabel('M');ylabel('graph error');xline(I.N-1+I.n);yline(I.obs_std);xline((I.N-1)*I.n);
%
% set(gcf,'Position',[100 100 1000 500])
% sgtitle(['Convergence with M, obs std = ', num2str(I.obs_std), ', regu = ', regu], 'Interpreter','none')





%%
function plot_scatter_2_v2(M_seq, obs_std_seq, err)
[~, L_o, ~] = size(err);


set(gca,'colororder',parula(L_o))
rgb = get(gca,'colororder');

for o = 1:L_o
    temp = squeeze(err(:, o, :));
    val = max(temp, [], 2);
    loglog(M_seq, val, 'Color', rgb(o, :), 'LineWidth', 3);
    hold on;
    val = min(temp, [], 2);
    loglog(M_seq, val, 'Color', rgb(o, :), 'LineWidth', 3);
    hold on;
    grid on;
end


Yticklabel = cell(L_o, 1);
for i = 1:L_o
    Yticklabel{i} = ['obs_std = ', num2str(obs_std_seq(i))];
end
cbh = colorbar; %Create Colorbar
cbh.Ticks = linspace(0, 1, L_o);
cbh.TickLabels = Yticklabel;

end


%%
function plot_scatter_2(M_seq, obs_std_seq, err)
[L_M, L_o, ~] = size(err);


set(gca,'colororder',parula(L_o))
rgb = get(gca,'colororder');

for m = 1:L_M
    for o = 1:L_o
        M = M_seq(m);
        obs_std = obs_std_seq(o);
        temp = squeeze(err(m, o, :));
        plot3(log10(M), log10(obs_std), log10(max(temp)), '.', 'MarkerSize', 15, 'MarkerFaceColor', rgb(o, :), 'MarkerEdgeColor', rgb(o, :));
        plot3(log10(M), log10(obs_std), log10(min(temp)), '.', 'MarkerSize', 15, 'MarkerFaceColor', rgb(o, :), 'MarkerEdgeColor', rgb(o, :));
        hold on;
        grid on;
    end
end


Yticklabel = cell(L_o, 1);
for i = 1:L_o
    Yticklabel{i} = ['obs_std = ', num2str(obs_std_seq(i))];
end
cbh = colorbar; %Create Colorbar
cbh.Ticks = linspace(0, 1, L_o);
cbh.TickLabels = Yticklabel;

end

%%
function plot_scatter_3(M_seq, obs_std_seq, err)
[L_M, L_o, test_num] = size(err);


set(gca,'colororder',parula(L_o))
rgb = get(gca,'colororder');

for m = 1:L_M
    for o = 1:L_o
        M = M_seq(m);
        obs_std = obs_std_seq(o);
        for b = 1:test_num
            %             plot3(log10(M), log10(obs_std), log10(err(m, o, b)), 'MarkerFaceColor', rgb(o, :), 'MarkerEdgeColor', rgb(o, :));
            plot3(log10(M), log10(obs_std), log10(err(m, o, b)), '.', 'MarkerSize', 15, 'MarkerFaceColor', rgb(o, :), 'MarkerEdgeColor', rgb(o, :));
            hold on;
            grid on;
        end
    end
end


Yticklabel = cell(L_o, 1);
for i = 1:L_o
    Yticklabel{i} = ['obs_std = ', num2str(obs_std_seq(i))];
end
cbh = colorbar; %Create Colorbar
cbh.Ticks = linspace(0, 1, L_o);
cbh.TickLabels = Yticklabel;

end


%% 
