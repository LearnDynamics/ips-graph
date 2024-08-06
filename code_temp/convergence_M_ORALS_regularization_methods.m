% Learn the adjacency matrix and kernel parameter of an Interacting Particle System on graph.
close all;  clear all;   clc;    addpaths;

%% load settings
I = system_settings();  % setting of the IPS and graph and its integrator

%% Test Convergence in M for ORALS
% sequence of M
num1 = (I.N-1)*I.n;
L = 8; % length of Mseq;
Mseq = ceil(10.^linspace(log10(num1/2), log10(1e5), L));
MM = Mseq(end)*3;  % Generate more data and sample from these data when testing convergence
fprintf('\nThe sequence of M is');disp(Mseq);

% sequence of random tests
test_num = 5; % multiple tests  --- seems not necessary

% sequence of regu methods
all_regu = {'lsqminnorm', 'ID', 'RKHS'};
R = length(all_regu);
%% generate data: and update I: exploration measure, basis matrix
fprintf('Generating trajectories, M = %i, and update I ...', MM);tic    % generate data if does not exist
[all_xpath, I] = get_M_path(I, MM, 'ParforProgressON', 1,'saveON', 1,'loadON', 1);
fprintf('done (%.2f sec).\n',toc);

%% prepare for the loop
all_c.ORALS_seq = cell(L, test_num, R);
all_c.ORSVD_seq = cell(L, test_num, R);
all_E.ORALS_seq = cell(L, test_num, R);
all_E.ORSVD_seq = cell(L, test_num, R);
comput_time = zeros(L, test_num, R);
all_COND.A = cell(L, test_num, R);
all_COND.L = cell(L, test_num, R);

for i = 1:L
    for b = 1:test_num
        path_id   = randperm(MM, Mseq(end));
        test_path = all_xpath(path_id); % sample from a large collection of path
        
        fprintf('\n M-sequence:  %i out of %i : \n',i, L);
        M = Mseq(i);
        for r = 1:R
            regu = all_regu{r};
            [E_ORSVD, c_ORSVD, E_ORALS, c_ORALS, Z_OR, ~, condA_orals, condL, timeORALS] = ...
                learn_kernel_graph_ORALS(test_path(1:M), I, 'reg_method', regu);
            
            all_c.ORALS_seq{i, b, r} = c_ORALS;
            all_c.ORSVD_seq{i, b, r} = c_ORSVD;
            
            all_E.ORALS_seq{i, b, r} = E_ORALS;
            all_E.ORSVD_seq{i, b, r} = E_ORSVD;
            
            comput_time(i, b, r) = timeORALS;
            all_COND.A{i, b, r} = condA_orals;
            all_COND.L{i, b, r} = condL;
        end
    end
end



%% Compute error
error.k_orals = zeros(L, test_num, R);
error.g_orals = zeros(L, test_num, R);
error.k_orsvd = zeros(L, test_num, R);
error.g_orsvd = zeros(L, test_num, R);
for i = 1:L
    for b = 1:test_num
        for r = 1:R
            E_ORSVD = all_E.ORSVD_seq{i, b, r};
            E_ORALS = all_E.ORALS_seq{i, b, r};
            
            c_ORALS = all_c.ORALS_seq{i, b, r};
            c_ORSVD = all_c.ORSVD_seq{i, b, r};
            
            error.k_orals(i, b, r) = kernel_err(c_ORALS, I);
            error.k_orsvd(i, b, r) = kernel_err(c_ORSVD, I);
            
            error.g_orals(i, b, r) = graph_err(E_ORALS, I);
            error.g_orsvd(i, b, r) = graph_err(E_ORSVD, I);
            
            Z_ORSVD = get_Z_from_E_c(E_ORSVD, c_ORSVD);
            Z_ORALS = get_Z_from_E_c(E_ORALS, c_ORALS);
            
            error.Z_orsvd(i, b, r) = get_Z_error(Z_ORSVD, I);
            error.Z_orals(i, b, r) = get_Z_error(Z_ORALS, I);
            error.Z_or(i, b, r)    = get_Z_error(Z_OR, I);
        end
    end
end

%%
figure;
for r = 1:R
    subplot(1, R, r);
    loglog(Mseq, squeeze(error.k_orals(:, :, r)), 'r.', 'MarkerSize', 10);
    grid on
    title(all_regu{r})
    ylim([1e-5, 1e2])
end
sgtitle('Kernel')

figure;
for r = 1:R
    subplot(1, R, r);
    loglog(Mseq, squeeze(error.g_orals(:, :, r)), 'r.', 'MarkerSize', 10);
    grid on
    title(all_regu{r})
    ylim([1e-5, 1e2])
end
sgtitle('Graph')

figure;
for r = 1:R
    subplot(1, R, r);
    loglog(Mseq, squeeze(error.Z_orals(:, :, r)), 'r.', 'MarkerSize', 10);
    grid on
    title(all_regu{r})
    ylim([1e-5, 1e2])
end
sgtitle('Z')

%%
newcolors = {'red','black','blue'};
figure;
hold on;
for b = 1:test_num
    for m = 1:L
        for i = 1:3
            all_COND.A{m, b}(i)
            plot(log10(Mseq(m)), log10(all_COND.A{m, b}(i)), '.', 'MarkerSize', 10, 'Color', newcolors{i});
        end
    end
end

%%
figure;
hold on;
for b = 1:test_num
    for m = 1:L
        for i = 1:3
            plot(log10(Mseq(m)), log10(all_COND.L{m, b}(i)), '.', 'MarkerSize', 10, 'Color', newcolors{i});
        end
    end
end




