%%% test_M_vs_noise

%% Test Convergence in M for both ALS and ORALS
loadON = 1;     % Load previous results or not
saveON = 1;     % Save

% generate data of multple trajectories
num0 = (I.N-1+I.n); % /I.d; % /I.steps;
num1 = (I.N-1)*I.n; %/I.d; %/I.steps;

% Mseq consists of two parts.
% samll range and large range; merge together to compute
% Analyze separately
L = 8;   % length of Mseq in the large range;
%regu = 'RKHS'; % 'None'; %'lsqminnorm'; 'ID','RKHS'
regu0 = 'None';
regu1 = 'ID';
regu2 = 'RKHS';
regustr = {regu0, regu1, regu2};

Mseq_small = ceil(10.^linspace(log10(num0/2), log10(num1*4), L));
Mseq_large = ceil(10.^linspace(log10(num1*8), log10(2000), L));
Mseq = [Mseq_small, Mseq_large];

MM = Mseq(end)*3;  % Generate more data and sample from these data when testing convergence
fprintf('\nThe sequence of M is');disp(Mseq);

str_name      = sprintf('conv_in_M_N%i_n%i_M%i_L%i_obsnr',I.N,I.n,Mseq(end),I.steps);
str_viscosity = ['_visc',num2str(I.viscosity)];
str_name      = [str_name, num2str(I.obs_std),str_viscosity,'_ic',I.initial];




%% compute the estimator; the different regularization method

% generate data: and update I: exploration measure, basis matrix
fprintf('Generating trajectories, M = %i, and update I ...', MM);tic    % generate data if does not exist
[all_xpath,I] = get_M_path(I, MM, 'ParforProgressON', 1,'saveON', 1,'loadON', 1);
fprintf('done (%.2f sec).\n',toc);

test_num = 1; % multiple tests  --- seems not necessary?
Mseq_length = length(Mseq);

kernel_norm = sqrt(I.c_true'*I.dict_mat*I.c_true);
graph_norm  = norm(I.E,'fro');


all_c = cell(1,3);
all_E = cell(1,3);


for ii=1:3
    all_c{ii}.ALS_seq = cell(Mseq_length, test_num);
    all_E{ii}.ALS_seq = cell(Mseq_length, test_num);
    for b = 1:test_num                 % About parallel - parfor is used in ORALS & ALS, thus not here.
        path_id   = randperm(MM, Mseq(end)); % why do this? MM?
        test_path = all_xpath(path_id); % sample from a large collection of path
        for i = 1:Mseq_length
            fprintf('\n M-sequence:  %i out of %i : \n',i, Mseq_length);

            M = Mseq(i);
            [E_ALS, c_ALS,~,~,timeALS] ...
                = learn_kerne_graph_ALS(test_path(1:M), I, 'niter', 50, 'normalizeON', 1, 'plotON', 0,'reg_method', regustr{ii});

            all_c{ii}.ALS_seq{i, b} = c_ALS;

            all_E{ii}.ALS_seq{i, b} = E_ALS;
        end
    end
end

%% Compute error
error = cell(1,3);
for ii=1:3
    for i = 1:Mseq_length
        for b = 1:test_num
            E_ALS   = all_E{ii}.ALS_seq{i, b};
            c_ALS   = all_c{ii}.ALS_seq{i, b};
            error{ii}.k_als(i, b) = kernel_err(c_ALS, I);
            error{ii}.g_als(i, b) = graph_err(E_ALS, I);
            Z_ALS = get_Z_from_E_c(E_ALS, c_ALS);
            error{ii}.Z_als(i, b) = get_Z_error(Z_ALS, I);
        end
    end
end


%% All range of M plot
% for ii=1:3
[~, ind] = sort(Mseq);
ttl = 'Relative error as sample size increases';
draw_graph_kernel_M_ALS(Mseq, ind, I, error, regustr, ttl);
figname = [I.SAVE_DIR_fig,str_name,'regu',regustr];
% % set_positionFontsAll;
% end





%% Previous Code %%


% if ~exist(est_filename0,'file') || ~loadON

% est_filename0 = [I.SAVE_DIR,str_name,'regu',regu0,'.mat']; % saves estimator for each regu;
% est_filename1 = [I.SAVE_DIR,str_name,'regu',regu1,'.mat']; % saves estimator for each regu;
% est_filename2 = [I.SAVE_DIR,str_name,'regu',regu2,'.mat']; % saves estimator for each regu;



% compute the RIP and condition number if file does not exist
%     plotRIP = 1;
%     [RIP_seq,COND_seq]= compute_RIP_cond(rip_filename,all_xpath,I,Mseq,num0,num1,plotRIP);

%     all_c.ORALS_seq = cell(Mseq_length, test_num);
%     all_c.ORSVD_seq = cell(Mseq_length, test_num);


%     all_E.ORALS_seq = cell(Mseq_length, test_num);
%     all_E.ORSVD_seq = cell(Mseq_length, test_num);

%     comput_time = zeros(2,Mseq_length, test_num);

% compute the estimators

%             [E_ORSVD, c_ORSVD, E_ORALS, c_ORALS, Z_OR, ~, condA_orals,timeORALS] ...
%                 = learn_kernel_graph_ORALS(test_path(1:M), I, 'reg_method', regu0);

%             all_c.ORALS_seq{i, b} = c_ORALS;
%             all_c.ORSVD_seq{i, b} = c_ORSVD;


%             all_E.ORALS_seq{i, b} = E_ORALS;
%             all_E.ORSVD_seq{i, b} = E_ORSVD;

%             comput_time(1,i,b) = timeORALS;
%             comput_time(2,i,b) = timeALS;

%     save(est_filename0,'all_E','all_c','graph_norm','kernel_norm','I','comput_time');
%     save(data_filename);
% else
%     load(est_filename0,'all_E','all_c','graph_norm','kernel_norm','I','comput_time');
% end
%



%         E_ORSVD = all_E.ORSVD_seq{i, b};
%         E_ORALS = all_E.ORALS_seq{i, b};


%         c_ORALS = all_c.ORALS_seq{i, b};
%         c_ORSVD = all_c.ORSVD_seq{i, b};

%         error.k_orals(i, b) = kernel_err(c_ORALS, I);
%         error.k_orsvd(i, b) = kernel_err(c_ORSVD, I);

%         error.g_orals(i, b) = graph_err(E_ORALS, I);
%         error.g_orsvd(i, b) = graph_err(E_ORSVD, I);


%         Z_ORSVD = get_Z_from_E_c(E_ORSVD, c_ORSVD);
%         Z_ORALS = get_Z_from_E_c(E_ORALS, c_ORALS);

%         error.Z_orsvd(i, b) = get_Z_error(Z_ORSVD, I);
%         error.Z_orals(i, b) = get_Z_error(Z_ORALS, I);
%         error.Z_or(i, b)    = get_Z_error(Z_OR, I);

% plot_computTime(comput_time,Mseq,figname);


