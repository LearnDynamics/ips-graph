function [all_path, I] = get_M_path_social_range(I, M, varargin)
%% Input parser
p = inputParser;
addRequired(p, 'I');
addRequired(p, 'M');
addOptional(p, 'saveON', 0);
addOptional(p, 'ParforProgressON', 0);
addOptional(p, 'loadON', 0);


parse(p, I, M, varargin{:});
I = p.Results.I;
M = p.Results.M;
saveON = p.Results.saveON;
loadON = p.Results.loadON;
ParforProgressON = p.Results.ParforProgressON;
%% Check if we have local path file
if loadON
    DataPath = get_Data_Path(M, I);
    if exist(DataPath, 'file')
        fprintf('Loading local trajectory file...')
        load(DataPath, 'all_path','I');
        return
    end
    
end
%%
progressON = 0;
rk4ON = 0;
all_path = cell(M, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Why don't we use parfor for all cases? %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if M> 1e4  % use parfor is M is large
%     parfor i = 1:M
%         X0 = set_particle_initial_all_dim(N, d, initial);
%         xpath        = graph_forward_model(I, X0, progressON, rk4ON);      % using RK4, the result will not converge.
%         all_path{i} = xpath + I.obs_std* randn(size(xpath));
%     end
% else
%     for i = 1:M
%         X0 = set_particle_initial_all_dim(N, d, initial);
%         xpath        = graph_forward_model(I, X0, progressON, rk4ON);      % using RK4, the result will not converge.
%         all_path{i} = xpath + I.obs_std* randn(size(xpath));
%     end
% end

N = I.N;
d = I.d;
initial = I.initial;
obs_std = I.obs_std;


if ParforProgressON;WaitMessage = parfor_wait(M, 'Waitbar', true);
else; WaitMessage = 0;end

parfor i = 1:M
    X0 = set_particle_initial_all_dim(N, d, initial);
    xpath        = graph_forward_model_social_range(I, X0, progressON, rk4ON);      % using RK4, the result will not converge.
    all_path{i} = xpath + obs_std* randn(size(xpath));
    
    if ParforProgressON;WaitMessage.Send;end
    
end

if ParforProgressON;WaitMessage.Destroy;end

fprintf('Data generated; updating I ...');
%% update I: the exploration measure rho and the basis matrix using all paths
% I = update_dict_mat(I, all_path); % compute the exploration measure rho and the basis matrix: using all paths
% I.kernel_norm = sqrt(I.c_true'*I.dict_mat*I.c_true);
%% save data
if saveON
    %     setting_name = ['M_', num2str(M), '_N_', num2str(N), '_steps_', num2str(I.steps)];
    %     DataPath = [I.SAVE_DIR, 'all_path_', setting_name, '.mat'];
    DataPath = get_Data_Path(M, I);
    save(DataPath, 'all_path','I');
end

end


function DataPath = get_Data_Path(M, I)
% setting_name = ['M_', num2str(M), '_N_', num2str(I.N), '_steps_', num2str(I.steps),'_ic',I.initial];

str_name      = sprintf('social_range_M_%i_N_%i_steps_%i_basis_case%i_d%i_obsnr', M, I.N, I.steps, I.basis_case, I.d);
str_name      = [str_name, num2str(I.obs_std),'_ic',I.initial, '_vis_', num2str(I.viscosity)];




DataPath = [I.SAVE_DIR, 'all_path_', str_name, '.mat'];
end




