%% generate a figure
close all;  clear all;   clc;    addpaths;   rng(10);

%% System settings
% Specify the settings here
% In order to fix the figures for paper
I = system_settings(); % setting of the IPS and graph and its integrator
I.viscosity = 0;    % viscosity (noise)
I.N = 12;               % number of agents
I.d = 2;               % dim of state vectors
I.t0       = 0;
I.dt       = 1e-3;     % time steps
I.steps    = 500000;      % time steps
I.obs_std  = 1e-7;     % observation noise
I.sparsity = 1; 
I.E = set_graph(I.N, 'sparsity', I.sparsity, 'plotON', 0);
I.initial = 'Unif_0_5';
I.basis_case = 1;     % repulsion-attraction---- with N=12, sparsity=1, forming rotating clusters; 
I = update_system_settings(I);


% %% load settings
% I.viscosity = 0;       % viscosity (noise)
% I.N = 20;              % number of agents
% I.d = 2;               % dim of state vectors
% I.t0       = 0;
% I.dt       = 1e-3;
% I.steps    = 5000;
% I.obs_std  = 1e-3;     % observation noise
%
% %% Graph
% sparsity = 0.4;
% I.E = set_graph(I.N, 'sparsity', sparsity, 'plotON', 0);
%
% plot_graph(I.E);
%
% % This is to construct the 'circle' graph
% N = I.N;
% I.E = [zeros(1, N-1), 1;eye(N-1, N-1), zeros(N-1, 1)];
% I.E = ones(N, N) - eye(N);
% %% Initial
% % I.initial = 'DN_1_1_0.5';  % or 'Unif' *********
% %  I.initial = 'N_0_1';  % or 'Unif'
% I.initial = 'Unif_0_0.1';


%% % %  dictionary for the kernel
% basis_case = 6;
% switch basis_case
%     case 1         % % % % dictionary1
%         p = 8;q = 2;cut = 0.5;
%         I.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
%         I.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
%         I.dict{3} = @(x) abs(x)<=cut;
%         I.dict{4} = @(x) sin(x);
%         I.dict{5} = @(x) x.^2;
%         I.c_true = [-0.3333, 1.3333, -160, 10, -0.5]';
%     case 2           % % % % dictionary2
%         I.dict{1} = @(x) x.^2;
%         I.dict{2} = @(x) abs(x);
%         I.dict{3} = @(x) cos(x);
%         I.c_true = [.1, 0.1, 3]';
%     case 3           % % % dictionary3
%         n = 3;
%         dict = cell(n, 1);
%         for i = 1:n
%             dict{i} = @(x) sin(x*i + i);
%         end
%         I.dict = dict;
%         I.c_true =  ones(n, 1);% [.1, 0.1, 3]';
%     case 4         % % % % dictionary1
%         p = 2;q = .2;cut = 0.5;
%         I.dict{1} = @(x)  x.^(-p-1).*(abs(x)>=cut);
%         I.dict{2} = @(x) - x.^(-q-1).*(abs(x)>=cut);
%         I.dict{3} = @(x) - 1.*(abs(x)<=cut);
%         c1 = -1; c2= 1;
%         c3 = (c1*I.dict{1}(cut) + c2*I.dict{2}(cut))/I.dict{3}(cut);
%         I.c_true = [c1, c2,c3]';
%     case 5
%         I.dict{1} = @(x) x.^2;
%         I.dict{2} = @(x) x.^4;
%         I.dict{3} = @(x) 0*x+1;
%         I.c_true = [-1, -1, -0.01]';
%
%     case 6
%         p = 8;q = 2;cut = 0.5;
%         I.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
%         I.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
%         I.dict{3} = @(x) abs(x)<=cut;
%         %         I.dict{4} = @(x) sin(x);
%         %         I.dict{5} = @(x) x.^2;
%         I.c_true = [-0.3333, 1.3333, -160]';
% end
% I = update_system_settings(I);

 % figure;fplot(I.phi_kernel, [0, 5],'linewidth',2); title('Interaction kernel')

%% Demo (1 traj) generate particles
progressON = 0;
xpath = graph_forward_model(I, I.X0, progressON,1);

%% Plot trajectories and the graph
plotON = 1;
if plotON==1 &&  I.d == 2
    graph_plot_motion(xpath, I, plotON);
end
% TODO: plot a few time instants in 2D, with color indicating time.

%% trajectory in 3D: with time a the 3rd dimension
% when the particles "merge", need to zoom out by changing tseq(1) and tN;
figure;
subplot(2,2,[1,3])
tN = I.steps *0.8;
tseq = 1:ceil(tN /141):tN;
for i = 1:I.N
    xt  = squeeze(xpath(i,1,tseq));
    yt  = squeeze(xpath(i,2,tseq));
    plot3(xt,yt,tseq*I.dt,'linewidth',4);  hold on;
    xlabel('x(t)'); ylabel('y(t)'); zlabel('Time t')
end
title('Trajectories of particles')

subplot(2,2,2);
G = digraph(I.E,'omitselfloops');
plot(G,'Layout','force', 'EdgeLabel', round(G.Edges.Weight,1),'linewidth',1,'markersize',12); box off 
title('Graph of network')
subplot(2,2,4)
fplot(I.phi_kernel, [0, 5],'linewidth',2); title('Interaction kernel')

SAVE_DIR = [getenv('HOME'),'/IPS_Graph_data/'];
if ~exist(SAVE_DIR,'dir'); mkdir(SAVE_DIR); end
str_name      = sprintf('trajGraph_N%i_n%i_L%i_sparcity%1.1f',I.N,I.n,I.steps,I.sparsity);
str_viscosity = ['_visc',num2str(I.viscosity)];
str_name      = [str_name, 'obs',num2str(I.obs_std),str_viscosity,'_ic',I.initial];
figname = [SAVE_DIR,str_name];
savefig([figname,'.fig']);
set_positionFontsAll;

I.SAVE_DIR = [getenv('HOME'),'/IPS_Graph_data/'];
if ~exist(I.SAVE_DIR,'dir'); mkdir(I.SAVE_DIR); end

I.SAVE_DIR_fig = [getenv('HOME'),'/IPS_Graph_data/figures/'];
if ~exist(I.SAVE_DIR_fig,'dir'); mkdir(I.SAVE_DIR_fig); end

