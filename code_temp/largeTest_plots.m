%% Setup
loadDirs = {'/Volumes/LargeSSD/DataAnalyses/LIKSonGraphs_davincileonardo31','/Volumes/LargeSSD/DataAnalyses/LIKSonGraphs_davincileonardo32'};

%% Merge results in multiple directories if needed.
% Prepare the variables. The output of largeTest is in a cell array called run
for q = 1:length(loadDirs)
    load(sprintf('%s/LIKSonGraphs_largeTest_setup.mat',loadDirs{q}));
    for x_idx_1 = length(paramStructs):-1:1
        N{q}(x_idx_1) = paramStructs{x_idx_1}.N;
        d{q}(x_idx_1) = paramStructs{x_idx_1}.d;
        viscosity{q}(x_idx_1) = paramStructs{x_idx_1}.viscosity;
        A_sparsity{q}(x_idx_1) = paramStructs{x_idx_1}.A_sparsity;
        L{q}(x_idx_1) = paramStructs{x_idx_1}.L;
        obs_std{q}(x_idx_1) = paramStructs{x_idx_1}.obs_std;
        M{q}(x_idx_1) = paramStructs{x_idx_1}.M;
        n{q}(x_idx_1) = paramStructs{x_idx_1}.n;
    end

    uniqM{q} = unique(M{q});
    uniqL{q} = unique(L{q});
    uniqN{q} = unique(N{q});
    uniqn{q} = unique(n{q});
    uniqd{q} = unique(d{q});
    uniqviscosity{q} = unique(viscosity{q});
    uniqA_sparsity{q} = unique(A_sparsity{q});
    uniqobs_std{q} = unique(obs_std{q});
    ntrials{q} = expSetup.ntrials;

    fprintf('\nThese are the ranges of parameters in the experiment %d:',q)
    fprintf('\nM    :[');fprintf(' %d',uniqM{q});fprintf(' ]');
    fprintf('\nL    :[');fprintf(' %d',uniqL{q});fprintf(' ]');
    fprintf('\nN    :[');fprintf(' %d',uniqN{q});fprintf(' ]');
    fprintf('\nn    :[');fprintf(' %d',uniqn{q});fprintf(' ]');
    fprintf('\nd    :[');fprintf(' %d',uniqd{q});fprintf(' ]');
    fprintf('\nnu   :[');fprintf(' %e',uniqviscosity{q});fprintf(' ]');
    fprintf('\nAsp  :[');fprintf(' %e',uniqA_sparsity{q});fprintf(' ]');
    fprintf('\nsigma:[');fprintf(' %e',uniqobs_std{q});fprintf(' ]');
    fprintf('\ntrials :%d',ntrials{q});
end

%% Set up quantities for 1-D plots
ToPlot.xVarName1    = 'M';

ToPlot.yVarIdxs     = [1,2,4];
ToPlot.yVarName     = {'kernel_err','A_err','Z_err','pathTestErr'};
ToPlot.yVarLegend   = {'Int. Kernel error','Graph error','Z error','Traj. error'};
ToPlot.yVarColor    = {hex2rgb('f47a00'),hex2rgb('007191'),[0,1,0],hex2rgb('d31f11')};
ToPlot.yVarMarker   = {'x','o','d','s'};

ToPlot.L_values             = uniqL{1}(1);
ToPlot.n_values             = uniqn{1}(3);
ToPlot.N_values             = uniqN{1}(3);
ToPlot.d_values             = uniqd{1}(1);
ToPlot.A_sparsity_values    = uniqA_sparsity{1}(1);
ToPlot.viscosity_values     = uniqviscosity{1}(1);
ToPlot.obs_std_values       = uniqobs_std{1}(1);

% Load and create table for 1-D plots
ys_ALS = []; ys_ORALS = [];

% organize the variable of interest, and take union over multiple rounds of experiments
uniqxVarunion = [];
for q = 1:length(loadDirs)
    eval(sprintf('xVar{q}      = %s{q};',ToPlot.xVarName1));
    eval(sprintf('uniqxVar{q}  = uniq%s{q};',ToPlot.xVarName1));
    uniqxVarunion = union(uniqxVarunion,uniqxVar{q})';
end

for x_idx_1 = 1:length(uniqxVarunion)                                                                                           % Loop through all the unique values of the variable of interest
    for q = 1:length(loadDirs)
        idxs_1  = find( xVar{q}==uniqxVarunion(x_idx_1) );
        k_ok    = 1;
        for k = 1:length(idxs_1)                                                                                                % loop through all the experiments with the given value of the x variable
            if  ismember(L{q}(idxs_1(k)),ToPlot.L_values) && ...
                    ismember(n{q}(idxs_1(k)),ToPlot.n_values) && ...
                    ismember(N{q}(idxs_1(k)),ToPlot.N_values) && ...
                    ismember(d{q}(idxs_1(k)),ToPlot.d_values) && ...
                    ismember(A_sparsity{q}(idxs_1(k)),ToPlot.A_sparsity_values) && ...
                    ismember(viscosity{q}(idxs_1(k)),ToPlot.viscosity_values) && ...
                    ismember(obs_std{q}(idxs_1(k)),ToPlot.obs_std_values)

                try
                    load(sprintf('%s/LIKSonGraphs_largeTest_%d.mat',loadDirs{q},idxs_1(k)),'cur_run');
                catch
                    continue
                end
                for pp = 1:length(cur_run.estALS)                                                                               % loop through the multiple runs of the same experiment
                    %                if isempty(ys), ys = zeros( length(uniqxVar), length(ToPlot.depErrVar),length(idxs),length(cur_run.estALS) ); end
                    for errVars_idx = 1:length(ToPlot.yVarIdxs)                                                                 % loop through the y-variables to be reported
                        if ~isfield(cur_run.estORALS{pp},'kernel_err')
                            cur_run.estORALS{pp}.kernel_err = cur_run.estORALS{pp}.c_err;
                        end
                        ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp) = cur_run.estALS{pp}.(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)});
                        ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp) = cur_run.estORALS{pp}.(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)});
                        if strcmp(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)},'pathTestErr')                                  % normalize pathTestErr by dividing by L. MM:TBD: length of test paths should not need to be L
                            ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp) = ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp)/L{q}(idxs_1(k));
                            ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp) = ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp)/L{q}(idxs_1(k));
                        end
                        if strcmp(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)},'A_err')                                        % normalize A_err
                            ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp) = ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp)/(N{q}(idxs_1(k))^2);
                            ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp) = ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp)/(N{q}(idxs_1(k))^2);
                        end
                    end
                    %fprintf('.')
                end
                expIdx_ys(x_idx_1,k_ok)=idxs_1(k);
                k_ok = k_ok + 1;
            end
        end
    end
end

ys_ALS = reshape(ys_ALS,[size(ys_ALS,1),size(ys_ALS,2),size(ys_ALS,3),size(ys_ALS,4)*size(ys_ALS,5)]);
ys_ORALS = reshape(ys_ORALS,[size(ys_ORALS,1),size(ys_ORALS,2),size(ys_ORALS,3),size(ys_ORALS,4)*size(ys_ORALS,5)]);

%% Plot
plotHandle = [];
plotFaceAlpha = 0.02;
figure;
for errVars_idx = 1:length(ToPlot.yVarIdxs)                                                                                     % loop through the y-variables to be reported
    loc_ys_mean             = squeeze(mean(ys_ALS(:,errVars_idx,:,:),4));                                                       % ALS results
    loc_ys_std              = squeeze(std(ys_ALS(:,errVars_idx,:,:),[],4));
    plotHandle(end+1)       = plot(log10(uniqxVarunion),log10(loc_ys_mean),'LineWidth',2,'LineStyle','-','Marker',ToPlot.yVarMarker{ToPlot.yVarIdxs(errVars_idx)},'MarkerFaceColor',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)},'Color',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)});
    hold on;
    patch(log10([uniqxVarunion;flipud(uniqxVarunion)]),log10([loc_ys_mean+loc_ys_std;flipud(loc_ys_mean-loc_ys_std)]),get(plotHandle(errVars_idx),'color'), 'FaceAlpha',plotFaceAlpha,'EdgeColor',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)});hold on;
    plot(log10(uniqxVarunion),log10(loc_ys_mean+loc_ys_std),'-','LineWidth',0.25,'color',get(plotHandle(errVars_idx),'color'));
    plot(log10(uniqxVarunion),log10(max(loc_ys_mean-loc_ys_std,0)),'-','LineWidth',0.25,'color',get(plotHandle(errVars_idx),'color'));
end
for errVars_idx = 1:length(ToPlot.yVarIdxs)                                                                                     % loop through the y-variables to be reported
    loc_ys_mean             = squeeze(mean(ys_ORALS(:,errVars_idx,:,:),4));                                                     % ORALS results
    loc_ys_std              = squeeze(std(ys_ORALS(:,errVars_idx,:,:),[],4));
    plotHandle(end+1)       = plot(log10(uniqxVarunion),log10(loc_ys_mean),'LineWidth',2,'color',get(plotHandle(errVars_idx),'color'),'LineStyle','-.','Marker',ToPlot.yVarMarker{ToPlot.yVarIdxs(errVars_idx)},'MarkerFaceColor',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)},'Color',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)});
    patch(log10([uniqxVarunion;flipud(uniqxVarunion)]),log10([loc_ys_mean+loc_ys_std;flipud(loc_ys_mean-loc_ys_std)]),get(plotHandle(errVars_idx),'color'), 'FaceAlpha',plotFaceAlpha);hold on;
    plot(log10(uniqxVarunion),log10(loc_ys_mean+loc_ys_std),'-','LineWidth',0.25,'color',get(plotHandle(errVars_idx),'color'));
    plot(log10(uniqxVarunion),log10(max(loc_ys_mean-loc_ys_std,0)),'-','LineWidth',0.25,'color',get(plotHandle(errVars_idx),'color'));
end
axis tight;grid on
xlabel('log_{10}M','FontSize',14); ylabel('Estimation Error','FontSize',14);
legend_labels = {};
for k = 1:length(ToPlot.yVarIdxs), legend_labels{end+1} = [ToPlot.yVarLegend{ToPlot.yVarIdxs(k)}, ' (ALS)']; end
for k = 1:length(ToPlot.yVarIdxs), legend_labels{end+1} = [ToPlot.yVarLegend{ToPlot.yVarIdxs(k)}, ' (ORALS)']; end

curaxes = gca;
vertLine_x= (ToPlot.N_values^2+ToPlot.n_values)/(ToPlot.L_values*ToPlot.d_values*ToPlot.N_values);
line(log10([vertLine_x,vertLine_x]),(([curaxes.YLim(1),curaxes.YLim(2)])),'LineWidth',1,'color',[0.1,0.1,0.1],'LineStyle','--');
text(log10(vertLine_x),curaxes.YLim(1)+0.25,'$\frac{N^2+p}{LdN}$', "Interpreter", 'latex','FontSize',14);
%set(curaxes,'XTick',union(get(curaxes,'XTick'),log10(vertLine_x)));
%xTickLabels=get(curaxes,'XTickLabel'); xTickLabels{2}='$\frac{N^2+p}{LdN}$';set(gca,'XTickLabel',xTickLabels);set(curaxes,'TickLabelInterpreter','latex','FontSize',14);

vertLine_x = (ToPlot.N_values^2*ToPlot.n_values)/(ToPlot.L_values*ToPlot.d_values*ToPlot.N_values);
line(log10([vertLine_x,vertLine_x]),([curaxes.YLim(1),curaxes.YLim(2)]),'LineWidth',1,'color',[0.1,0.1,0.1],'LineStyle','--');
text(log10(vertLine_x),curaxes.YLim(1)+0.25,'$\frac{N^2p}{LdN}$', "Interpreter", 'latex','FontSize',14);
%set(curaxes,'XTick',union(get(gca,'XTick'),log10(vertLine_x)));
%xTickLabels=get(curaxes,'XTickLabel'); xTickLabels{5}='$\frac{N^2p}{LdN}$';set(gca,'XTickLabel',xTickLabels);set(curaxes,'TickLabelInterpreter','latex','FontSize',14);
%legend(plotHandle,legend_labels,'FontSize',12,'Location','NorthEast');
%gcfpos = get(gcf,'Position'); set(gcf,'Position',[gcfpos(1),gcfpos(2),gcfpos(3)*1.2,gcfpos(4)*1.2]);

%% Set up quantities for 1-D plots
ToPlot.L_values             = uniqL{1}(3);
ToPlot.n_values             = uniqn{1}(3);
ToPlot.N_values             = uniqN{1}(3);
ToPlot.d_values             = uniqd{1}(1);
ToPlot.A_sparsity_values    = uniqA_sparsity{1}(1);
ToPlot.viscosity_values     = uniqviscosity{1}(1);
ToPlot.obs_std_values       = uniqobs_std{1}(1);

% Load and create table for 1-D plots
ys_ALS = []; ys_ORALS = [];

% organize the variable of interest, and take union over multiple rounds of experiments
uniqxVarunion = [];
for q = 1:length(loadDirs)
    eval(sprintf('xVar{q}      = %s{q};',ToPlot.xVarName1));
    eval(sprintf('uniqxVar{q}  = uniq%s{q};',ToPlot.xVarName1));
    uniqxVarunion = union(uniqxVarunion,uniqxVar{q})';
end

for x_idx_1 = 1:length(uniqxVarunion)                                                                                           % Loop through all the unique values of the variable of interest
    for q = 1:length(loadDirs)
        idxs_1  = find( xVar{q}==uniqxVarunion(x_idx_1) );
        k_ok    = 1;
        for k = 1:length(idxs_1)                                                                                                % loop through all the experiments with the given value of the x variable
            if  ismember(L{q}(idxs_1(k)),ToPlot.L_values) && ...
                    ismember(n{q}(idxs_1(k)),ToPlot.n_values) && ...
                    ismember(N{q}(idxs_1(k)),ToPlot.N_values) && ...
                    ismember(d{q}(idxs_1(k)),ToPlot.d_values) && ...
                    ismember(A_sparsity{q}(idxs_1(k)),ToPlot.A_sparsity_values) && ...
                    ismember(viscosity{q}(idxs_1(k)),ToPlot.viscosity_values) && ...
                    ismember(obs_std{q}(idxs_1(k)),ToPlot.obs_std_values)

                try
                    load(sprintf('%s/LIKSonGraphs_largeTest_%d.mat',loadDirs{q},idxs_1(k)),'cur_run');
                catch
                    continue
                end
                for pp = 1:length(cur_run.estALS)                                                                               % loop through the multiple runs of the same experiment
                    %                if isempty(ys), ys = zeros( length(uniqxVar), length(ToPlot.depErrVar),length(idxs),length(cur_run.estALS) ); end
                    for errVars_idx = 1:length(ToPlot.yVarIdxs)                                                                 % loop through the y-variables to be reported
                        if ~isfield(cur_run.estORALS{pp},'kernel_err')
                            cur_run.estORALS{pp}.kernel_err = cur_run.estORALS{pp}.c_err;
                        end
                        ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp) = cur_run.estALS{pp}.(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)});
                        ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp) = cur_run.estORALS{pp}.(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)});
                        if strcmp(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)},'pathTestErr')                                  % normalize pathTestErr by dividing by L. MM:TBD: length of test paths should not need to be L
                            ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp) = ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp)/L{q}(idxs_1(k));
                            ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp) = ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp)/L{q}(idxs_1(k));
                        end
                        if strcmp(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)},'A_err')                                        % normalize A_err
                            ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp) = ys_ALS(x_idx_1,errVars_idx,k_ok,q,pp)/(N{q}(idxs_1(k))^2);
                            ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp) = ys_ORALS(x_idx_1,errVars_idx,k_ok,q,pp)/(N{q}(idxs_1(k))^2);
                        end
                    end
                    %fprintf('.')
                end
                expIdx_ys(x_idx_1,k_ok)=idxs_1(k);
                k_ok = k_ok + 1;
            end
        end
    end
end

ys_ALS = reshape(ys_ALS,[size(ys_ALS,1),size(ys_ALS,2),size(ys_ALS,3),size(ys_ALS,4)*size(ys_ALS,5)]);
ys_ORALS = reshape(ys_ORALS,[size(ys_ORALS,1),size(ys_ORALS,2),size(ys_ORALS,3),size(ys_ORALS,4)*size(ys_ORALS,5)]);

%% Plot
plotHandle = [];
plotFaceAlpha = 0.02;
figure;
for errVars_idx = 1:length(ToPlot.yVarIdxs)                                                                                     % loop through the y-variables to be reported
    loc_ys_mean             = squeeze(mean(ys_ALS(:,errVars_idx,:,:),4));                                                       % ALS results
    loc_ys_std              = squeeze(std(ys_ALS(:,errVars_idx,:,:),[],4));
    plotHandle(end+1)       = plot(log10(uniqxVarunion),log10(loc_ys_mean),'LineWidth',2,'LineStyle','-','Marker',ToPlot.yVarMarker{ToPlot.yVarIdxs(errVars_idx)},'MarkerFaceColor',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)},'Color',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)});
    hold on;
    patch(log10([uniqxVarunion;flipud(uniqxVarunion)]),log10([loc_ys_mean+loc_ys_std;flipud(loc_ys_mean-loc_ys_std)]),get(plotHandle(errVars_idx),'color'), 'FaceAlpha',plotFaceAlpha,'EdgeColor',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)});hold on;
    plot(log10(uniqxVarunion),log10(loc_ys_mean+loc_ys_std),'-','LineWidth',0.25,'color',get(plotHandle(errVars_idx),'color'));
    plot(log10(uniqxVarunion),log10(max(loc_ys_mean-loc_ys_std,0)),'-','LineWidth',0.25,'color',get(plotHandle(errVars_idx),'color'));
end
for errVars_idx = 1:length(ToPlot.yVarIdxs)                                                                                     % loop through the y-variables to be reported
    loc_ys_mean             = squeeze(mean(ys_ORALS(:,errVars_idx,:,:),4));                                                     % ORALS results
    loc_ys_std              = squeeze(std(ys_ORALS(:,errVars_idx,:,:),[],4));
    plotHandle(end+1)       = plot(log10(uniqxVarunion),log10(loc_ys_mean),'LineWidth',2,'color',get(plotHandle(errVars_idx),'color'),'LineStyle','-.','Marker',ToPlot.yVarMarker{ToPlot.yVarIdxs(errVars_idx)},'MarkerFaceColor',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)},'Color',ToPlot.yVarColor{ToPlot.yVarIdxs(errVars_idx)});
    patch(log10([uniqxVarunion;flipud(uniqxVarunion)]),log10([loc_ys_mean+loc_ys_std;flipud(loc_ys_mean-loc_ys_std)]),get(plotHandle(errVars_idx),'color'), 'FaceAlpha',plotFaceAlpha);hold on;
    plot(log10(uniqxVarunion),log10(loc_ys_mean+loc_ys_std),'-','LineWidth',0.25,'color',get(plotHandle(errVars_idx),'color'));
    plot(log10(uniqxVarunion),log10(max(loc_ys_mean-loc_ys_std,0)),'-','LineWidth',0.25,'color',get(plotHandle(errVars_idx),'color'));
end
axis tight;grid on
xlabel('log_{10}M','FontSize',14); ylabel('Estimation Error','FontSize',14);
legend_labels = {};
for k = 1:length(ToPlot.yVarIdxs), legend_labels{end+1} = [ToPlot.yVarLegend{ToPlot.yVarIdxs(k)}, ' (ALS)']; end
for k = 1:length(ToPlot.yVarIdxs), legend_labels{end+1} = [ToPlot.yVarLegend{ToPlot.yVarIdxs(k)}, ' (ORALS)']; end

curaxes = gca;
vertLine_x= (ToPlot.N_values^2+ToPlot.n_values)/(ToPlot.L_values*ToPlot.d_values*ToPlot.N_values);
line(log10([vertLine_x,vertLine_x]),(([curaxes.YLim(1),curaxes.YLim(2)])),'LineWidth',1,'color',[0.1,0.1,0.1],'LineStyle','--');
text(log10(vertLine_x),curaxes.YLim(1)+0.25,'$\frac{N^2+p}{LdN}$', "Interpreter", 'latex','FontSize',14);
%set(curaxes,'XTick',union(get(curaxes,'XTick'),log10(vertLine_x)));
%xTickLabels=get(curaxes,'XTickLabel'); xTickLabels{2}='$\frac{N^2+p}{LdN}$';set(gca,'XTickLabel',xTickLabels);set(curaxes,'TickLabelInterpreter','latex','FontSize',14);

vertLine_x = (ToPlot.N_values^2*ToPlot.n_values)/(ToPlot.L_values*ToPlot.d_values*ToPlot.N_values);
line(log10([vertLine_x,vertLine_x]),([curaxes.YLim(1),curaxes.YLim(2)]),'LineWidth',1,'color',[0.1,0.1,0.1],'LineStyle','--');
text(log10(vertLine_x),curaxes.YLim(1)+0.25,'$\frac{N^2p}{LdN}$', "Interpreter", 'latex','FontSize',14);
%set(curaxes,'XTick',union(get(gca,'XTick'),log10(vertLine_x)));
%xTickLabels=get(curaxes,'XTickLabel'); xTickLabels{5}='$\frac{N^2p}{LdN}$';set(gca,'XTickLabel',xTickLabels);set(curaxes,'TickLabelInterpreter','latex','FontSize',14);
%legend(plotHandle,legend_labels,'FontSize',12,'Location','NorthEast');
%gcfpos = get(gcf,'Position'); set(gcf,'Position',[gcfpos(1),gcfpos(2),gcfpos(3)*1.2,gcfpos(4)*1.2]);


leg=legend(plotHandle,legend_labels,'FontSize',12,'Location','NorthEast');
leg.Position(1)=leg.Position(1)+0.05;
leg.Position(2)=leg.Position(2)+0.05;
%gcfpos = get(gcf,'Position'); set(gcf,'Position',[gcfpos(1),gcfpos(2),gcfpos(3)*1.2,gcfpos(4)*1.2]);

%% Set up quantities for 2-D plots
ToPlot.xVarName1    = 'M';
ToPlot.xVarName2    = 'L';

ToPlot.yVarIdxs     = [1,2,4];
ToPlot.yVarName     = {'kernel_err','A_err','Z_err','pathTestErr'};

% Load and create table for 2-D plots
ys_ALS = [];

uniqxVarunion1 = []; uniqxVarunion2 = [];
for q = 1:length(loadDirs)
    eval(sprintf('xVar1{q}      = %s{q};',ToPlot.xVarName1));
    eval(sprintf('xVar2{q}      = %s{q};',ToPlot.xVarName2));
    eval(sprintf('uniqxVar1{q}  = uniq%s{q};',ToPlot.xVarName1));
    eval(sprintf('uniqxVar2{q}  = uniq%s{q};',ToPlot.xVarName2));
    uniqxVarunion1 = union(uniqxVarunion1,uniqxVar1{q})';
    uniqxVarunion2 = union(uniqxVarunion2,uniqxVar2{q})';
end

for q = 1:length(loadDirs)
    for x_idx_1 = 1:length(uniqxVarunion1)
        idxs_1  = find( xVar1{q}==uniqxVarunion1(x_idx_1) );
        if isempty(idxs_1), continue; end
        for x_idx_2 = 1:length(uniqxVarunion2)
            idxs_2  = find( xVar2{q}==uniqxVarunion2(x_idx_2) );
            idxs = intersect( idxs_1,idxs_2 );                                                                                  % indices of the experiments with the current values of uniqxVar1 and uniqxVar2
            for k = 1:length(idxs)                                                                                              % loop through all the experiments with the given value of the x,y variables
                uniqxVar1idx = find(uniqxVarunion1==xVar1{q}(idxs(k)));
                uniqxVar2idx = find(uniqxVarunion2==xVar2{q}(idxs(k)));
                if  ismember(n{q}(idxs(k)),ToPlot.n_values) && ...                                                              % check if the values of other variables satisfy the request
                        ismember(N{q}(idxs(k)),ToPlot.N_values) && ...
                        ismember(d{q}(idxs(k)),ToPlot.d_values) && ...
                        ismember(A_sparsity{q}(idxs(k)),ToPlot.A_sparsity_values) && ...
                        ismember(viscosity{q}(idxs(k)),ToPlot.viscosity_values) && ...
                        ismember(obs_std{q}(idxs(k)),ToPlot.obs_std_values)
                    try
                        load(sprintf('%s/LIKSonGraphs_largeTest_%d.mat',loadDirs{q},idxs(k)),'cur_run');
                    catch
                        continue
                    end

                    for pp = 1:length(cur_run.estALS)                                                                           % loop through the multiple runs of the same experiment
                        for errVars_idx = 1:length(ToPlot.yVarIdxs)                                                             % loop through the y-variables to be reported
                            ys_ALS(x_idx_1,x_idx_2,errVars_idx,q,pp) = cur_run.estALS{pp}.(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)});
                        end
                        if strcmp(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)},'pathTestErr')                                  % normalize pathTestErr by dividing by L. MM:TBD: length of test paths should not need to be L
                            ys_ALS(x_idx_1,x_idx_2,errVars_idx,q,pp) = ys_ALS(x_idx_1,x_idx_2,errVars_idx,q,pp)/L{q}(idxs(k));
                        end
                        if strcmp(ToPlot.yVarName{ToPlot.yVarIdxs(errVars_idx)},'A_err')                                        % normalize A_err
                            ys_ALS(x_idx_1,x_idx_2,errVars_idx,q,pp) = ys_ALS(x_idx_1,x_idx_2,errVars_idx,q,pp)/(N{q}(idxs(k))^2);
                        end
                        % fprintf('.')
                    end
                    expIdx_ys(x_idx_1,x_idx_2)=idxs(k);
                end
            end
        end
    end
end

ys_ALS = reshape(ys_ALS,[size(ys_ALS,1),size(ys_ALS,2),size(ys_ALS,3),size(ys_ALS,4)*size(ys_ALS,5)]);
ys_ORALS = reshape(ys_ORALS,[size(ys_ORALS,1),size(ys_ORALS,2),size(ys_ORALS,3),size(ys_ORALS,4)*size(ys_ORALS,5)]);

%% Plot
figure;
set(gcf,'Position',[1160 1230 1453 188]);
for errVars_idx = 1:length(ToPlot.yVarIdxs)                                                                                     % loop through the y-variables to be reported
    subplot(1,length(ToPlot.yVarIdxs),errVars_idx);
    loc_ys_mean             = squeeze(mean(ys_ALS(:,:,errVars_idx,:),4));
    loc_ys_std              = squeeze(std(ys_ALS(:,:,errVars_idx,:),[],4));
    imagesc(log10(flipud(loc_ys_mean')));
    xticks(1:length(uniqxVarunion1)); xticklabels(uniqxVarunion1);
    yticks(1:length(uniqxVarunion2)); yticklabels(flipud(uniqxVarunion2));
    colorbar; axis tight;
    xlabel(ToPlot.xVarName1,'FontSize',16); ylabel(ToPlot.xVarName2,'FontSize',16);
    title(ToPlot.yVarLegend(ToPlot.yVarIdxs(errVars_idx)),'FontSize',16);
    axis equal;
end
