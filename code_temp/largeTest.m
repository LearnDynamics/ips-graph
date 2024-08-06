% Make sure to run using parallelism
delete(gcp('nocreate'));
gcp;

% Data where to save the experiment's results; if it exists and non-empty, it will load the parameters for that experiment
% and continue it
expSetup.saveDir = '~/DataAnalyses/LIKSonGraphs_try';

% If true, removes existing directory, its experimental setup, and its results. Overrides conflicting options that may follows (e.g. COMPLETE_PARTIAL_RESULTS)
expSetup.START_ANEW = true;

% When true, if saveDir includes a file with the partial results of an experiment, tries to re-do the experiment; if false,
% it will skip that experiment
expSetup.COMPLETE_PARTIAL_RESULTS = true;

if expSetup.START_ANEW
    system(sprintf('rm -rf %s',expSetup.saveDir));
    expSetup.COMPLETE_PARTIAL_RESULTS = false;
end

if ~exist(expSetup.saveDir,"dir"), mkdir(expSetup.saveDir); end

%% Initialize and randomize
if ~exist(sprintf('%s/LIKSonGraphs_largeTest_setup.mat',expSetup.saveDir),'file')
    expSetup.M_test          = 50;
    expSetup.ntrials         = 8;
    expSetup.runALS          = true;
    expSetup.runORALS        = true;                                                                                            % use B-tensor and long matrix
    expSetup.runORSVD        = false;
    expSetup.runORALS_normalMat  = false;                                                                                        % use normal matrix in ORALS --- the original approach
    expSetup.saveALSstats    = false;                                                                                           % ALS stats during iterations

    fprintf('\nGenerating parameters...')
    params{1}.Name  = 'N';          params{1}.Values = {4,8,32,64};
    params{2}.Name  = 'd';          params{2}.Values = {1};
    params{3}.Name  = 'viscosity';  params{3}.Values = {1e-4,1e-2};
    params{4}.Name  = 'A_sparsity'; params{4}.Values = {0.1,0.8};
    params{5}.Name  = 'L';          params{5}.Values = {2,4,8,16};   % MM: dt is fixed, and we let T=dt*L
    params{6}.Name  = 'obs_std';    params{6}.Values = {1e-4,1e-2};
    params{7}.Name  = 'M';          params{7}.Values = {16,32,64,256,1024,4096,8192};
    params{8}.Name  = 'n';          params{8}.Values = {4,8,16,32};

    paramStructs    = GenerateStructuresWithVariedParameters( params );
    rng(1);
    expSetup.p_perm = randperm(length(paramStructs));
    save(sprintf('%s/LIKSonGraphs_largeTest_setup.mat',expSetup.saveDir),'params','paramStructs','expSetup');
else
    fprintf('\nLoading existing sets of parameters...')
    load(sprintf('%s/LIKSonGraphs_largeTest_setup.mat',expSetup.saveDir));
end

runExperimentBatch( expSetup, paramStructs );

fprintf('\n');