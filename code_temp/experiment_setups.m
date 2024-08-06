switch experimentOption
    case 1
        %% Good example with rather small random graph, small M, moderate viscosity and observation noise, good ALS performance
        M                       = 16;                                                                                           % number of trajectories
        dyn_sys.L               = 10;                                                                                           % number observations equispaced in time [0,T]

        dyn_sys.N               = 10;                                                                                           % number of agents
        dyn_sys.d               = 2;                                                                                            % dim of state vectors
        dyn_sys.viscosity       = 1e-3;                                                                                         % viscosity (forcing noise in the dynamics)
        dyn_sys.G               = createGraph('random_sparse',struct('N',dyn_sys.N, 'sparsity', 0.4, 'normalized', 1,'selfLoops',0) ); % create influence graph
        dyn_sys.A               = dyn_sys.G.W;
        dyn_sys.initial         = 'Unif_0_5';
        dyn_sys.T               = 1e-1;
        dyn_sys.dt              = dyn_sys.T/dyn_sys.L;                                                                          % this gets reduced (within graph_forward_model) to T/L so that L observations fit in 0:dt:T
        dyn_sys.obs_std         = 1e-3;                                                                                         % observation noise
        dyn_sys                 = system_settings( dyn_sys );                                                                   % settings of the IPS and graph and its integrator

        learning_set            = learning_settings( 'randomSmoothFourierWithDecay' );
        learning_set.Z_true     = get_Z_from_E_c( dyn_sys.A, learning_set.c );                                                  % Z is the product of E and c
        dyn_sys.phi_kernel      = learning_set.phi_kernel;                                                                      % MM: logical inconsistency here and in the following because of this
    case 2
        %% Good example with a circle graph, small M, moderate viscosity and observation noise, good ALS performance
        M                       = 16;                                                                                           % number of trajectories; phase transition for ORALS is between 16 and 32
        dyn_sys.L               = 10;                                                                                           % number observations equispaced in time [0,T]

        dyn_sys.N               = 10;                                                                                           % number of agents
        dyn_sys.d               = 2;                                                                                            % dim of state vectors
        dyn_sys.viscosity       = 1e-4;                                                                                         % viscosity (forcing noise in the dynamics)
        dyn_sys.G               = createGraph('circle',struct('N',dyn_sys.N, 'normalized', 1,'selfLoops',0) );                  % create influence graph
        dyn_sys.A               = dyn_sys.G.W;
        dyn_sys.initial         = 'Unif_0_5';
        dyn_sys.T               = 1e-1;
        dyn_sys.dt              = dyn_sys.T/dyn_sys.L;                                                                          % this gets reduced (within graph_forward_model) to T/L so that L observations fit in 0:dt:T
        dyn_sys.obs_std         = 1e-4;                                                                                         % observation noise
        dyn_sys                 = system_settings( dyn_sys );                                                                   % settings of the IPS and graph and its integrator

        learning_set            = learning_settings( 'randomSmoothFourierWithDecay' );
        learning_set.Z_true     = get_Z_from_E_c( dyn_sys.A, learning_set.c );                                                  % Z is the product of E and c
        dyn_sys.phi_kernel      = learning_set.phi_kernel;                                                                      % MM: logical inconsistency here and in the following because of this    
    case 3
        %% Good example with a circle graph, with sparse-in-time observations small M, moderate viscosity and observation noise
        %  L is 10 times smaller than the number of steps simulated in [0,T], yet ALS recovers a good graph, but the kernel is (of course!) off
        %  This is pretty robust to noise (both viscosity and in the obs). Increasing M increases the path prediction accuracy,
        %   the graph estimation accuracy by ORALS, but not the kernel accuracy (of course!)
        M                       = 64;                                                                                           % number of trajectories; phase transition for ORALS is between 16 and 32
        dyn_sys.L               = 40;                                                                                           % number observations equispaced in time [0,T]

        dyn_sys.N               = 10;                                                                                           % number of agents
        dyn_sys.d               = 2;                                                                                            % dim of state vectors
        dyn_sys.viscosity       = 1e-4;                                                                                         % viscosity (forcing noise in the dynamics)
        dyn_sys.G               = createGraph('circle',struct('N',dyn_sys.N, 'normalized', 1,'selfLoops',0) );                  % create influence graph
        dyn_sys.A               = dyn_sys.G.W;
        dyn_sys.initial         = 'Unif_0_5';
        dyn_sys.T               = 1e-1;
        dyn_sys.obs_std         = 1e-4;                                                                                         % observation noise
        dyn_sys                 = system_settings( dyn_sys );                                                                   % settings of the IPS and graph and its integrator

        learning_set            = learning_settings( 'randomSmoothFourierWithDecay' );
        learning_set.Z_true     = get_Z_from_E_c( dyn_sys.A, learning_set.c );                                                  % Z is the product of E and c
        dyn_sys.phi_kernel      = learning_set.phi_kernel;                                                                      % MM: logical inconsistency here and in the following because of this
    case 4
        %% Example with a circle graph, with small M, large L and T moderate viscosity and observation noise
        %  L is 10 times smaller than the number of steps simulated in [0,T], yet ALS recovers a good graph, but the kernel is (of course!) off
        %  This is pretty robust to noise (both viscosity and in the obs). Increasing M increases the path prediction accuracy,
        %   the graph estimation accuracy by ORALS, but not the kernel accuracy (of course!)
        %  Increasing M is quite interesting, as is decreasing L.
        M                       = 32;                                                                                           % number of trajectories; phase transition for ORALS is between 16 and 32
        dyn_sys.L               = 64;                                                                                           % number observations equispaced in time [0,T]

        dyn_sys.N               = 10;                                                                                           % number of agents
        dyn_sys.d               = 2;                                                                                            % dim of state vectors
        dyn_sys.viscosity       = 1e-2;                                                                                         % viscosity (forcing noise in the dynamics)
        dyn_sys.G               = createGraph('circle',struct('N',dyn_sys.N, 'normalized', 1,'selfLoops',0) );                  % create influence graph
        dyn_sys.A               = dyn_sys.G.W;
        dyn_sys.initial         = 'Unif_0_5';
        dyn_sys.dt              = 1e-1;                                                                                         % this gets reduced (within graph_forward_model) to T/L if needed so that L observations fit in 0:dt:T
        dyn_sys.T               = dyn_sys.dt*dyn_sys.L;
        dyn_sys.obs_std         = 1e-10;                                                                                         % observation noise
        dyn_sys                 = system_settings( dyn_sys );                                                                   % settings of the IPS and graph and its integrator

        learning_set            = learning_settings( 'randomSmoothFourierWithDecay' );
        learning_set.Z_true     = get_Z_from_E_c( dyn_sys.A, learning_set.c );                                                  % Z is the product of E and c
        dyn_sys.phi_kernel      = learning_set.phi_kernel;                                                                      % MM: logical inconsistency here and in the following because of this
    case 5
        %% Example with a circle graph, with small M, large L and T moderate viscosity and observation noise
        % ORALS does not seem to converge as L grows (?)
        M                       = 32;                                                                                           % number of trajectories; phase transition for ORALS is between 16 and 32
        dyn_sys.L               = 2048;                                                                                           % number observations equispaced in time [0,T]
        n                       = 40;
        dyn_sys.N               = 40;                                                                                           % number of agents
        dyn_sys.d               = 2;                                                                                            % dim of state vectors
        
        dyn_sys.viscosity       = 1e-2;                                                                                         % viscosity (forcing noise in the dynamics)
        dyn_sys.G               = createGraph('circle',struct('N',dyn_sys.N, 'normalized', 1,'selfLoops',0) );                  % create influence graph
        dyn_sys.A               = dyn_sys.G.W;
        dyn_sys.initial         = 'Unif_0_5';
        dyn_sys.dt              = 1e-1;                                                                                         % this gets reduced (within graph_forward_model) to T/L if needed so that L observations fit in 0:dt:T
        dyn_sys.T               = dyn_sys.dt*dyn_sys.L;
        dyn_sys.obs_std         = 1e-10;                                                                                         % observation noise
        dyn_sys                 = system_settings( dyn_sys );                                                                   % settings of the IPS and graph and its integrator

        learning_set            = learning_settings( 'randomSmoothFourierWithDecay', dyn_sys, struct('n',n) );
        learning_set.Z_true     = get_Z_from_E_c( dyn_sys.A, learning_set.c );                                                  % Z is the product of E and c
        dyn_sys.phi_kernel      = learning_set.phi_kernel;                                                                      % MM: logical inconsistency here and in the following because of this

end


