function [err_relmeanL2,estPathObj,meanL2traj] = getPathTestError( dyn_sys, learning_set, A_hat, c_hat, testingPathsObj )

err_relmeanL2   = Inf;
estPathObj       = [];
meanL2traj      = Inf;

if ~isempty( testingPathsObj.paths )
    %dyn_sys.viscosity           = 0;                                                                                            % No point in adding stochastic forcing term to test path (? this is not quite right, the situation seems complex, and what I think it the right thing is computationally expensive, so settling for this now)
    dyn_sys.A                   = A_hat;
    [dyn_sys.phi_kernel,...
       dyn_sys.phi_kernel_cheb] = get_kernel_from_c( c_hat, learning_set.dict );
    dyn_sys.phi_kernel          = dyn_sys.phi_kernel_cheb;
    dyn_sys.initial             = getPathsNoObsNoise( testingPathsObj );                                                        % The initial conditions are the same as those of the true paths
    dyn_sys.obs_std             = 0;
    estPathObj                  = get_paths( dyn_sys, length(testingPathsObj.paths), 'ParforProgressON', 0, ...
                                                'saveON', 0, 'loadON', 0, 'forcing_noise', testingPathsObj.forcing_noise );     % Compare with "true" paths, with the same forcing noise
    [err_relmeanL2,meanL2traj]  = traj_err( testingPathsObj.paths, estPathObj.paths );
end

end