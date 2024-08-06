function [err_relmeanL2, meanL2traj, err_meanL2, err_relL2] = traj_err( traj_true, traj2_est )

err_relL2   = zeros(length(traj_true),1);
err_L2      = zeros(length(traj_true),1);
for m = 1:length(traj_true)
    L2sqtraj(m)     = sum((traj_true{m}).^2,'all');
    err_relL2(m)    = sqrt(sum((traj_true{m}-traj2_est{m}).^2,'all')/L2sqtraj(m));
    err_L2(m)       = sqrt(sum((traj_true{m}-traj2_est{m}).^2,'all'));
end

meanL2traj      = mean(sqrt(L2sqtraj));
err_relmeanL2   = mean( err_relL2 );
err_meanL2      = mean(err_L2);
end