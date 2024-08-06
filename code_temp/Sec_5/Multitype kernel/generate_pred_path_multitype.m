function pred_path = generate_pred_path_multitype(estResult, dyn_sys, initial)


dyn_sys.A                   = estResult.A;
% dyn_sys.phi_kernel          = estResult.phi_kernel;
dyn_sys.obs_std             = 0;
dyn_sys.phi_kernel_choices  = estResult.sorted_multi_kernel;
dyn_sys.kernel_idx          = estResult.sorted_id;

initial_type = class(initial);

switch initial_type
    case 'double'
        X0                          = initial;                                               
        pred_path                   = graph_forward_model_multitype( dyn_sys, X0, 1, 1);      
    case 'cell'
        M = length(initial);
        pred_path = cell(M, 1);
        parfor i = 1:M
            X0 = initial{i};
            pred_path{i}            = graph_forward_model_multitype( dyn_sys, X0, 0, 1);      
        end
end