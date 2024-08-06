function err = graph_err(A, dyn_system)

% Returns the relative Frobenius norm error between the adjacency matrix A and the matrix dyn_system.A

err = norm(A - dyn_system.A, 'fro')./norm(dyn_system.A,'fro');

end