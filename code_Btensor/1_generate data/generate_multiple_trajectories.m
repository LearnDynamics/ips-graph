function all_xpath = generate_multiple_trajectories(I, M)
all_xpath = cell(M, 1);
for i = 1:M
    I.X0 = set_particle_initial_all_dim(I.N, I.d, I.initial);
    all_xpath{i} = graph_forward_model(I, 0);
end
end