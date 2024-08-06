function X0 = set_particle_initial_all_dim( N, d, InitialConditions )

% This code generate the initial distribution of the interacting particles
% To avoid difficulty in higher dimensions, we only take 3 essential distributions

% IN:
%   InitialConditions: can be a string:
%       N_mu_sigma:       Normal distribution with mean (mu, ..., mu) and Sigma = Id*sigma
%       DN_mu_sigma:      Average of two normal distributions with mean +-(mu,..., mu) and Sigma = Id*sigma
%       Unif_a_b:         Uniform distribution in the region [a,b]^d
%   If InitialConditions is an, it should be and N by d array and X0 will be simply equal to InitialConditions
%

% OUT:
%   X0 of size (N,d)
%

if isstring( InitialConditions ) || ischar( InitialConditions )                                                                 % If InitialConditions is a string, draw from the corresponding distribution
    Names = strsplit(InitialConditions,'_');
    switch Names{1}
        case 'N'      % normal(mu,sigma)
            mu = str2double(Names(2));
            sigma = str2double(Names(3));
            X0 = normrnd(mu,sigma,[N,d]);
        case 'DN'     % Double Normal: two normal pdfs with symmetric centers
            mu = str2double(Names(2));
            sigma = str2double(Names(3));
            N1 = floor(N/2);
            X0 = [normrnd(mu,sigma,[N1,d]); normrnd(-mu, sigma, [N-N1,d])];
            ind = randperm(N);
            X0 = X0(ind, :);
        case 'Unif'   % uniform
            a = str2double(Names(2));
            b = str2double(Names(3));
            X0 = unifrnd(a, b, N, d);
    end
elseif isfloat( InitialConditions )                                                                                             % If InitialConditions is an array, these are the initial conditions
    X0 = InitialConditions;
end

end