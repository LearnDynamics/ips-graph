function E = reform_E(E)
% The problem we previous solved are N linear systems of N-1 parameters
% We have to transform the result into the actual adjacent matrix 
% Simply putting 0 diagonals.

N = length(E);
U = triu(E,1);
D = E - U;

E = [U;zeros(1, N)] + [zeros(1, N);D];
end