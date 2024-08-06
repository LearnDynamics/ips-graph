function [xpath,forcing_noise,xpath_allt,Lidxs] = graph_forward_model( dyn_sys, X0, progressON, rk4ON, forcing_noise )

% Generate a path of the graph IPS, on the time interval [0,T], with stepsize dt, observed at the times 0:T/L:T
%
% OUT:
%   xpath           : N x d x L path for N agents, each d-dimensional, observed at the times 0:T/L:T
%   forcing_noise   : the Wiener path in the stochastic forcing term
%   xpathallt       : N x d x T/dt path for N agents, each d-dimensional, observed at the times 0:dt:T
%   Lidxs           : L vector of indices into the time dimension of xpathallt indicating when the observations of xpath are (in xpathallt)
%

%% load parameters
T       = dyn_sys.T;
dt      = dyn_sys.dt;
L       = dyn_sys.L;
d       = dyn_sys.d;
N       = dyn_sys.N;
v       = dyn_sys.viscosity;
A       = dyn_sys.A;

if nargin < 3,      progressON = 1;         end
if nargin < 4,      rk4ON = 0;              end
if nargin < 5,      forcing_noise = [];     end

tN      = L; % floor(T/dt);    % this tN = floor(T/dt) is totally unexpected. It leads to dX != aBc. Spent several days to figure it out.  
if tN < L                                                                                             % If dt is too large, yielding tN<L, reduce dt so that we can ensure L observations. This choice of action is debatable, but it's chosen to make the choice of L observations a priority
    dt  = T/L;
    tN  = floor(T/dt);
end
if isempty( forcing_noise ) 
    forcing_noise       = sqrt(2*v)/sqrt(dt)* randn(N, d, tN);
end

phi_kernel          = dyn_sys.phi_kernel;
phi_kernel_for_PS   = @(x) (x~=0).*phi_kernel(x)./abs(x);

RHS         = @(x) RHSfn(x, A, phi_kernel_for_PS);                                                                              % right-hand side of the system of stochastic ODEs
xpath_allt  = myrk4( X0, dt, RHS, tN, forcing_noise, progressON, rk4ON );                                                       % integrated the system

Lidxs = unique(round(linspace(1,tN,L)));

xpath = xpath_allt(:,:,Lidxs);

end

% The RHS function of the ODE
function dx = RHSfn(x, E, kernel)

[N, d]  = size(x);
dx      = zeros(N, d);

for i = 1:N
    dif                 = x(i, :) - x;                                     % N x d 
    dis                 = sqrt(sum(dif.^2, 2));  
    Kphi                = kernel(dis).*dif;                                 % Note kernel is phi(x)./abs(x)
    Kphi(isnan(Kphi))   = 0;
    dx(i, :)            = E(:,i)'* Kphi; % sum(Kphi.*E(:, i), 1);         % it uses E(:,i)!!! Different from paper that use E(i,:)
                                                                           % Also, Kphi*E(:,i) with size(Khpi) = N x d,   
end

end


%%
function x = myrk4(x0, dt, RHS, tN, noise, progressON, rk4ON)

[N, d] = size(x0);
if ~exist('noise','var')
    noise = zeros(N, d, tN);
end

x = zeros(N, d, tN+1);
x(:, :, 1) = x0;

if progressON;reverseStr = '';end

for t = 1:tN
    % deltax = RK4(x0,dt,RHS,0, rk4ON);
    deltax = RHS(x0); 
    x0     = x0 + dt* (deltax+ noise(:, :, t));
    x(:, :, t+1) = x0;
    if progressON
        percentDone = 100 * t / tN;
        msg = sprintf('Generating particles : %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end
if progressON
    fprintf('\n')
end
end


function deltax = RK4(x0,dt,RHS_x,force, rk4ON)
if rk4ON  % This RK4 is not good for SDEs  
    k1 = RHS_x(x0)             + force;
    k2 = RHS_x(x0+1/2.*dt.*k1) + force;
    k3 = RHS_x(x0+1/2.*dt.*k2) + force;
    k4 = RHS_x(x0+dt.*k3)      + force;
    deltax= 1/6*(k1+2*k2+2*k3+k4);
else
    deltax = RHS_x(x0);
end
end


