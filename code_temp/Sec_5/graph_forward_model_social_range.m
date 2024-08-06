function xpath = graph_forward_model_social_range(I, X0, progressON, rk4ON )
% solve dynamical system
%% load parameters
steps = I.steps;
d = I.d;
N = I.N;
dt = I.dt;
v = I.viscosity;
E = I.E;

id_long = I.id_long;
id_short = I.id_short;

if nargin == 2
    progressON = 1;
end

if nargin == 3
    rk4ON = 0;
end

noise = sqrt(2*v)/sqrt(dt)* randn(N, d, steps);

phi_kernel_long = I.phi_kernel_long;
phi_kernel_short = I.phi_kernel_short;

phi_kernel_long_for_PS = @(x) (x~=0).*phi_kernel_long(x)./abs(x);
phi_kernel_short_for_PS = @(x) (x~=0).*phi_kernel_short(x)./abs(x);

% myRK4 integration

RHS   = @(x) RHSfn(x, E, phi_kernel_long_for_PS, phi_kernel_short_for_PS, id_long, id_short);
xpath = myrk4(X0, dt, RHS, steps, noise, progressON, rk4ON);  %


end







function dx = RHSfn(x, E, kernel_long, kernel_short, id_long, id_short)
% The RHS function of the ODE
[N, d] = size(x);
% y = reshape(x,d,N); % [x_1^T \vdots x_N^T]  d x N
dx = zeros(N, d);

for k = 1:length(id_long)
    i = id_long(k);
    
    dif  = x(i, :) - x;
    dis  = sqrt(sum(dif.^2, 2));
    temp = kernel_long(dis).*dif;
    temp(isnan(temp)) = 0;
    %     aa(:, i) = sum(temp, 2);
    dx(i, :) = sum(temp.*E(:, i), 1);
end

for k = 1:length(id_short)
    i = id_short(k);
    
    dif  = x(i, :) - x;
    dis  = sqrt(sum(dif.^2, 2));
    temp = kernel_short(dis).*dif;
    temp(isnan(temp)) = 0;
    %     aa(:, i) = sum(temp, 2);
    dx(i, :) = sum(temp.*E(:, i), 1);
end



% for i = 1:N
%     dif  = x(i, :) - x;
%     dis  = sqrt(sum(dif.^2, 2));
%     temp = kernel(dis).*dif;
%     temp(isnan(temp)) = 0;
%     %     aa(:, i) = sum(temp, 2);
%     dx(i, :) = sum(temp.*E(:, i), 1);
% end
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
    deltax = RK4(x0,dt,RHS,0, rk4ON);
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
if rk4ON
    k1 = RHS_x(x0)             + force;
    k2 = RHS_x(x0+1/2.*dt.*k1) + force;
    k3 = RHS_x(x0+1/2.*dt.*k2) + force;
    k4 = RHS_x(x0+dt.*k3)      + force;
    deltax= 1/6*(k1+2*k2+2*k3+k4);
else
    deltax = RHS_x(x0);
end
end
