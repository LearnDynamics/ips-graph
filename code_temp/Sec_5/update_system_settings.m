function I = update_system_settings(I)
% This code update the system settings give the fundmental ones

%% Update the basis functions for kernel
switch I.basis_case
    case 1         % % % % dictionary1
        p = 8;q = 2;cut = 0.5;
        I.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
        I.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
        I.dict{3} = @(x) abs(x)<=cut;
        I.dict{4} = @(x) 0*x+1;
        I.dict{5} = @(x) x.^2;
        I.c = [-0.3333, 1.3333, -160, 3, -2]';
    case 2           % % % % dictionary2
        I.dict{1} = @(x) x.^2;
        I.dict{2} = @(x) abs(x);
        I.dict{3} = @(x) cos(x);
        I.c = [.1, 0.1, 3]';
    case 3           % % % dictionary3
        n = 10;
        dict = cell(n, 1);
        for i = 1:n
            dict{i} = @(x) sin(x*i + i);
        end
        I.dict = dict;
        I.c =  ones(n, 1);% [.1, 0.1, 3]';
    case 4         % % % % dictionary1
        p = 2;q = .2;cut = 0.5;
        I.dict{1} = @(x)  x.^(-p-1).*(abs(x)>=cut);
        I.dict{2} = @(x) - x.^(-q-1).*(abs(x)>=cut);
        I.dict{3} = @(x) - 1.*(abs(x)<=cut);
        c1 = -1; c2= 1;
        c3 = (c1*I.dict{1}(cut) + c2*I.dict{2}(cut))/I.dict{3}(cut);
        I.c = [c1, c2, c3]';
    case 5
        I.dict{1} = @(x) x.^2;
        I.dict{2} = @(x) x.^4;
        I.dict{3} = @(x) 0*x+1;
        I.c = [-10, -10, -0.1]';
    case 6
        p = 8;q = 2;cut = 0.5;
        I.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
        I.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
        I.dict{3} = @(x) abs(x)<=cut;
        I.c = [-0.3333, 1.3333, -160]';
    case 7         % % % % dictionary1
        p = 8;q = 2;cut = 0.5;
        I.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
        I.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
        I.dict{3} = @(x) abs(x)<=cut;
        I.dict{4} = @(x) sin(x);
        I.dict{5} = @(x) x.^2;
        
        for l = 6:20
            I.dict{l} = @(x) sin(x/(l-4));
        end
        
        I.c(1:5) = [-0.3333, 1.3333, -160, 10, -0.5];    %
        I.c(6:20) = randn(15, 1)*0.1;
        %         I.c = I.c';
    case 8           % % % dictionary3
        n = 6;
        dict = cell(n, 1);
        for i = 1:n
            dict{i} = @(x) sin(x*i + i);
        end
        I.dict = dict;
        I.c =  ones(n, 1);% [.1, 0.1, 3]';
    case 'LF'
        n = 2;
        dict = cell(n, 1);
        threshold = 1;
        dict{1} = @(x) (x<=threshold);
        dict{2} = @(x) (x<=threshold*1.5).*(x>threshold);
        I.dict = dict;
        I.c = [-1, -0.1]';
    case 'randomSmoothFourier'  
        try I = rmfield(I,'dict'); catch end
        try I = rmfield(I,'c_true'); catch end
        n = 20;
        I.c = zeros(n,1);
        for k = 1:n
            I.dict{k}     = @(x) sin(2*pi*x*k);
            I.c(k)  = 1/(k)*randn(1,1);
        end
    case 'randomSmoothFourierWithDecay'  
        try I = rmfield(I,'dict'); catch end
        try I = rmfield(I,'c_true'); catch end
        n = 20;
        I.c = zeros(n,1);
        for k = 1:n
            I.dict{k}     = @(x) sin(2*pi*x*k)./(x+eps);
            I.c(k)  = 1/(k)*randn(1,1);
        end
end

I.n = length(I.dict);       % Dimension of the hypothesis space
I.phi_kernel = get_kernel_from_c(I.c, I.dict);

%% Special examples

if strcmp(I.basis_case, 'Kuramoto')
    %     basis_choice = 'Hermite';
    basis_choice = 'trigonometric';
    switch basis_choice
        case 'Hermite'
            syms x
            n = 10;
            dict_sym = hermiteH([0:1:n-1], x);
            
            dict = cell(n, 1);
            dict{1} = @(x) 0*x + 1;
            for i = 2:n
                dict{i} = matlabFunction(dict_sym(i));
            end
        case 'trigonometric'
            n = 7;
            dict = cell(2*n+1, 1);
            dict{1} = @(x) cos(x);
            for i = 1:n
                dict{2*i} = @(x) sin((i+1)*x);
                dict{2*i+1} = @(x) cos((i+1)*x);
            end
    end
    I.dict = dict;
    I.n = length(I.dict);       % Dimension of the hypothesis space
    n = I.n;
    I.c = ones(I.n, 1);% [.1, 0.1, 3]';
    
    I.phi_kernel = @(x) sin(x);
    
    % find the truc coefficient according to the given basis
    xgrid = 0:0.01:5;
    N = length(xgrid);
    b = sin(xgrid);
    A = zeros(N, n);
    for i = 1:I.n
        A(:, i) = dict{i}(xgrid);
    end
    I.c = A\b';
end

if strcmp(I.basis_case, 'social_range')
    %     lb = 0;
    %     rb = 2;
    %     knot_num = 10;
    %     deg = 3;
    %     free_bdry = 1;
    %     short_ratio = 0.3;
    %
    %     [basis, knots] = spline_basis_MATALB(lb, rb, knot_num, deg, free_bdry);
    %     I.dict = basis;
    %     %     I.c_short = [1, 1, 1, 1, 1, 0.4, 0.2, 0, 0, 0, 0, 0];
    %     %     I.c_long = [-10, -10, -10, 0, 0, 2, 1, 1, 0, 0, 0, 0];
    %
    %     I.c_short = [1, 1, 1, 1, 1, 0.4, -10, -10, -10, -10, 0, 0];
    %     I.c_long = [-10, -10, -10, 0, 0, 2, 1, 1, 0, 0, 0, 0];
    
    
    lb = 0;
    rb = 2;
    knot_num = 5;
    deg = 2;
    free_bdry = 1;
    short_ratio = 0.3;
    
    [basis, knots] = spline_basis_MATALB(lb, rb, knot_num, deg, free_bdry);
    I.dict = basis;
    
    I.c_short = [2, 2, 1, 0, 0, 0];
    I.c_long = [0, 0, 0, 1, 2, 2];
    
    
    
    
    
    
    
    I.n = length(I.dict);
    I.phi_kernel_short = get_kernel_from_c(I.c_short, I.dict);
    I.phi_kernel_long  = get_kernel_from_c(I.c_long, I.dict);
    
    plotON = 1;
    if plotON
        figure;
        fplot(I.phi_kernel_long, [0, 10]);
        hold on;
        fplot(I.phi_kernel_short, [0, 10]);
        legend('long', 'short')
    end
    %     I.phi_kernel = I.phi_kernel_long;
    %     I.c = I.c_long';
    
    id = randperm(I.N);
    short_num = floor(I.N*short_ratio);
    I.id_short = id(1:short_num);
    I.id_long = id(short_num+1:end);
    
    I.coef_mat = zeros(I.N, I.n);
    for j = 1:length(I.id_short)
        i = I.id_short(j);
        I.coef_mat(i, :) = I.c_short;
    end
    for j = 1:length(I.id_long)
        i = I.id_long(j);
        I.coef_mat(i, :) = I.c_long;
    end
    
    [U, S, V] = svd(I.coef_mat);
    
    UU = U(:, 1:2);
    VV = V(:, 1:2);
    SS = S(1:2, 1:2);
    
    I.UU = UU;
    I.VV = VV*SS;
end


%% other updates
I.TN = I.steps + 1;                                         % Total time steps
I.tgrid = I.t0:I.dt:(I.steps)*I.dt;
I.X0 = set_particle_initial_all_dim(I.N, I.d, I.initial);   % Initial condition
I.Z_true = get_Z_from_E_c(I.A, I.c);     % Z is the product of E and c
I.graph_norm  = norm(I.A,'fro');
%% dictionary for the kernel
% p = 8;q = 2;cut = 0.5;
% I.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
% I.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
% I.dict{3} = @(x) abs(x)<=cut;
% I.dict{4} = @(x) sin(x);
% I.dict{5} = @(x) x.^2;
% I.c = [-0.3333, 1.3333, -160, 10, -0.5]';


%% dictionary for the kernel
% % % First setting, 20 basis
% p = 8;q = 2;cut = 0.5;
% I.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
% I.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
% I.dict{3} = @(x) abs(x)<=cut;
% I.dict{4} = @(x) sin(x);
% I.dict{5} = @(x) x.^2;
%
% for l = 6:20
%     I.dict{l} = @(x) sin((l-4)*x);
% end
%
% I.c(1:5) = [-0.3333, 1.3333, -160, 10, -0.5];    %
% I.c(6:20) = randn(15, 1)*0.1;
% I.c = I.c';

%% dictionary for the kernel
% Second setting

% I.dict{1} = @(x) x.^2;
% I.dict{2} = @(x) abs(x);
% I.dict{3} = @(x) cos(x);
% I.c = [1, 0.1, 3]';

%% dictionary for the kernel

% I.dict{1} = @(x) x.^2;
% I.dict{2} = @(x) abs(x);
% I.dict{3} = @(x) (x-1).^2;
% I.c = [-1, 0.1, -3]';

%% dictionary for the kernel
% Second setting

% I.dict{1} = @(x) x.^2;
% I.dict{2} = @(x) abs(x);
% I.dict{3} = @(x) cos(x);



%% ===== This is a choice of basis
% n = 3;
% dict = cell(n, 1);
% for i = 1:n
%     dict{i} = @(x) sin(x*i + i);
% end
% I.dict = dict;
% I.c = ones(n, 1);
%






%% Old stuff
% I.dict_mat = zeros(I.n, I.n);
% for i = 1:I.n
%     for j = 1:I.n
%         I.dict_mat(i, j) = integral(@(x) I.dict{i}(x).*I.dict{j}(x), 0, I.L);
%     end
% end
% I.nlfn = 'LJnew_8_2'; %***********************
% [I.phi_kernel, I.Phi_potential] = get_nlfn(I.nlfn, I.L, I.dx);
%% update sysInfo settings
% I.M = 2*I.L/I.dx;                     % The number of space mesh for each dimension
% I.mid = floor(I.L/I.dx) + 1;
% I.T = I.dt * I.steps;

% I.rgrid = 0:I.dx:I.L;

% xgrid = -I.L:I.dx:I.L;
% [I.XX, I.YY] = meshgrid(xgrid, xgrid);


% [I.phi_kernel, I.Phi_potential] = get_nlfn(I.nlfn, I.L, I.dx);

% [sysInfo.U0, sysInfo.X0] = set_Initial(sysInfo.L, sysInfo.dx, sysInfo.N, sysInfo.initial, sysInfo.d);

%% oberved particle index
% particleIndex   = 1:sysInfo.observation_gap:sysInfo.N;      % index of observed particles
% N = I.N;
% gap = I.observation_gap;
%
% I.K = floor(N/gap);
% particleIndex = randsample(I.N, I.K);
% d   = I.d;
%
% xindex   = zeros(1,d*I.K);
% for i=1:I.K
%     ind1 = d*(i-1)+1:d*i;
%     ind2 = (particleIndex(i)-1)*d + (1:d);
%     xindex(ind1) = ind2;
% end
% I.xindex = xindex;
%% RBM
% if isfield(I,'RBM')
%     [I.U0, I.X0] = set_Initial(I.L, I.dx, I.RBM.N, I.initial, I.d);
%     [~, I.RBM.X0_traj] = set_Initial(I.L, I.dx, I.RBM.K, I.initial, I.d);
% end
end