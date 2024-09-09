function learning_setup = learning_settings( basisType, dyn_sys, opts )

% Settings for learning the interaction kernels

if isnumeric(basisType)
    fprintf('Parametric inference: kernel-type = %d \n', basisType); 
else 
    fprintf('Nonparametric inference: kernel-type = %s\n', basisType); 
end

%% Set basis functions for kernel estimation
switch basisType
    case 1         % % % % dictionary1
        p = 8;q = 2;cut = 0.5;
        learning_setup.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
        learning_setup.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
        learning_setup.dict{3} = @(x) abs(x)<=cut;
        learning_setup.dict{4} = @(x) 0*x+1;
        learning_setup.dict{5} = @(x) x.^2;
        learning_setup.c = [-0.3333, 1.3333, -160, 3, -2]';
    case 2           % % % % dictionary2
        learning_setup.dict{1} = @(x) x.^2;
        learning_setup.dict{2} = @(x) abs(x);
        learning_setup.dict{3} = @(x) cos(x);
        learning_setup.c = [.1, 0.1, 3]';
    case 3           % % % dictionary3
        n = 10;
        dict = cell(n, 1);
        for i = 1:n
            dict{i} = @(x) sin(x*i + i)*2^(i-2);
        end
        learning_setup.dict = dict;
        learning_setup.c =  ones(n, 1);% [.1, 0.1, 3]';
    case 4         % % % % dictionary1
        p = 2;q = .2;cut = 0.5;
        learning_setup.dict{1} = @(x)  x.^(-p-1).*(abs(x)>=cut);
        learning_setup.dict{2} = @(x) - x.^(-q-1).*(abs(x)>=cut);
        learning_setup.dict{3} = @(x) - 1.*(abs(x)<=cut);
        c1 = -1; c2= 1;
        c3 = (c1*learning_setup.dict{1}(cut) + c2*learning_setup.dict{2}(cut))/learning_setup.dict{3}(cut);
        learning_setup.c = [c1, c2, c3]';
    case 5
        learning_setup.dict{1} = @(x) x.^2;
        learning_setup.dict{2} = @(x) x.^4;
        learning_setup.dict{3} = @(x) 0*x+1;
        learning_setup.c = [-10, -10, -0.1]';
    case 6  % LJ model 
        p = 8;q = 2;cut = 0.5;
        learning_setup.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
        learning_setup.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
        learning_setup.dict{3} = @(x) abs(x)<=cut;
        learning_setup.c = [-0.3333, 1.3333, -160]';
    case 61  % same as 
        p = 8;q = 2;cut = 0.5; 
        learning_setup.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
        learning_setup.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
        learning_setup.dict{3} = @(x) abs(x)<=cut;
        learning_setup.c = [-0.3333, 1.3333, -160]';
    case 'typical_example_Lenard_Jones'
        % p = 8;q = 2;cut = 0.5;
        learning_setup.dict = cell(10, 1);
        for k = 0:2
            learning_setup.dict{1+k} =  @(x) x.^(-9).*(abs(x)>0.5 + 0.25*k);
        end

        for k = 0:2
            learning_setup.dict{4+k} =  @(x) x.^(-3).*(abs(x)>0.5 + 0.25*k);
        end

        for k = 0:3
            learning_setup.dict{7+k} =  @(x) abs(x)<=(0.5 + 0.25*k);
        end

        learning_setup.c = zeros(10, 1);
        learning_setup.c(1) = -1/3;
        learning_setup.c(4) = 4/3;
        learning_setup.c(7) = -160;
        %% Test the shape of basis
        % xgrid = linspace(0, 3, 1000)';
        % for k = 1:10
        %     y(:, k) = learning_setup.dict{k}(xgrid);
        % end
        % figure;subplot(131)
        % plot(xgrid, log10(y(:, 1:4)+0.001));legend('1','2','3','4')
        % subplot(132)
        % plot(xgrid, log10(y(:, 5:7)+0.001));legend('5','6','7')
        % subplot(133)
        % plot(xgrid, log10(y(:, 8:10)+0.001));legend('8','9','10')

        
    case 7         % % % % dictionary1
        p = 8;q = 2;cut = 0.5;
        learning_setup.dict{1} = @(x) x.^(-p-1).*(abs(x)>cut);
        learning_setup.dict{2} = @(x) x.^(-q-1).*(abs(x)>cut);
        learning_setup.dict{3} = @(x) abs(x)<=cut;
        learning_setup.dict{4} = @(x) sin(x);
        learning_setup.dict{5} = @(x) x.^2;

        for l = 6:12
            learning_setup.dict{l} = @(x) sin(x*(l-4));
        end
        learning_setup.c = zeros(12,1); 
        learning_setup.c(1:5) = [-0.3333, 1.3333, -160, 10, -0.5]';    %
        learning_setup.c(6:12) = randn(7, 1)*1;
        %         I.c = I.c';
    case 8           % % % dictionary3
        n = 6;
        dict = cell(n, 1);
        for i = 1:n
            dict{i} = @(x) sin(x*i + i);
        end
        learning_setup.dict = dict;
        learning_setup.c =  ones(n, 1);% [.1, 0.1, 3]';
    case 'LF'
        n = 2;
        dict = cell(n, 1);
        threshold = 1;
        dict{1} = @(x) (x<=threshold);
        dict{2} = @(x) (x<=threshold*1.5).*(x>threshold);
        learning_setup.dict = dict;
        learning_setup.c = [-1, -0.1]';
    case 'randomSmoothFourier'
        try learning_setup = rmfield(learning_setup,'dict'); catch end
        try learning_setup = rmfield(learning_setup,'c_true'); catch end
        n = 20;
        learning_setup.c = zeros(n,1);
        for k = 1:n
            learning_setup.dict{k}  = @(x) sin(2*pi*x*k);
            learning_setup.c(k)     = 1/(k)*randn(1,1);
        end
    case 'randomSmoothFourierWithDecay'
        try learning_setup = rmfield(learning_setup,'dict'); catch,     end                                                     % clear up existing fields that need to be create anew in a moment
        try learning_setup = rmfield(learning_setup,'c_true'); catch,   end
        if nargin<3 || isempty(opts) || ~isfield(opts,'n')
            opts.n = 20;
        end
        learning_setup.c = zeros(opts.n,1);
        for k = opts.n:-1:1
            learning_setup.dict{k}  = @(x) sin(2*pi*x*k)./(x+0.1);
            learning_setup.c(k)     = 1/(k)*randn(1,1);
        end
end

learning_setup.kernelSupp   = [0,10];                                                                                           %MM:TBD:this should done better, system-specific, etc...
 try
     learning_setup.n             = length(learning_setup.dict);                                                                     % Dimension of the hypothesis space
 %   [learning_setup.phi_kernel,learning_setup.phi_kernel_cheb] = get_kernel_from_c(learning_setup.c, learning_setup.dict, learning_setup.kernelSupp);
     learning_setup.phi_kernel = get_kernel_from_c(learning_setup.c, learning_setup.dict, learning_setup.kernelSupp);
 end



%% Special examples
if strcmp(basisType, 'Kuramoto')
    %     basis_choice = 'Hermite';
    basis_choice = 'trigonometric';


    basis_choice = opts;
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
        case 'inHypo_Space'
            n = 7;
            dict = cell(2*n+2, 1);
            dict{1} = @(x) cos(x);
            for i = 1:n
                dict{2*i} = @(x) sin((i+1)*x);
                dict{2*i+1} = @(x) cos((i+1)*x);
            end
            dict{2*n+2} = @(x) sin(x); 
    end
    learning_setup.dict = dict;
    learning_setup.n = length(learning_setup.dict);       % Dimension of the hypothesis space
    n = learning_setup.n;
    learning_setup.c = ones(learning_setup.n, 1);% [.1, 0.1, 3]';

    learning_setup.phi_kernel = @(x) sin(x);

    % find the truc coefficient according to the given basis
    xgrid = 0:0.01:5;
    N = length(xgrid);
    b = sin(xgrid);
    A = zeros(N, n);
    for i = 1:learning_setup.n
        A(:, i) = dict{i}(xgrid);
    end
    learning_setup.c = A\b';
end

if strcmp(basisType, 'multitype')
    basis_choice = opts;
    % basis_choice = 'spline_long_short';


    % First we determine the basis for inference 

    switch basis_choice
        case 'spline_long_short'      % two kernels, long range and short range
            lb = 0;rb = 5;knot_num = 15;deg = 2;free_bdry = 1;
            
            [basis, ~] = spline_basis_MATALB(lb, rb, knot_num, deg, free_bdry);
            learning_setup.dict = basis;

            learning_setup.c_choices{1} = [0, -1, -2, -2, -1, 0, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0]';
            learning_setup.c_choices{2} = -[0, 0, 0, 0, 0, 0, -1, -2, -2, -1, 0, 1, 2, 2, 1, 0]';

            % learning_setup.c_choices{1} = [3, 3, 3, 2, 2, 1, 1, 0, 0, 0, 0]';
            % learning_setup.c_choices{2} = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3]';

            % learning_setup.c_choices{1} = [3, 1, 0]';
            % learning_setup.c_choices{2} = [0, 1, 2]';


            [learning_setup.phi_kernel_choices{1}, ~] = get_kernel_from_c(learning_setup.c_choices{1}, learning_setup.dict);
            [learning_setup.phi_kernel_choices{2}, ~]  = get_kernel_from_c(learning_setup.c_choices{2}, learning_setup.dict);

            type_weight = [1, 1]; % This is the weight of different kernel types. Will normalize later so it's easy to use.
        case 'fake_multitype'      % two kernels, long range and short range
            lb = 0;rb = 5;knot_num = 15;deg = 2;free_bdry = 1;

            [basis, ~] = spline_basis_MATALB(lb, rb, knot_num, deg, free_bdry);
            learning_setup.dict = basis;

            learning_setup.c_choices{1} = [0, -1, -2, -2, -1, 0, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0]';
            learning_setup.c_choices{2} = -[0, 0, 0, 0, 0, 0, -1, -2, -2, -1, 0, 1, 2, 2, 1, 0]';

            % learning_setup.c_choices{1} = [3, 3, 3, 2, 2, 1, 1, 0, 0, 0, 0]';
            % learning_setup.c_choices{2} = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3]';

            % learning_setup.c_choices{1} = [3, 1, 0]';
            % learning_setup.c_choices{2} = [0, 1, 2]';


            [learning_setup.phi_kernel_choices{1}, ~] = get_kernel_from_c(learning_setup.c_choices{1}, learning_setup.dict);
            [learning_setup.phi_kernel_choices{2}, ~]  = get_kernel_from_c(learning_setup.c_choices{2}, learning_setup.dict);

            type_weight = [1, 0]; % This is the weight of different kernel types. Will normalize later so it's easy to use.
        case 'Fourier_random'
            try learning_setup = rmfield(learning_setup,'dict'); catch end
            try learning_setup = rmfield(learning_setup,'c_true'); catch end

            n = 20;
            type_nums = 3;

            learning_setup.c = zeros(n,1);
            for k = 1:n
                learning_setup.dict{k}  = @(x) sin(2*pi*x*k);
                learning_setup.c(k)     = 1/(k)*randn(1,1);
            end
    

            for i = 1:type_nums
                learning_setup.c_choices{i} = randn(1, n);
                [learning_setup.phi_kernel_choices{i}, ~] = get_kernel_from_c(learning_setup.c_choices{i}, learning_setup.dict);
            end
            
            type_weight = rand(type_nums, 1); % This is the weight of different kernel types. Will normalize later so it's easy to use. 
    end
    learning_setup.n = length(learning_setup.dict);
    learning_setup.num_kernel_choices = length(learning_setup.phi_kernel_choices);
    
    
    plotON = 0;         % plot them if necessary
    if plotON   
        figure;hold on;
        for i = 1:learning_setup.num_kernel_choices
            fplot(learning_setup.phi_kernel_choices{i}, [0, 5], 'LineWidth',2,'DisplayName',['Kernel No.', num2str(i)]);
        end
        legend();
    end



    % Secondly we distribute those kernels to particles                                                     
    type_weight_normalized = type_weight./sum(type_weight);     % Will normalize later so it's easy to use. 
    kernel_idx = randsample(1:learning_setup.num_kernel_choices, dyn_sys.N, true, type_weight_normalized);
    
    coef_mat = zeros(learning_setup.n, dyn_sys.N);
    for i = 1:dyn_sys.N
        coef_mat(:, i) = learning_setup.c_choices{kernel_idx(i)};
    end
    

    learning_setup.kernel_idx = kernel_idx;
    learning_setup.coef_mat = coef_mat;

    [U, S, V] = svd(learning_setup.coef_mat);

    UU = U(:, 1:2);
    VV = V(:, 1:2);
    SS = S(1:2, 1:2);

    learning_setup.u = UU*SS;         % Let norm(UU) = 1 to fix the scale
    learning_setup.v = VV;

    % id = randperm(dyn_sys.N);
    % 
    % short_num = floor(dyn_sys.N*short_ratio);
    % learning_setup.id_short = id(1:short_num);
    % learning_setup.id_long = id(short_num+1:end);
    % 
    % 
    % 
    % id = randperm(dyn_sys.N);
    % short_num = floor(dyn_sys.N*short_ratio);
    % learning_setup.id_short = id(1:short_num);
    % learning_setup.id_long = id(short_num+1:end);
    % 
    % learning_setup.coef_mat = zeros(dyn_sys.N, learning_setup.n);
    % for j = 1:length(learning_setup.id_short)
    %     i = learning_setup.id_short(j);
    %     learning_setup.coef_mat(i, :) = learning_setup.c_short;
    % end
    % for j = 1:length(learning_setup.id_long)
    %     i = learning_setup.id_long(j);
    %     learning_setup.coef_mat(i, :) = learning_setup.c_long;
    % end
    % 
    % [U, S, V] = svd(learning_setup.coef_mat);
    % 
    % UU = U(:, 1:2);
    % VV = V(:, 1:2);
    % SS = S(1:2, 1:2);
    % 
    % learning_setup.UU = UU;
    % learning_setup.VV = VV*SS;
end

return
