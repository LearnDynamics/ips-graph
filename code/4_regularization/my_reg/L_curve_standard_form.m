function [lambda_opt, c] = L_curve_standard_form(A, b, plotON)
%{
In this code, we solve the following standard form optimization
        c'A*c - 2b'*c + b'A^{-1}b + lambda c'c
We used explicit computation of the L-curvature from Hansen's code. In his
setting (classical LSE), the goal is minimize

|| A1c - b1 ||^2 + lambda^2 ||c||^2

We need to pay attention that
1. the parameter lambda has ^2 in Hansen's setting, but we only take lambda
2. we need to rephrase all the matrices in order to use this code.

We could recover the same setting as Hansen's, with new parameter matrices

A1 = sqrt(A)
b1 = sqrt(A^{-1})b

by noting that  c'*A1'*A1*c - 2c'*A1'*b1 , thus, A1'*A1= A; A1'*b1= b;
2022-10-03 Quanjun Lang
%}

% initialize
N  = 100;
lb = -16;
rb = 2;

lambda_seq = 10.^(linspace(lb, rb, N));
E = zeros(N, 1); R = E;
%% Spectral Analysis of A
[U, S, ~] = svd(A);
s = diag(S);
b_coord = U'*b;
s2 = sqrt(s);



for ll = 1:N
    lambda = lambda_seq(ll);
    c_coord = b_coord./(s + lambda);
    E(ll) = norm(s2.*c_coord - b_coord./s2);
    R(ll) = norm(c_coord);
end

% Compute g = - curvature of L-curve.
beta = b_coord./s2;
xi = b_coord./s;
ss = s2;



g = lcfun(sqrt(lambda_seq),ss,beta,xi);
% Locate the corner.  If the curvature is negative everywhere,
% then define the leftmost point of the L-curve as the corner.
[~,gi] = min(g);
reg_c = fminbnd('lcfun',sqrt(lambda_seq(max(gi-1,1))),sqrt(lambda_seq(min(gi+1,N))),...
    optimset('Display','on'),ss,beta,xi); % Minimizer.

kappa_max = - lcfun(reg_c,ss,beta,xi); % Maximum curvature.
reg_c = reg_c^2;
if (kappa_max < 0)
    reg_c = lambda_seq(1); E_opt = E(1); R_opt = R(1);
    c_coord = b_coord./(s + reg_c);
else
    c_coord = b_coord./(s + reg_c);
    E_opt = norm(s2.*c_coord - b_coord./s2);
    R_opt = norm(c_coord);
end


lambda_opt = reg_c;
c = U*c_coord;

%%
if plotON
    figure;
    subplot(321);plot(log10(lambda_seq), -g, 'LineWidth', 3);hold on;plot(log10(reg_c), kappa_max, '+', 'LineWidth', 5);grid on
    legend(['\kappa_{max} = ', num2str(kappa_max)],['\lambda_{opt} = ', num2str(reg_c)],'Location','best')
    title('curvature \kappa of the L-curve');xlabel('log_{10}(\lambda)');ylabel('Curvature \kappa')
    
    subplot(322);scatter(log10(E), log10(R));hold on;plot(log10(E_opt), log10(R_opt), '+', 'LineWidth', 5);grid on
    title('L-curve');xlabel('E');ylabel('R');
    legend('L-curve',['\lambda_{opt} = ', num2str(reg_c)])
    subplot(323);scatter(log10(lambda_seq), log10(E));hold on;xlabel('\lambda');ylabel('E');plot(log10(lambda_opt), log10(E_opt), '+', 'LineWidth', 5);
    subplot(324);scatter(log10(lambda_seq), log10(R));hold on;xlabel('\lambda');ylabel('R');plot(log10(lambda_opt), log10(R_opt), '+', 'LineWidth', 5);
    subplot(325);plot(log10(s),'o');hold on;plot(log10(abs(U'*b)),'.');title('Picard Ratio');legend('\sigma_i', 'u_i^T b')
    set(gcf,'Position',[100 100 800 550])
end


if false
    %% Generate figure for the paper
    f = figure;
    g1 = subplot(1,2,1);
    line1 = plot(log10(E), log10(R), 'o');hold on;
    line2 = plot(log10(E_opt), log10(R_opt), '+', 'LineWidth', 5);grid on
    title('L-curve');
%     xlabel('$\text{log}_{10}(\mathcal{E})$', 'Interpreter', 'latex');
    xlabel('$\log_{10}\mathcal{E}$', 'Interpreter', 'latex');
    ylabel('$\log_{10}\mathcal{R}$', 'Interpreter', 'latex');
%     y_t = -10:2:10;
%     y_label = cell(length(y_t), 1);
%     for i = 1:length(y_t)
%         y_label{i} = ['10^{', num2str(y_t(i)), '}'];
%     end
%     yticks(-10:2:10);
%     yticklabels(y_label);
    
    g2 = subplot(1,2,2);
    line3 = plot(log10(lambda_seq), -g, 'LineWidth', 3);hold on;
    line4 = plot(log10(reg_c), kappa_max, '+', 'LineWidth', 5);grid on
    xlim([lb, rb]);
    %     legend(['\kappa_{max} = ', num2str(kappa_max)],['\lambda_{opt} = ', num2str(reg_c)],'Location','best')
    % legend(['\kappa_{max} = ', num2str(kappa_max)], 'Location','best')
    title('Curvature');
    xlabel('$\log_{10} \lambda$', 'Interpreter', 'latex');
    ylabel('Curvature \kappa')
    
    
    %     legend('L-curve',['\lambda_{opt} = ', num2str(reg_c)])
    %     subplot(223);scatter(log10(lambda_seq), log10(E));hold on;
    %     xlim([lb, rb]);
    %     xlabel('log_{10}(\lambda)');;ylabel('E');plot(log10(lambda_opt), log10(E_opt), '+', 'LineWidth', 5);
    %     subplot(224);scatter(log10(lambda_seq), log10(R));hold on;
    %     xlim([lb, rb]);
    %     xlabel('log_{10}(\lambda)');;ylabel('R');plot(log10(lambda_opt), log10(R_opt), '+', 'LineWidth', 5);
    % %     subplot(325);plot(log10(s),'o');hold on;plot(log10(abs(U'*b)),'.');title('Picard Ratio');legend('\sigma_i', 'u_i^T b')
    %     set(gcf,'Position',[100 100 800 550])
    
    
    % Create a tile on the right column to get its position
    %     ax = subplot(1,3,3,'Visible','off');
    %     axPos = ax.Position;
    %     delete(ax)
    % Construct a Legend with the data from the sub-plots
    %     lg = legend([line1,line2,line3]);
    % Move the legend to the position of the extra axes
    % lg.Position(1:2) = axPos(1:2);
    
    
    
    lg = legend([line1,line2,line3, line4], {'L-curve', ['\lambda_{opt} = ', num2str(reg_c, '%.4e')], 'curvature', ['\kappa_{max} = ', num2str(kappa_max, '%.4e')]}, 'Orientation','horizontal','Box','off');
    lg.Position(1:2) = [0.16, 0.014];
    tightfig;
    set_positionFontsAll;
    % f.PaperPosition = [0, 0, 24, 10]
    pos=get(g2,'position');  % retrieve the current values
    pos(2)=1.49*pos(2);        % try reducing width 10%
    pos(4)=0.91*pos(4);        % try reducing width 10%
    set(g2, 'position', pos);
    %
    pos=get(g1,'position');  % retrieve the current values
    pos(2)=1.49*pos(2);        % try reducing width 10%
    pos(4)=0.91*pos(4);        % try reducing width 10%
    set(g1, 'position', pos);
    
    
    
end
end

