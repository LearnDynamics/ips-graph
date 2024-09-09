function a_test_Btensor(dyn_sys,learning_setup,Btmp,dX,n,N,B,testN)
%% test if dX /dt == a*B*c
A     = dyn_sys.A; a1 = A(1,:);  a1 = A(:,1)';
ctrue = learning_setup.c;

switch testN
    case 1    % B1 that uses 
        B1    = Btmp;
        Bdx   = 0; for k=1:n; Bdx  = Bdx + B1(:,:,:, :,k)*ctrue(k); end
        Bdx1  =0;  for i=1:N; Bdx1 = Bdx1+ a1(i)*Bdx(i,:,:,:); end
        norm(Bdx1 - dX(N,:,:,:),'fro')
        figure;
        plot(squeeze(Bdx1(1,1,1,:)),'LineWidth',1); hold on;
        plot(squeeze(dX(N,1,1,:)),'--x'); legend('aBc','dX')
    case 2 % B{i} has size N x dLM x n; dX has size N x dLM
        B1 =B{1};
        Bdx   = 0; for k=1:n; Bdx  = Bdx + B1(:,:,k)*ctrue(k); end
        Bdx1  = a1*Bdx;
        norm(Bdx1 - dX(1,:),'fro')
        figure;
        plot(Bdx1(1,:),'LineWidth',1); hold on;
        plot(dX(1,:),'--x'); legend('aBc','dX')
        mean(Bdx1(1,:)./dX(1,:))  
    case 3
        B1 =B{1};
        Bdx   = 0; for k=1:n; Bdx  = Bdx + B1(:,:,:, :,k)*ctrue(k); end
        Bdx1  =0;  for i=1:N; Bdx1 = Bdx1+ a1(i)*Bdx(i,:,:,:); end
        norm(Bdx1 - dX(1,:,:,:),'fro')
        figure;
        plot(squeeze(Bdx1(1,1,1,:)),'LineWidth',1); hold on;
        plot(squeeze(dX(1,1,1,:)),'--x'); legend('aBc','dX')
end


end