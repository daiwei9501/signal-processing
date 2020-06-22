%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deyu Meng, Fernando De la Torre. Robust matrix factorization %
% with unknown noise, ICCV, 2013                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo on data with mixture noise experiments                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:\Users\70916\Desktop\ICCV-MoG\Weighted L2 MF\Efficient')
clear;
clc;

m = 40;
n = 20;                     %Data size
r = 4;                      %Rank
num = 10;
for kk = 1:30
    kk
    RU = randn(m,r);
    RV = randn(r,n);
    X_Ori = RU * RV;        %Original data
    Ind = randperm(m*n);
    p1 = floor(m*n*0.2);
    W = ones(m,n);
    W(Ind(1:p1)) = 0;       %Indicator matrix
    X_Noi = X_Ori;
    X_Noi = W.*X_Noi;       %Add missing components
    p2 = floor(m*n*0.2);
    X_Noi(Ind(p1+1:p1+p2)) = X_Noi(Ind(p1+1:p1+p2)) + rand(1,p2)* 10 - 5; %Add uniform noise
    p3 = floor(m*n*0.2);
    X_Noi(Ind(p1+p2+1:p1+p2+p3)) = X_Noi(Ind(p1+p2+1:p1+p2+p3)) + randn(1,p3)* 0.2; %Add Gaussian noise
    X_Noi(Ind(p1+p2+p3+1:end)) = X_Noi(Ind(p1+p2+p3+1:end)) + randn(1,m*n-p1-p2-p3)* 0.01; %Add Gaussian noise
    
    aa = median(abs(X_Noi(Ind(p1+1:end))));
    aa = sqrt(aa/r);
    for i = 1:num
        U0 = rand(m,r)*aa*2-aa;
        V0 = rand(n,r)*aa*2-aa;
        
        %%%%%%%%%%%%%%%%%%GMM method %%%%%%%%%%%%%%%%%%%%%%%%
        param.maxiter = 100;
        param.OriX = X_Ori;
        param.InU = U0;
        param.InV = V0;
        param.k = 3;
        param.display = 0;
        param.NumIter = 100;
        param.tol = 1.0e-50;
        param.method = 2;
        [label, model,TW, A,B,llh] =  MLGMDN(W,X_Noi,r,param);
        E1G(kk,i) = sum(sum(abs(W.*(X_Noi-A*B'))));
        E2G(kk,i) = sum(sum((W.*(X_Noi - A*B')).^2));
        E3G(kk,i) = sum(sum(((X_Ori - A*B')).^2));
        E4G(kk,i) = sum(sum(abs((X_Ori-A*B'))));
        E5G(kk,i) = subspace(RU,A);
        E6G(kk,i) = subspace(RV',B);
        SG(kk,i) = llh(end);
        
    end
end

for kk = 1:30
    %%%%%%%%%%%%%%%%%%GMM method %%%%%%%%%%%%%%%%%%%%%%%%
    [a ii] = max(SG(kk,1:num));
    EE1G(kk) = E1G(kk,ii);
    EE2G(kk) = E2G(kk,ii);
    EE3G(kk) = E3G(kk,ii);
    EE4G(kk) = E4G(kk,ii);
    EE5G(kk) = E5G(kk,ii);
    EE6G(kk) = E6G(kk,ii);
end

i = 0;
figure;
i = i + 1;
subplot(2,3,i);
bar(EE1G);
axis([0.3 30.7 0 max(EE1G)*1.1]);
title('|W*(X_{Noise}-UV^T)|_1');
i = i + 1;
subplot(2,3,i);
bar(EE2G);
axis([0.3 30.7 0 max(EE2G)*1.1]);
title('|W*(X_{Noise}-UV^T)|_2');
i = i + 1;
subplot(2,3,i);
bar(EE3G);
axis([0.3 30.7 0 max(EE3G)*1.1]);
title('|X_{Clean}-UV^T|_2');
i = i + 1;
subplot(2,3,i);
bar(EE4G);
axis([0.3 30.7 0 max(EE4G)*1.1]);
title('|X_{Clean}-UV^T|_1');
i = i + 1;
subplot(2,3,i);
bar(EE5G);
axis([0.3 30.7 0 max(EE5G)*1.1]);
title('Subspace(U,U_{Groundtruth})');
i = i + 1;
subplot(2,3,i);
bar(EE6G);
axis([0.3 30.7 0 max(EE6G)*1.1]);
title('Subspace(V,V_{Groundtruth})');

