clc
clear
close all
addpath('.\Gabor')
addpath('.\RPCA')
%% simulate the signal
L = 1e4;        % signal length
T = 300;        % mean period of impulses

impulses = zeros(L,1);
I = 0:T:L;
I = round(I + 1e-2*T*randn(size(I)));
I = I(I<=L);
I = I(I>0);
impulses(I) = sign(randn(size(I))) + .1*randn(size(I));

r = 0.95;
a = [1 -2*cos(2*pi*.2)*r r^2];
b = [1 -1];
x = filter(b,a,impulses);
% noise = randn(L,1);
% NSR = 3;  % NSR = 3 is close to working limit...
SNR =  -5; %dB
% signal = x + NSR*noise*std(x);
power_signal =  std(x); % std(x,1) normalizes by N and produces the square root of the second
   %  moment of the sample about its mean
%power_signal = std(x)^2;
%power_noise = (power_signal/(10^(0.1*SNR)))^0.5
%power_signal1 = sqrt(sum(x.^2-mean(x))/L);
power_noise=power_signal/( 10^(0.1*SNR));
signal  = x  + randn(L,1)*power_noise;

figure(1);
plot(signal,'b');
hold on
plot(x,'g');
%%
% Gabor analysis without zero-padding
Nw = 2^6;
R = Nw*2/8;
% R = Nw*1/8;
g = hanning(Nw,'periodic')/sqrt(1.5);
c = cgt(signal,g,R);

figure(2);
imagesc(abs(c))
%% construct the rank-1 matrix
M = abs(c).^2;

figure(3);
imagesc(M)

% rank = 1;
% lambda = 10; power =0;
% [LR,SR,RMSE,error]=SSGoDec(M,rank,lambda,power);

% LR = LR.';
% SR = SR.';
% Accelerated Proximal Gradient: APG
% RPCA for low rank and sparse matrix seperation
% LR is the low rank part and SR is sparse part

 % lambda = 0.01;
 % mu = 10e-3;
 % [LR, SR] = proximal_gradient_rpca(M, lambda,10000,1e-7,0,1,0.9,mu);
  %lambda = 0.01;
  %[LR, SR] = proximal_gradient_rpca(M, lambda);

D_measured = M;
Li = 0.5;
N_iter = 10;
SC = 1e-2;
mu = norm(D_measured,'fro'); mu_final = 1e-16*mu; alp = 0.02;
[LR,SR,para_val]  = RPCA_FISTA(Li,N_iter,D_measured,SC,mu_final,mu,alp);


%% plot the reconstructed sparse and low rank matrix
figure(4)
imagesc(10*log10(abs(LR)))
title('Low rank decomposition L','FontSize',14)
%clim = get(gca, 'clim');
cmap=colormap('hot');
colormap(cmap(end:-1:1, :));
% set(gca, 'clim',[0 max(max(max(LR)),max(max(SR)))]);
set(gca, 'clim', max(max(max(10*log10(abs(LR)))),max(max(10*log10(abs(SR))))) + [-30 0])
xlabel('Time(seconds)','FontSize',14)
ylabel('Frequency(Hz)','FontSize',14)

figure(5)
imagesc(10*log10(abs(SR)))
set(gca, 'clim', max(max(max(10*log10(abs(LR)))),max(max(10*log10(abs(SR)))))+ [-30 0])
%imagesc(x_axis,[],SR)
%set(gca, 'clim',[0 max(max(max(LR)),max(max(SR)))]);
%set(gca, 'clim',clim);
cmap=colormap('hot');
colormap(cmap(end:-1:1, :));
title('Sparse decomposition S (lambda = 0.01)','FontSize',14);
xlabel('Time(seconds)','FontSize',14);
ylabel('Frequency(Hz)','FontSize',14);

%% reconstruct the complex low rank and sparse matrix
PHASE = angle(c);
LRC = sqrt(LR).*exp(1i.*PHASE);
SRC = sqrt(SR).*exp(1i.*PHASE);

% refine the results and filtering
gain = 6;
mask = double(abs(SRC)> (gain*abs(LRC))); 
SRC = SRC.*mask;

% plot the scree diagram
D1 = svd(LR);
D2 = svd(SR);
D3 = svd(SRC);
D4 = svd(LRC);

figure(6);
subplot(411);
stem(D1);
xlabel('sigular values of low rank');
subplot(412);
stem(D2);
xlabel('sigular values of real sparse matrix');
subplot(413);
stem(D3);
xlabel('sigular values of complex sparse matrix after masking');
subplot(414);
stem(D4);
xlabel('sigular values of complex low rank');
%
%%
figure(7)
imagesc(abs(SRC));
cmap=colormap('hot');
colormap(cmap(end:-1:1, :));
% set(gca, 'clim',[0 max(max(max(abs(LRC))),max(max(abs(SRC))))]);
title('Sparse decomposition X (Nw = 256 lambda = 0.02)','FontSize',14)
xlabel('Time(seconds)','FontSize',14)
ylabel('Frequency(Hz)','FontSize',14)

% inverse Gabor transform
[xx_SRC,gg] = icrgt(SRC,g,R);

figure(8)
plot(signal)
hold on
plot(x,'g'),
hold on
plot(real(xx_SRC),'r')
title('Nw = 256  lambda = 0.4','FontSize',14)
%legend('Measurements','Transient Signal','Separated Transient Signal','fontsize',16)
xlabel('Sample Number','FontSize',14)
ylabel('Amplitude','FontSize',14)


