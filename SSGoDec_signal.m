clc
clear
close all
addpath('.\Gabor');
addpath('.\GoDec\SSGoDec');
%% simulate the signal
L = 1e4;        % signal length
T = 600;        % mean period of impulses
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
SNR =  -10; %dB
power_signal =  std(x)^2;

%A = rand(L,1)*3;
%k = 1;
%noise = zeros(L,1);
%for i =1:L 
%   noise(i,1) = (randn(1)+k)*(A(i,1)<=1) +  (randn(1)-k)*(A(i,1)<=2 && A(i,1)>2 )+(randn(1))*(A(i,1)<=3 && A(i,1)>3) ;
%end
power_noise=power_signal/( 10^(0.1*SNR));
%signal = x + noise *power_noise;
signal  = x  + randn(L,1)*power_noise^0.5;
t = linspace(0,1,10000);
%x = x*1e3;
%signal = signal*1e3;
%%
figure;
subplot(211)
plot(t,signal,'b');
xlabel('Time (s)','FontSize',12)
ylabel('Amplitude','FontSize',12)
subplot(212)
plot(t,x,'g');
xlabel('Time (s)','FontSize',12)
ylabel('Amplitude','FontSize',12)
%%
Fs = 1e4;
Nw = 50;
%window = gausswin(256);
window = hanning(Nw,'periodic')/sqrt(1.5);
%window = hanning(Nw,'periodic');
%%同步压缩变换
[c,f] = fsst(signal,Fs,window);
% Gabor 
%R = Nw/8;
%c = cgt(signal,window,R);
%%
% 显示分析数据和基本的频谱分析
%Fs = 1e4;
%Nw =  Fs;
%Nv =  ceil(3/4*Nw);
%Nfft = Nw;
%[S,f] = pwelch(signal,Nw,Nv,Nfft,Fs);
%figure,
%plot(f,10*log10(S)) % 只显示到6000Hz
%xlabel('频率 (赫兹)'),ylabel('dB'),
%title(['信号频谱 ; 频率分辨率 = ',num2str(diff(f(1:2))),'Hz'])
%%
%sax1 = hilbert(signal);
%envelopofsax1 = abs(sax1).^2;
%E1 = fft(envelopofsax1);
%N = 1e4;
%figure,
%subplot(211)
%plot(abs(E1)/N);
%xlim([0 600])
%subplot(212)
%plot(abs(E2)/N);
%xlim([0 600])
%% construct the rank-1 matrix 
M = abs(c).^2;

D = svd(M);
figure;
stem(D);
%M = abs(Sx).^2;
%%
%D_measured = M;
%Li = 0.5;
%N_iter = 20;
%SC = 1e-2;
%mu = norm(D_measured,'fro'); mu_final = 1e-16*mu; alp = 0.02;
%[LR,SR,para_val]  = RPCA_FISTA(Li,N_iter,D_measured,SC,mu_final,mu,alp);
%SSGoDec
tau = 0.02;
power = 6;
rank = 1;
[LR,SR,RMSE,error] = SSGoDec(M,rank,tau,power);
%%
figure,
%imagesc(t,f,10*log10(abs(M)));
imagesc(10*log10(abs(M)));
xlabel('Time (s)','FontSize',14)
ylabel('Frequency(Hz)','FontSize',14)
%colormap(1-gray)
cmap=colormap('hot');
colormap(cmap(end:-1:1, :));
%set(gca, 'clim',[0 max(max(max(LR)),max(max(SR)))]);
set(gca, 'clim', max(max(max(10*log10(abs(LR)))),max(max(10*log10(abs(SR))))) + [-30 0])
colorbar()
%% plot the reconstructed sparse and low rank matrix
figure,
%subplot(211)
imagesc(t,f,10*log10(abs(LR)))
clim = get(gca, 'clim');
cmap=colormap('hot');
colormap(cmap(end:-1:1, :));
%set(gca, 'clim',[0 max(max(max(LR)),max(max(SR)))]);
set(gca, 'clim', max(max(max(10*log10(abs(LR)))),max(max(10*log10(abs(SR))))) + [-30 0])
colorbar()
%colormap(1-gray)
%title('Low rank decomposition L','FontSize',14)
xlabel('Time (s)','FontSize',14)
ylabel('Frequency(Hz)','FontSize',14)
title('L')

figure,
%subplot(212)
%imagesc(abs(SR))
imagesc(t,f,10*log10(abs(SR)))
cmap=colormap('hot');
colormap(cmap(end:-1:1, :));
%set(gca, 'clim',[0 max(max(max(LR)),max(max(SR)))]);
set(gca, 'clim', max(max(max(10*log10(abs(LR)))),max(max(10*log10(abs(SR))))) + [-30 0])
colorbar()
%colormap(1-gray)
%title('Sparse decomposition S','FontSize',14);
xlabel('Time (s)','FontSize',14)
ylabel('Frequency(Hz)','FontSize',14)
title('S')
%%
% plot the scree diagram
D1 = svd(LR);
D2 = svd(SR);
figure;
subplot(211);
stem(D1);
xlabel('sigular values of low rank matrix');
subplot(212);
stem(D2);
xlabel('sigular values of sparse matrix');

%% reconstruct the complex low rank and sparse matrix
PHASE = angle(c);
LRC = sqrt(LR).*exp(1i.*PHASE);
SRC = sqrt(SR).*exp(1i.*PHASE);

% refine the results and filtering
%gain = 0.1;
%mask = double(abs(SRC)< (gain*abs(LRC))); 
Gain = 13.5;
%mask = double(abs(LRC)>Gain*mean(mean(abs(LRC))));
mask = double(abs(LRC)>Gain*mean(mean(abs(LRC))));
LRC = LRC.*mask;
%%
figure
imagesc(t,f,10*log10(abs(LRC)));
colormap(cmap(end:-1:1, :));
colorbar()
xlabel('Time(seconds)','FontSize',14)
ylabel('Frequency(Hz)','FontSize',14)
%set(gca, 'clim',[0 max(max(max(LR)),max(max(SR)))]);
%colormap(1-gray);
%set(gca, 'clim',[0 max(max(max(abs(LRC))),max(max(abs(SRC))))]);
set(gca, 'clim', max(max(max(10*log10(abs(LR)))),max(max(10*log10(abs(SR)))))+ [-30 0])
%title('Sparse decomposition LRC','FontSize',14)
%xlabel('Time(seconds)','FontSize',14)
%ylabel('Frequency(Hz)','FontSize',14)

%还原低秩部分
xx_LRC = ifsst(LRC,window);

%Gabor
%[xx_LRC,gg] = icrgt(LRC,window,R);
%%
figure
plot(t,signal)
hold on
plot(t,x,'g'),
hold on
plot(t,xx_LRC,'r')
xlabel('Time','FontSize',12)
ylabel('Amplitude','FontSize',12)

figure,
subplot(311)
plot(t,signal,'b');
xlabel('Time (s)','FontSize',12)
ylabel('Amplitude','FontSize',12)
%ylim([-2,2])
%text(0,2,8,'(A)','Color','k','FontSize',16)
subplot(312)
plot(t,x,'g');
xlabel('Time (s)','FontSize',12)
ylabel('Amplitude','FontSize',12)
%ylim([-2,2])
%text(0,2,8,'(B)','Color','k','FontSize',16)
subplot(313)
plot(t,xx_LRC,'r');
xlabel('Time (s)','FontSize',12)
ylabel('Amplitude','FontSize',12)
ylim([-2,2])
%text(0,2,8,'(C)','Color','k','FontSize',16)
