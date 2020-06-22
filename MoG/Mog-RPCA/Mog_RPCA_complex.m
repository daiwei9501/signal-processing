clc
clear
close all
addpath(genpath('C:\Yuliang\MATLAB experment\insa project\Rerange the program\signal dection\surveillence 7 paper\apg'));
addpath(genpath('C:\Yuliang\MATLAB experment\insa project\Rerange the program\signal dection\surveillence 7 paper\inexact_alm_rpca'));
addpath(genpath('C:\Yuliang\MATLAB experment\insa project\Rerange the program\signal dection\surveillence 7 paper\dual'));
addpath(genpath('C:\Yuliang\MATLAB experment\insa project\Rerange the program\signal dection\surveillence 7 paper\apg_partial'));
addpath(genpath('C:\Yuliang\MATLAB experment\insa project\Rerange the program\signal dection\surveillence 7 paper\svt'));
addpath(genpath('C:\Yuliang\MATLAB experment\insa project\Rerange the program\signal dection\surveillence 7 paper\exact_alm_rpca'));
addpath('E:\reseach\MATLAB Codes and data ensemble\research codes\matlab codes signal detection\Gabor')
addpath(genpath('C:\Yuliang\Research\MATLAB Codes and data ensemble\research codes\Spetral_matrix_denoising\apg\'));
%addpath(genpath('E:\reseach\MATLAB Codes and data ensemble\research codes\matlab codes signal detection\Gabor\'));
addpath('C:\Users\70916\Desktop\signal\低秩分解\Gabor')
addpath('C:\Users\Administrator\Documents\MATLAB\mog')
addpath('C:\Users\Administrator\Documents\MATLAB\matlab codes signal detection\Gabor')

%{
     This program is to detect the transient signal using the RPCA
     programmer : Jerome antoni and Yu Liang 08/02/2012
                  xiaobawang1984@gmail.com
     Supervisors : Jerome Antoni and Quentin Leclere
%}

%% simulate the signal
L = 1e4;        % signal length
T = 300;        % mean period of impulses

impulses = zeros(L,1); % creat signal
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
noise=zeros(L,1);
SNR =  -5; %SIGNAL-NOISE RATIO- dB                
% signal = x + NSR*noise*std(x);
%power_signal =  std(x); % standard derviation -std(x,1) normalizes by N and produces the square root of the second
   %  moment of the sample about its mean
%power_noise=power_signal/( 10^(0.05*SNR));
power_signal = std(x)^2;
power_noise = (power_signal/(10^(0.1*SNR)))^0.5
%% i.i.d. Gaussian Noise
if 1 == 1
    noise = randn(L,1)*power_noise;
end
%% Non-i.i.d. Gaussian Noise- zero-mean Gaussian noise 
%and The signal noise ratio (SNR) value of each band is generated from uniform distribution
%with value in the range of [-5, 5]dB
if 1 == 0
    SNR =  -5 + rand(1,L)*10; 
    sigma_signal = power_signal;
    sigma_noi = sigma_signal./(10.^(SNR./20));
    for i=1:L
       noise(i,1) = noise(i,1) + randn *sqrt(sigma_noi(i));
    end
end

%% Rayleigh noise ,expectation is  sqrt(pi/2)*sigma
if 1 == 0
    sigma = sqrt(2/(4-pi));  %std of a rayleigh noise is 1 ,when sigma = sqrt(2/(4-pi))
    noise = raylrnd(sigma,L,1) - sqrt(pi / 2) * sigma;
    noise = noise * power_noise;
end

 %% Gaussian + Stripe Noise
%     SNR = 10 + rand(1,B)*10;
%     SNR1 = 10.^(SNR./10);
%     sigma_noi = sigma_signal./SNR1;
%     for i=1:B
%         Noi_H(:,:,i) = Ori_H(:,:,i) + randn(M,N)*sqrt(sigma_noi(i));
%     end
%     band = ceil((B-20)*rand(1,40)+10);              % bands chose to add impluse noise
%     rate = 0.5+rand(1,length(band))*0.2;            % rate of impluse noise in these bands

%% colored nois-Specify the color of the noise as one of 
%'pink' | 'white' | 'brown' | 'blue' | 'purple' | 'custom'
%https://en.wikipedia.org/wiki/Noise_(electronics)
% 
% hcn=dsp.ColoredNoise('Color','pink','SamplesPerFrame',L);
% noise=hcn();

%% gaussians overlay
% SNR = -5;
% SNR2 = 1; 
% power_signal =  std(x); % standard derviation -std(x,1) normalizes by N and produces the square root of the second
% power_noise=power_signal/( 10^(0.1*SNR));
% power_noise2=power_signal/( 10^(0.1*SNR2));
% noise=0.3*randn(L,1)*power_noise+0.5*randn(L,1)*power_noise+0.2*randn(L,1)*power_noise;
% noise=0.3*randn(L,1)*power_noise+0.5*randn(L,1)*power_noise+randn(L,1)*power_noise2;
%% mutimodel gaussian
if 1 == 0
    A = rand(L,1)*3;
    k=1;
    noise = zeros(L,1);
    for i =1:L 
        noise(i,1) = (randn(1)+k)*(A(i,1)<=1) +  (randn(1)-k)*(A(i,1)<=2 && A(i,1)>2 )+(randn(1))*(A(i,1)<=3 && A(i,1)>2) ;
    end
    noise =noise *power_noise;
end
%% gamma noise/Erlang noise ,also Exponential Distribution Noise when b = 1 
if 1 == 0
    a = 15;
    b = a^2;
    n_Erlang = zeros(L,1); 

    for j=1:b
        n_Erlang = n_Erlang + (-1/a)*log(1 - rand(L,1));
    end
    n_Erlang = n_Erlang - b/a;
    noise = n_Erlang * power_noise;
end

%% Uniform noise
if 1 == 0
    a = -sqrt(12)/2;
    b = sqrt(12)/2;
    n_Uniform = a + (b-a)*rand(L,1);
    noise = n_Uniform * power_noise ;
end

 
%% 
   
signal  = x  + noise;

figure;
plot(signal,'b');
hold on
plot(x,'g');

% Gabor analysis without zero-padding
Nw = 2^6;
R = Nw*2/8;
% R = Nw*1/8;
g = hanning(Nw,'periodic')/sqrt(1.5);
c = cgt(signal,g,R);

figure;
imagesc(abs(c))

%% construct the rank-1 matrix
M = real(c);
M2 = imag(c);


% M = abs(c).^2;
% figure;
% imagesc(M)
% title('powerspectrum of a time domain signal + noise','FontSize',14);
% xlabel('Time(seconds)','FontSize',14);
% ylabel('Frequency(Hz)','FontSize',14);

%% 谱上加噪声 噪声实部虚部分别符合高斯分布
if 1 ==0

    c1 = ( randn(size(c))  + randn(size(c))*i )* power_noise;
    M1 = abs(c1).^2 ;
    figure;
    imagesc(M1)
    title('Noise on PowerSpectrum','FontSize',14);
    xlabel('Time(seconds)','FontSize',14);
    ylabel('Frequency(Hz)','FontSize',14);
    M = M + M1;
    figure;
    imagesc(M)
    title('Signal and Noise on PowerSpectrum','FontSize',14);
    [signal,gg] = icrgt(M,g,R);
end




%% 
% Perform MoG-RPCA
param.mog_k = 3;                         %k个高斯分部混合
param.lr_init = 'SVD';
param.maxiter = 100;
param.initial_rank = 20;                 %low rank的初始rank 数，不要设太小
param.tol = 1e-3;
    
lr_prior.a0 = 1e-6;
lr_prior.b0 = 1e-6;

mog_prior.mu0 = 0;
mog_prior.c0 = 1e-3;
mog_prior.d0 = 1e-3;
mog_prior.alpha0 = 1e-3;
mog_prior.beta0 = 1e-3;

[LR, SR, r] = mog_rpca(M, param, lr_prior, mog_prior);
SR
[LR2, SR2, r2] = mog_rpca(M2, param, lr_prior, mog_prior);
SR2
%% plot the reconstructed sparse and low rank matrix
LR = LR.U*LR.V'; % low rank matrix
LR2 = LR2.U*LR2.V';

for i=1:param.mog_k
    figure;
    imagesc(10*log10(abs(SR.R(:,:,i))));
    title(strcat('Real Part Sparse decomposition ',num2str(i)),'FontSize',14);
    xlabel('Time(seconds)','FontSize',14);
    ylabel('Frequency(Hz)','FontSize',14);
end

for i=1:param.mog_k
    figure;
    imagesc(10*log10(abs(SR2.R(:,:,i))));
    title(strcat('Imaginary Part Sparse decomposition',num2str(i)),'FontSize',14);
    xlabel('Time(seconds)','FontSize',14);
    ylabel('Frequency(Hz)','FontSize',14);
end

%% plot the reconstructed sparse and low rank matrix
figure
imagesc(10*log10(abs(LR)))
title('Real Part Low rank decomposition L','FontSize',14)

figure
imagesc(10*log10(abs(LR2)))
title('Imaginary Part Low rank decomposition L','FontSize',14)
%clim = get(gca, 'clim');
cmap=colormap('hot');
colormap(cmap(end:-1:1, :));
% set(gca, 'clim',[0 max(max(max(LR)),max(max(SR)))]);
% set(gca, 'clim', max(max(max(10*log10(abs(LR)))),max(max(10*log10(abs(SR))))) + [-30 0])
xlabel('Time(seconds)','FontSize',14)
ylabel('Frequency(Hz)','FontSize',14)



%% reconstruct the complex low rank and sparse matrix
SR_R= SR.R;
SR2_R = SR2.R;
size_SR_R=size(SR_R);
size_SR2_R=size(SR2_R);
SR = 0;

for i=1:size_SR_R(3)
SR = SR +SR_R(:,:,i) ; % 1, 2, 3
end

SR2 = 0;
for i=1:size_SR2_R(3)
SR2 = SR2 + SR2_R(:,:,i) ; % 1, 2, 3
end

LRC = LR + LR2*j;
gain = 2/6;
mask = double(abs(LRC)> (gain*max(max(abs(LRC))))); 
LRC = LRC.*mask;

figure
imagesc(abs(LRC));
cmap=colormap('hot');
colormap(cmap(end:-1:1, :));
% set(gca, 'clim',[0 max(max(max(abs(LRC))),max(max(abs(SRC))))]);
title('Sparse decomposition (after the mask)','FontSize',14)
xlabel('Time(seconds)','FontSize',14)
ylabel('Frequency(Hz)','FontSize',14)

[xx_LRC,gg] = icrgt(LRC,g,R);

figure
plot(signal)
hold on
plot(x,'g'),
hold on
plot(real(xx_LRC),'r')
title('Nw = 256  lambda = 0.4','FontSize',14)
legend('Measurements','Transient Signal','Separated Transient Signal')
xlabel('Sample Number','FontSize',14)
ylabel('Amplitude','FontSize',14)


