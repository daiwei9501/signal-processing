clc
clear
close all

% This MATLAB code is to simulate the conventional beamforming and DAMAS,
% with a version that can work in a strong noise

%% 
m = 121;
n = 1000;
r = 3;

%% physical parameters
rhoAir = 1.2;  % air density
c = 340;  % speed of sound
Nsnap = 1000;  % number of snapshots
fre_index = 6000;  % Hz working frequency
Ns = 3; % # of source

% source modeling (equivelent sources model)
d_calcul = .05;  
xMin = -.5;
xMax = .5;
yMin = -.5;
yMax = .5;

xD = xMin:d_calcul:xMax;
yD = yMin:d_calcul:yMax;
Nx = length(xD);
Ny = length(yD);

[Xm,Ym] = meshgrid(xD,yD);
Zm = zeros(size(Xm(:)));
eqSources_Coord = [Xm(:), Ym(:) Zm];

% real sources (in the grids of the equivelent sources model)
realSrc_Coord = [-0.3 0.25 0;0.25 0.25 0;0.35 -0.25 0];
source_pos = [16 5;16 16;6 18];

%Src_Src = [-0.3 0.25;0.25 0.25;0.35,-0.25];

% microphone coordinates modeling
dxy = .2*ones(1,2); % space between microphones x and y (m)
xymin = -1*ones(1,2); % regular array
xymax = 1*ones(1,2); % regular array
zm = 1;

[Xm,Ym,Zm] = Antenne_rect(dxy,xymin,xymax,zm);% Renvoie les positions Xm,Ym,Zm d'une antenne rectangulaire
Xm = Xm - mean(Xm);
Ym = Ym - mean(Ym);
array_Coord = [Xm,Ym,Zm] ;

%show the microphone and sources topology
figure;
plot3(eqSources_Coord(:,1),eqSources_Coord(:,2),eqSources_Coord(:,3),'.b');
hold on;
plot3(array_Coord(:,1),array_Coord(:,2),array_Coord(:,3),'Or','LineWidth',1.5);
plot3(realSrc_Coord(1,1),realSrc_Coord(1,2),realSrc_Coord(1,3), '.r','markersize',30);
plot3(realSrc_Coord(2,1),realSrc_Coord(2,2),realSrc_Coord(2,3), '.r','markersize',30);
plot3(realSrc_Coord(3,1),realSrc_Coord(3,2),realSrc_Coord(3,3), '.r','markersize',30);
%title('Simulation','Fontname','Times New Roman','Fontsize',24);

%%
% the distance matrix between the microphones and sources
dist_RealSource_Micro = getDistanceSrcMic(realSrc_Coord,array_Coord); % 3*121
dist_eqSource_Micro = getDistanceSrcMic(eqSources_Coord,array_Coord); % 441*121
%%
% algorithm paramerters initiation
numberOfMicro = size(array_Coord,1);
numberOfeqSources = size(eqSources_Coord,1);
pdirect = zeros(numberOfMicro,1);
p = pdirect;

% propogation function
H_realbeam = transferMatrixbeamforming(fre_index,c,dist_RealSource_Micro); % 3*121
H_discretebeam = transferMatrixbeamforming(fre_index,c,dist_eqSource_Micro); % 441*121
src_Strength = (randn(Ns,Nsnap) + 1i*randn(Ns,Nsnap))/sqrt(2); % 3*1000
%%
Q_ref = 4*10^(-10);
power_source = var(src_Strength.');
power_dB = 10*log10(power_source/Q_ref);

pdirect = H_realbeam * src_Strength; %
M = size(array_Coord,1);

% SNR1 = -20;; % signal to nosie ratio
% SNR2 = -20; % signal to nosie ratio
% n = Gauss_noise( SNR1,SNR2,power_signal,numberOfMicro,Nsnap );
% p = pdirect + n; % 30*100

SNR = -20;
power_signal = mean((std(pdirect',0,1)).^2);
power_noise = power_signal/(10^(0.1*SNR));

%%
%加噪声
Ind = randperm(m*n);
p1 = floor(m*n*0.2);
p2 = floor(m*n*0.2);
p3 = floor(m*n*0.2);
p = pdirect(:).';

%p = p + (rand(1,m*n)* 1 - 0.5 + (1i*(rand(1,m*n)* 1 - 0.5) + randn(1,m*n)* 0.02 + 1i*randn(1,m*n)*0.02 + randn(1,m*n)*0.01 + 1i*randn(1,m*n)*0.01)*sqrt(power_noise));

p(Ind(p1+1:p1+p2)) = p(Ind(p1+1:p1+p2)) + (rand(1,p2)* 10 - 5 + 1i*(rand(1,p2)* 10 - 5))*sqrt(power_noise);

p(Ind(p1+p2+1:p1+p2+p3)) = p(Ind(p1+p2+1:p1+p2+p3)) + (randn(1,p3)* 0.2 + 1i*randn(1,p3))*sqrt(power_noise);
p(Ind(p1+p2+p3+1:end)) = p(Ind(p1+p2+p3+1:end)) + (randn(1,m*n-p1-p2-p3)*0.01 + 1i*randn(1,m*n-p1-p2-p3)*0.01)*sqrt(power_noise); 
p = reshape(p, m, n);

%p = pdirect+(randn(M,Nsnap)+1i*randn(M,Nsnap))*sqrt(power_noise)/sqrt(2);

%p_noise1 = (0.5 + 1i*(rand(m,1)* 1 - 0.5))*(sqrt(2))*sqrt(power_noise);
%p_noise2 = (1 + 1i*randn(m,1))/sqrt(2)*sqrt(power_noise);
%p = pdirect;
%p(:,3) = p(:,3) + p_noise1;
%p(:,10) = p(:,10) + p_noise2;
%%
% spectral matrix of the microphone measurements
Spp = p*p'/Nsnap; % 30*30
Spp = (Spp + Spp')./2;
%
Sss = pdirect*pdirect'/Nsnap; % 30*30
Sss = (Sss + Sss')./2;

% 对p矩阵进行去噪声操作

aa = median(real(Spp(:)));
aa = sqrt(aa/r);
bb = median(imag(Spp(:)));
bb = sqrt(bb/r);
U0 = rand(m,r)*aa*2-aa + 1i*(rand(m,r)*bb*2-bb);
V0 = rand(n,r)*aa*2-aa + 1i*(rand(n,r)*bb*2-bb);
        
%%%%%%%%%%%%%%%%%%GMM method %%%%%%%%%%%%%%%%%%%%%%%%
param.maxiter = 100;
param.OriX = pdirect;
param.InU = U0;
param.InV = V0;
param.k = 3;
param.display = 0;
param.NumIter = 100;
param.tol = 1.0e-5;%param.tol = 1.0e-5;
param.method = 4;
W = ones(size(p));
[label, model,TW, A,B,llh] =  MLGMDN(W,p,r,param);
%%
% Perform MoG-RPCA
lambda = 0.1; % 
[LR, SR] = proximal_gradient_rpca(Spp, lambda);
Spp_RPCA = LR;
%%
p_denoise = A*B';
Spp_Mog = p_denoise*p_denoise'/Nsnap;
Spp_Mog = (Spp_Mog + Spp_Mog')./2;
%%
figure
%subplot(311);
imagesc(abs(Sss));
colorbar()
title('without noise')
%%
figure,
%subplot(312);
imagesc(abs(Spp));
colorbar()
title('with strong noise')
%%
%subplot(313);
figure,
imagesc(abs(Spp_Mog));
colorbar()
title('de-noised version');
%%
figure,
imagesc(abs(Spp_RPCA));
colorbar()
title('RPCA denoised')
%%
% figure;
% plot(realSrc_Coord(1,1),realSrc_Coord(1,2), '.r','markersize',30);
% hold on
% plot(realSrc_Coord(2,1),realSrc_Coord(2,2), '.g','markersize',30);
% hold on
% plot(realSrc_Coord(3,1),realSrc_Coord(3,2), '.m','markersize',30);
% axis([-0.6 0.6 -0.6 0.6]);
% axis ij
% title('real sources')

% measurement error with strong noise
percent_org = sum(abs(Spp-Sss).^2)/sum(abs(Sss).^2) 

% measurement error with denosing method
percent_denoy = sum(abs(Spp_Mog-Sss).^2)/sum(abs(Sss).^2)

% RPCA
percent_RPCA = sum(abs(Spp_RPCA-Sss).^2)/sum(abs(Sss).^2)
%% the following codes can be ingored 

%% frequency beamforming with strong noise
q_sources = size(numberOfeqSources,1);
steering_matrix = zeros(size(H_discretebeam));
for index_sources = 1:numberOfeqSources
     % direct steering vector is vector length of M
     norm_vector = norm(H_discretebeam(:,index_sources),2);
     % steering_vector = H_discretebeam(:,index_sources); % H_discretebeam: 30 * 483 
     steering_vector = H_discretebeam(:,index_sources)/(norm_vector);
     steering_matrix(:,index_sources) = steering_vector;
     q_beamforming(index_sources) =  steering_vector'*Spp*steering_vector;
end

q_recon_beamforming = reshape(abs(q_beamforming),length(yD),length(xD));% sum for each row

figure;
% Q_ref = max(max(q_recon_beamforming));% 50*10^-9 = (50 nm^3/s) reference value for dB unit
imagesc(xD,yD,10*log10(q_recon_beamforming./Q_ref));% plot in dB unit,hi is the 
% axis xy
dynamique = 10; % set the dynamic range
max_Value = max(max(10*log10(q_recon_beamforming./Q_ref)));
set(gca, 'clim', [-dynamique 0]+ max_Value);
title('direct beamforming')


% DAMAS method to deconvolve the blurry  beamforming results for high spatial resolution

C = (abs(steering_matrix'*H_discretebeam)).^2;

options.TolX = 0.05;
q_damas = lsqnonneg(C,abs(q_beamforming.'), options);

q_recon_damas = reshape(abs(q_damas),length(yD),length(xD));% sum for each row

figure;
% Q_ref = max(max(q_recon_damas));% 50*10^-9 = (50 nm^3/s) reference value for dB unit
imagesc(xD,yD,10*log10(q_recon_damas./Q_ref));% plot in dB unit,hi is the 
colorbar()
% axis xy
dynamique = 10; % set the dynamic range
max_Value = max(max(10*log10(q_recon_damas./Q_ref)));
set(gca, 'clim', [-dynamique 0]+ max_Value);
%title('direct DAMAS')

q_recon_damas_outDN = 10*log10(q_recon_damas./Q_ref);
Q_rec1_outDN = q_recon_damas_outDN(source_pos(1,1),source_pos(1,2));
Q_rec2_outDN = q_recon_damas_outDN(source_pos(2,1),source_pos(2,2));
Q_rec3_outDN = q_recon_damas_outDN(source_pos(3,1),source_pos(3,2));

error1 = [Q_rec1_outDN Q_rec2_outDN Q_rec3_outDN] - power_dB

%% frequency beamforming with denoising method
q_sources = size(numberOfeqSources,1);
steering_matrix = zeros(size(H_discretebeam));
for index_sources = 1:numberOfeqSources
     % direct steering vector is vector length of M
     norm_vector = norm(H_discretebeam(:,index_sources),2);
     % steering_vector = H_discretebeam(:,index_sources); % H_discretebeam: 30 * 483 
     steering_vector = H_discretebeam(:,index_sources)/(norm_vector);
     steering_matrix(:,index_sources) = steering_vector;
     q_beamforming(index_sources) =  steering_vector'*Spp_Mog*steering_vector;
end

q_recon_beamforming = reshape(abs(q_beamforming),length(yD),length(xD));% sum for each row

figure;
% Q_ref = max(max(q_recon_beamforming));% 50*10^-9 = (50 nm^3/s) reference value for dB unit
imagesc(xD,yD,10*log10(q_recon_beamforming./Q_ref));% plot in dB unit,hi is the 

% axis xy
dynamique = 10; % set the dynamic range
max_Value = max(max(10*log10(q_recon_beamforming./Q_ref)));
set(gca, 'clim', [-dynamique 0]+ max_Value);
title('de-noised beamforming')
colorbar()
% DAMAS method to deconvolve the blurry  beamforming results for high spatial resolution

C = (abs(steering_matrix'*H_discretebeam)).^2;

options.TolX = 0.05;
q_damas = lsqnonneg(C,abs(q_beamforming.'), options);

q_recon_damas = reshape(abs(q_damas),length(yD),length(xD));% sum for each row

figure;
% Q_ref = max(max(q_recon_damas));% 50*10^-9 = (50 nm^3/s) reference value for dB unit
imagesc(xD,yD,10*log10(q_recon_damas./Q_ref));% plot in dB unit,hi is the 
% axis xy
dynamique = 10; % set the dynamic range
max_Value = max(max(10*log10(q_recon_damas./Q_ref)));
set(gca, 'clim', [-dynamique 0]+ max_Value);
%title('de-noised DAMAS')
colorbar()
q_recon_damas_DN = 10*log10(q_recon_damas./Q_ref);
Q_rec1_DN = q_recon_damas_DN(source_pos(1,1),source_pos(1,2));
Q_rec2_DN = q_recon_damas_DN(source_pos(2,1),source_pos(2,2));
Q_rec3_DN = q_recon_damas_DN(source_pos(3,1),source_pos(3,2));

error2 = [Q_rec1_DN Q_rec2_DN Q_rec3_DN] - power_dB


%%
%% frequency beamforming with denoising method
q_sources = size(numberOfeqSources,1);
steering_matrix = zeros(size(H_discretebeam));
for index_sources = 1:numberOfeqSources
     % direct steering vector is vector length of M
     norm_vector = norm(H_discretebeam(:,index_sources),2);
     % steering_vector = H_discretebeam(:,index_sources); % H_discretebeam: 30 * 483 
     steering_vector = H_discretebeam(:,index_sources)/(norm_vector);
     steering_matrix(:,index_sources) = steering_vector;
     q_beamforming(index_sources) =  steering_vector'*Spp_RPCA*steering_vector;
end

q_recon_beamforming = reshape(abs(q_beamforming),length(yD),length(xD));% sum for each row

figure;
% Q_ref = max(max(q_recon_beamforming));% 50*10^-9 = (50 nm^3/s) reference value for dB unit
imagesc(xD,yD,10*log10(q_recon_beamforming./Q_ref));% plot in dB unit,hi is the 

% axis xy
dynamique = 10; % set the dynamic range
max_Value = max(max(10*log10(q_recon_beamforming./Q_ref)));
set(gca, 'clim', [-dynamique 0]+ max_Value);
title('de-noised RPCA beamforming')
colorbar()
% DAMAS method to deconvolve the blurry  beamforming results for high spatial resolution

C = (abs(steering_matrix'*H_discretebeam)).^2;

options.TolX = 0.05;
q_damas = lsqnonneg(C,abs(q_beamforming.'), options);

q_recon_damas = reshape(abs(q_damas),length(yD),length(xD));% sum for each row

figure;
% Q_ref = max(max(q_recon_damas));% 50*10^-9 = (50 nm^3/s) reference value for dB unit
imagesc(xD,yD,10*log10(q_recon_damas./Q_ref));% plot in dB unit,hi is the 
% axis xy
dynamique = 10; % set the dynamic range
max_Value = max(max(10*log10(q_recon_damas./Q_ref)));
set(gca, 'clim', [-dynamique 0]+ max_Value);
%title('de-noised RPCA DAMAS')
colorbar()
q_recon_damas_DN = 10*log10(q_recon_damas./Q_ref);
Q_rec1_DN = q_recon_damas_DN(source_pos(1,1),source_pos(1,2));
Q_rec2_DN = q_recon_damas_DN(source_pos(2,1),source_pos(2,2));
Q_rec3_DN = q_recon_damas_DN(source_pos(3,1),source_pos(3,2));

error2 = [Q_rec1_DN Q_rec2_DN Q_rec3_DN] - power_dB

















