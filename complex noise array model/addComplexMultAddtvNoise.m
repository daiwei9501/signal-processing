function [pNoisy] = addComplexMultAddtvNoise(p,SNR)
% pNoisy = addComplexMultAddtvNoise(p,SNR) 
% Routine to add multiplicative and additive perturbations to a given
% signal p, with a given signal-to-noise ratio (SNR).
% 
%     p   = noise-free signal (e.g. simulated acoustic pressure)
%     SNR = Signal-to-Noise Ratio in [dB]
% 
% 
% References: 
%      
%       - Q. Leclère, "Acoustic imaging using under-determined inverse
%       approaches: Frequency limitations and optimal regularization,"
%       Journal of Sound and Vibration, vol. 321, no. 3-5, pp. 605 – 619,
%       2009.
%
%       - A. Pereira, J. Antoni, Q. Leclere, "Empirical Bayesian
%       regularization of the inverse acoustic problem", manuscript
%       submitted to Applied Acoustics.
%
%  
% A. Pereira, LVA-Insa de Lyon, 20/01/2013.
% antonio.pereira@insa-lyon.fr

numberOfmics = size(p,1);

alpha = randn(numberOfmics,1); %zero mean Gaussian variable with Variance = 1
delta = randn(numberOfmics,1); %zero mean Gaussian variable with Variance = 1

theta = 0+(2*pi-0).*rand(numberOfmics,1); %uniformly random distributed variable between [ 0 2pi ] 
phi = 0+(2*pi-0).*rand(numberOfmics,1);   %uniformly random distributed variable between [ 0 2pi ]

term1 = 10^(-SNR/20);
term2 = alpha.*exp(1i*theta).*p + delta.*exp(1i*phi)*sqrt(norm(p)^2/numberOfmics);

pNoisy = p + term1*term2;