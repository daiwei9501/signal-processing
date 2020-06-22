function c = cgt_freq(X,g,R)
% c = cgt_freq(X,g,R)
% Gabor transform of signal x based on analysis window g with time shift R.
% Coefficients are arranged in a matrix c with raw indices corresponding
% to (positive) frequencies and column indices to time shifts.
% Frequency domain computation
%
%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni : Nov. 2010
%%%%%%%%%%%%%%%%%%%%%%%

X = X(:);
g = g(:);
Nfft = length(X);
Nw = length(g);
if rem(Nw,2 > 0)
    error('the length of window g must be even !');
end
K = floor((Nfft-Nw)/R+1);

f = (0:Nfft-1)'/Nfft;
c = zeros(Nw/2+1,K);
for i = 0:Nw/2
    Gi = conj(fft(g.*exp(2i*pi*i*(0:Nw-1)'/Nw),Nfft));
    for k = 0:K-1
        c(i+1,k+1) = mean(Gi.*X.*exp(2i*pi*k*R*f));
    end
    c(i+1,:) = c(i+1,:).*exp(-2i*pi*(0:K-1)*R*i/Nw);
end