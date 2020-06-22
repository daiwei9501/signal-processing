function c = cgt(x,g,R)
% c = cgt(x,g,R)
% Gabor transform of signal x based on analysis window g with time shift R.
% Coefficients are arranged in a matrix c with raw indices corresponding
% to (positive) frequencies and column indices to time shifts.
%
%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni : Août 2010
%%%%%%%%%%%%%%%%%%%%%%%

x = x(:);
g = g(:);%hanning window
L = length(x);%the length of the signal
Nw = length(g);
if rem(Nw,2 > 0)
    error('the length of window g must be even !');
end
K = floor((L-Nw)/R+1);%the number of move times

c = zeros(Nw/2+1,K);
% c = zeros(Nw,K);
ind = 1;
wi = 2i*pi*(0:Nw-1)'/Nw;
for k = 1:K
    temp = fft(x(ind:ind+Nw-1).*conj(g)).*exp(-wi*(k-1)*R);
    c(:,k) = temp(1:Nw/2+1);% number of the frequency (sysmetric)
%     c(:,k) = temp(1:Nw);
    ind = ind + R;
end
