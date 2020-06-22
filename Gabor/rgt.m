function c = rgt(x,g,R)
% c = rgt(x,g,R)
% Gabor transform of signal x based on a real analysis window g with time shift R.
% Coefficients are arranged in a matrix c with column index corresponding
% to (positive) frequencies and raw index to time shifts.
%
%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni : Août 2010
%%%%%%%%%%%%%%%%%%%%%%%

x = x(:);
g = g(:);
L = length(x);
Nw = length(g);
if rem(Nw,2 > 0)
    error('the length of window g must be even !');
end
K = floor((L-Nw)/R+1);

c = zeros(Nw/2+1,K);
ind = 1;
wi = 2i*pi*(0:Nw-1)'/Nw;
for k = 1:K
    temp = fft(x(ind:ind+Nw-1).*g).*exp(-wi*(k-1)*R);
    c(:,k) = temp(1:Nw/2+1);
    ind = ind + R;
end