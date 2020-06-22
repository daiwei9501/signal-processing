function [x,gg] = icgt(c,g,R,I)
% [x,gg] = icgt(c,g,R,I)
% Inverse Gabor transform of coefficients c based on a complex analysis window
% g with time shift R.
% If nargin > 3, x is reconstructed from the vector of frequency bins I only.
% (gg is the resulting tappering near the signal edges)
%
%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni : Août 2010
%%%%%%%%%%%%%%%%%%%%%%%

[Nw,K] = size(c);
g = g(:);
L = (K-1)*R + Nw;

if nargin < 4
    c = ifft(c);
else
    C = zeros(Nw,K);
    n = [0:Nw-1]';
    for k = I
        C = C + exp(2i*pi*(k-1)*n/Nw)*c(k,:)/Nw;
    end
    c = C;
    clear C
end

x = zeros(L,1);
ind = 1:Nw;size(g),size(c),Nw
for k = 1:K
    x(ind) = x(ind) + g.*c(mod(ind-1,Nw)+1,k);
    ind = ind + R;
end
if nargout > 1
    gg = zeros(L,1);
    ind = 1:Nw;
    for k = 1:K
        gg(ind) = gg(ind) + abs(g).^2;
        ind = ind + R;
    end
end
