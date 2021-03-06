function [x,gg] = irgt_ij(c,g,R,I,J)
% [x,gg] = irgt(c,g,R,I)
% Inverse Gabor transform of ONE coefficient c in frequency bins I and time
% index J, based on a real analysis window g with time shift R.
% (gg is the resulting tappering near the signal edges)
%
%%%%%%%%%%%%%%%%%%%%%%%
% J. Antoni : Ao�t 2010
%%%%%%%%%%%%%%%%%%%%%%%

[Nw,K] = size(c);
Nw = 2*(Nw-1);
g = g(:);
L = (K-1)*R + Nw;


n = [0:Nw-1]';
c = exp(2i*pi*(I-1)*n/Nw)*c(I,J)/Nw;

x = zeros(L,1);
ind = [1:Nw] + (J-1)*R;
x(ind) = g.*c(mod(ind-1,Nw)+1);

if nargout > 1
    gg = zeros(L,1);
    ind = 1:Nw;
    for k = 1:K
        gg(ind) = gg(ind) + g.^2;
        ind = ind + R;
    end
end
