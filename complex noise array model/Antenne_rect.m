function [Xm,Ym,Zm] = Antenne_rect(dxy,xymin,xymax,zm)
% [Xm,Ym,Zm] = Antenne_rect(dxy,xymin,xymax,zm)
% Renvoie les positions Xm,Ym,Zm d'une antenne rectangulaire

% Position des microphones
x = xymin(1):dxy(1):xymax(1);
y = xymin(2):dxy(2):xymax(2);
Mx = length(x);
My = length(y);
Xm = repmat(x,My,1);
Xm = Xm(:);
Ym = repmat(y',Mx,1);
% [Xm,Ym] = meshgrid(x,y);
%  Xm = Xm(:);
%  Ym = Ym(:);
Zm = zm*ones(Mx*My,1);