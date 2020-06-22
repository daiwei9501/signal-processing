% Teste la fonction cgt

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
close all
clear all

L = 1e4;
Nw = 2^10;
window = 1;
switch window
    case 1
        R = Nw/4;
        g = hanning(Nw,'periodic')/sqrt(1.5);
    case 2
        R = Nw/2;
        g = sqrt(hanning(Nw,'periodic'));
end
% g = hilbert(g);

K = floor((L-Nw)/R+1);
G = zeros(L,1);
ind = 1;
for k = 1:K
    G(ind:ind+Nw-1) = G(ind:ind+Nw-1) + g.^2;
    ind = ind + R;
end
figure
plot(G)

t = 0:L-1;     
% x = chirp(t,0,L-1,.5);
x = chirp(t,0,L-1,.5,'q');
% x = ones(L,1);
x = x(:);
disp('cgt:'),tic,
c = cgt(x,g,R);
toc
figure,
imagesc(abs(c)),axis xy

disp('igt:'),tic,
[xr,w] = igt(c,g,R);
xr = real(xr);
w = real(w);
% [xr,w] = igt(c,g,R,10:50);
% xr = 2*real(xr);
toc
figure
subplot(211),plot(x(1:length(xr))),hold on,plot(x(1:length(xr)).*w/max(w),'g'),plot(real(xr),':r')
% max(imag(xr))/max(real(xr))
% subplot(212),plot(abs(fft(x))),hold on,plot(abs(fft(real(xr))),':r')
subplot(212),plot(real(xr)-x(1:length(xr)).*w/max(w)),
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('dgt:'),tic,c=dgt(x,g,R,Nw);toc
disp('idgt:'),tic,xr=idgt(c,g,R,Nw);toc
