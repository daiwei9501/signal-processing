function [ n ] = Gauss_noise( SNR1,SNR2,power_signal,numberOfMicro,Nsnap )
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%
SNR=SNR1+(SNR2-SNR1)*rand(numberOfMicro,1);
sigB=sqrt(power_signal./(10.^(0.1.*SNR)));
for i=1:numberOfMicro
    noise(i,:)=sigB(i).*(randn(1,Nsnap) + 1i*randn(1,Nsnap))/sqrt(2);
end
n=noise;
end

