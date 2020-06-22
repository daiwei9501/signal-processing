function [H] = transferMatrixbeamforming(freq,c,distSourceMicro)


H = exp(-1i*2*pi*freq*(distSourceMicro./c))./(distSourceMicro);