function [distanceSourceMicro] = getDistanceSrcMic(source_Coordinates,array_Coordinates)
%
%
%
%
numberOfmicro = size(array_Coordinates,1);
distanceSourceMicro = zeros(numberOfmicro,size(source_Coordinates,1));
for ii = 1 : numberOfmicro
    distanceSourceMicro(ii,:) = sqrt( (source_Coordinates(:,1)-array_Coordinates(ii,1)).^2 +...
        (source_Coordinates(:,2)-array_Coordinates(ii,2)).^2 + ...
        (source_Coordinates(:,3)-array_Coordinates(ii,3)).^2);
end