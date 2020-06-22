function [ val ] = norm_fusedl1(X,beta)

val = sum(sum(abs(X)));
L=0;
for i=1:size(X,2)-1
        L = L + abs(X(:,i)-X(:,i+1));
end
val = val + beta* L;
end

