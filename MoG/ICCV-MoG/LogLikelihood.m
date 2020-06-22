function [L gamma] = LogLikelihood(Weight,mu,sigma,e)
k = length(Weight);
n = length(e);
L = 0;
for i = 1:n
    s = 0;
    for j = 1:k
        t1 = Weight(j)/sqrt(2*pi*sigma(j));
        t2 = max(0.0001,exp(-(e(i)-mu(j))^2/(2*sigma(j))));
        s = s + t1*t2;
        g(i,j) = Weight(j)*t2/sqrt(2*pi*sigma(j));
    end
    L = L+log(s);
end

for i = 1:n
    t = 0;
    for j = 1:k
        t = t + g(i,j);
    end
    gamma(i,:) = g(i,:)/t;
end