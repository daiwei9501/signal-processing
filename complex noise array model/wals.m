function [OutU,OutV] = wals(X,W,InU,InV,MaxIter,tol)
% matrix factorization with weighted alternating least square method
% Input:
%   X: d x n input data matrix
%   W: d x n weight matrix
%   InU: initialization of one factor
%   InV: initialization of the other factor
%   MaxIter: maximum iteration
%   tol: the threshold for stop
% Output:
%   OutU: one predicted factor
%   OutV: the other predicted factor
% author: licong(mathelc@qq.com)20190819 

%% regularization parameter
lambda_u = 0.0005;
lambda_v = 0.0005;
MaxIter = 100;

[r] = size(InU,2);
[d,n] = size(X);
OutU = InU;
OutV = InV;
for i = 1:MaxIter
    for jj = 1:n
        temp1 = OutU.';
        temp2 = reshape(kron(W(:,jj),ones(r,1)), [1, d*r]).*reshape(temp1(:),[1,d*r]);
        OutV(jj,:) = reshape((W(:,jj).*W(:,jj).*X(:,jj)+lambda_u),[1,d])*OutU/dot(temp2,temp2);
    end
    for ii = 1:d
        temp3 = OutV.';
        temp4 = kron(W(ii,:),ones(1,r)).*reshape(temp3(:),[1,n*r]);
        OutU(ii,:) = (W(ii,:).*W(ii,:).*X(ii,:)+lambda_v)*OutV/dot(temp4,temp4);
    end    

    if norm(X - OutU*OutV') < tol
        disp('hel');
        break;
    else
        InU = OutU;
    end
end
end


