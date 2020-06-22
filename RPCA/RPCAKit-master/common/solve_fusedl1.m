function [ E ] = solve_fusedl1( W, lambda,beta)
%MYSOLVE_L1L2 Summary of this function goes here
%   Detailed explanation goes here

[m, n] = size(W);
E = W;
for i = 1 : n 
    norm_col = norm_fusedl1(W(:,i),beta);
    if (norm_col > lambda)
        E(:,i) = (norm_col - lambda) * W(:,i) / norm_col;
    else
        E(:,i) = zeros(m, 1);
    end 
end
end
