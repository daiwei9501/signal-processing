clc
clear
paths ='C:\Users\70916\Desktop\signal\RPCA\RPCAKit-master\common';
addpath(paths);

rng(1);

rows = 100;
n_space = 5;
cluster_size = 50;

A = rand(rows, n_space) * rand(n_space, n_space);

permute_inds = reshape(repmat(1:n_space, cluster_size, 1), 1, n_space * cluster_size );
A = A(:, permute_inds);

corruption = 0.3;

N = randn(size(A)) * corruption;

X = A + N;

[A_est,N_est] = rpca_l1(X, 0.1);

rmpath(paths);