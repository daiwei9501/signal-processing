function [RL_matrix,RS_matrix,para_val]  = RPCA_FISTA(Li,N_iter,D_measured,SC,mu_final,mu,alp);

% This code can solve the problem
% minimize     \|L\|_{*} + \alpha\|S\|_{l1}  
% subject to   L + S = M
% this code can be used to detect the transient signal in the strong noise
% (Yu liang prepared in 09/02/2015)
    
% input parameters
% Li is the step size Li = 1.2 or 0.5
% N_iter is maximum iteration for each regularization step N_iter = 10 or 20
% D_measured is measured spectral matrix
% SC stopping criterion SC = 10^-2
% mu is initial regularization parameter mu = norm(D_measured,'fro')
% mu_final is final regularization parameter mu_final = 1e-16*mu
% alp is the regularization parameter between sparsity and low rank alp = 0.01

Lk = zeros(size(D_measured,1),size(D_measured,2)); % initial Low rank matrix
LYk = Lk;

Sk = zeros(size(D_measured,1),size(D_measured,2)); % initial sparse matrix
SYk = Sk;

mu_cal = mu;
tk = 1;
outer_iter = 1;
tic
while mu > mu_final
    mu = max(mu * 0.7, mu_final);
for i = 1:N_iter
    LGk = LYk + Li.* (D_measured - (LYk + SYk)); % gradient descent of low rank matrix
                  
    %LGk = (LGk + LGk')/2; % impose the hermitian property
    [U,S,V]= svd(LGk);
    Lkk = U*[diag(max(diag(abs(S)) - mu*Li,0)) zeros(size(D_measured,1),size(D_measured,2)-size(D_measured,1))]*V'; % low rank matrix estimation
    
    SGk = SYk + Li.* (D_measured - (LYk + SYk)); % gradient descent of sparse matrix
    Skk = max(abs(SGk)-alp*mu*Li,0); % sparse matrix estimation
    
    tk1 = 0.5 + 0.5*sqrt(1+4*tk^2); % update intermediate argument
    
    LYk = Lkk + ((tk - 1)/tk1)*(Lkk - Lk); % low rank matrix update
    Lk = Lkk;
    SYk = Skk + ((tk - 1)/tk1)*(Skk - Sk); % sparse matrix update
    Sk = Skk;
    
    tk = tk1; % update accelerate parameters
    
    err = norm(D_measured - (LYk + SYk),'fro')/norm(D_measured,'fro');
 
end
    total_iter = N_iter * outer_iter; 
    outer_iter = outer_iter + 1;
    if  err < SC; % stopping criterion
        break
    end
end

diagS = svd(Lk);
K = sum(diagS > mu_cal*Li);
disp(['rank of low rank matrix = ',num2str(K)])
disp(['number of iterations = ',num2str(total_iter)])

% output parameters
RL_matrix = Lk; % low rank matrix
RS_matrix = Sk; % sparse matrix
para_val.cost_time = toc; % cost time
para_val.mu = mu; % output parameter
para_val.err1 = err; % residual error
para_val.total_iter = total_iter; % total iteration steps
    