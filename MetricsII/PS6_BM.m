clear
clc

N = 1000; % Number of observations.

beta0 = 0.4; % Beta 0
beta1 = 0.5; % Beta 1
beta2 = 0.1; % Beta 2
beta3 = 0.6; % Beta 3
beta4 = 0.2; % Beta 4
beta5 = 0.3; % Beta 5

alpha = 1; % Alpha for Fuller.

sigma1 = 1; % Sigma 1
sigma2 = 2; % Sigma 2
rho = 0.8; % Rho

mu = [0 0]; % Mean of error terms.
sigma = [sigma1^2 rho*sigma1*sigma2; rho*sigma1*sigma2 sigma2^2]; % Var matrix

sim = 50; % Number of simus

bhat_ols = zeros(sim,3);  % OLS estimates
bhat_iv = zeros(sim,3);  % IV estimates
bhat_liml = zeros(sim,3); % LIML estimates
bhat_full = zeros(sim,3); % Fuller estimates

%% Estimation

z  = ((1:N)'./N).*normrnd(0,1,N,1); % Exog regressor
w  = ((1:N)'./N).*normrnd(0,1,N,1); % Instrument

for s = 1:sim
  
    e = mvnrnd(mu,sigma,N); % var cov matrix
    
    u = e(:,1); % Error terms for y                                                
    v = e(:,2); % Error terms for x
    
    x = beta3 + beta4*z + beta5*w + v; % Endog regressor
    y = beta0 + beta1*z + beta2*x + u;
    
    % OLS estimates.
    
    x_ols = [ones(N,1) z x];
    bols = inv(x_ols'*x_ols)*x_ols'*y; % OLS estimate of beta
    bhat_ols(s,:) = bols';

    % Results: OLS 

    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_ols(:,1)), std(bhat_ols(:,1)))
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_ols(:,2)), std(bhat_ols(:,2)))
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_ols(:,3)), std(bhat_ols(:,3)))
        
    % IV estimates.
    
    x_fs = [ones(N,1) z w];
    xhat = x_fs*inv(x_fs'*x_fs)*x_fs'*x;
    
    x_ss = [ones(N,1) z xhat];
    biv = inv(x_ss'*x_ss)*x_ss'*y; % IV estimate of beta
    bhat_iv(s,:) = biv';
    
    % Results: 2SLS
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_iv(:,1)), std(bhat_iv(:,1)))
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_iv(:,2)), std(bhat_iv(:,2)))
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_iv(:,3)), std(bhat_iv(:,3)))
   
    % LIML estimates.
    
    z_liml = [ones(N,1) z];
    w_liml = [ones(N,1) z w];
    
    Mz = eye(N) - z_liml*inv(z_liml'*z_liml)*z_liml';
    Mw = eye(N) - w_liml*inv(w_liml'*w_liml)*w_liml';
    
    Wz = [y x]'*Mz*[y x];
    Ww = [y x]'*Mw*[y x];
    kappa = min(eig((Ww^(-1/2))*Wz*(Ww^(-1/2))));
    bliml = inv(x_ols'*(eye(N)-(kappa*Mw))*x_ols)*(x_ols'*(eye(N)-(kappa*Mw))*y);
   
    bhat_liml(s,:) = bliml'; % LIML estimate of beta
    
    % Results: LIML
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_liml(:,1)), std(bhat_liml(:,1)))
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_liml(:,2)), std(bhat_liml(:,2)))
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_liml(:,3)), std(bhat_liml(:,3)))

    % Fuller LIML estimates
    
    z_liml = [ones(N,1) z];
    w_liml = [ones(N,1) z w];
    
    Mz = eye(N) - z_liml*inv(z_liml'*z_liml)*z_liml';
    Mw = eye(N) - w_liml*inv(w_liml'*w_liml)*w_liml';
    
    Wz = [y x]'*Mz*[y x];
    Ww = [y x]'*Mw*[y x];
    kappa = min(eig((Ww^(-1/2))*Wz*(Ww^(-1/2)))) - alpha/(N-3);
    bfull = inv(x_ols'*(eye(N)-(kappa*Mw))*x_ols)*(x_ols'*(eye(N)-(kappa*Mw))*y);
    
    bhat_full(s,:) = bfull'; % LIML estimate of beta with Fuller
    
    % Results: Fuller
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_full(:,1)), std(bhat_full(:,1)))
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_full(:,2)), std(bhat_full(:,2)))
    
    fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_full(:,3)), std(bhat_full(:,3)))

end


