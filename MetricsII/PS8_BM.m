%% NOTE: 
% This code includes both question 1 and 2. Each question has a clear 
% clc command, which would get rid of the results. To get the results for 
% each question, highlight the selected area and run it.

%% QUESTION 1 %%

clear
clc

N = 1000; % Number of observations.

beta0 = 0.4; % Beta 0.
beta1 = 0.5; % Beta 1.
beta2 = 0.1; % Beta 2.
beta3 = 0.6; % Beta 3.
beta4 = 0.2; % Beta 4.
beta5_extra_weak_instrument = 0.1 ; % Additional Instrument 

sigma_1 = 1; % Sigma 1.
sigma_2 = 2; % Sigma 2.
rho = 0.2;

mu = [0 0]; % Mean of error terms.
sigma = [sigma_1^2 rho*sigma_1*sigma_2; rho*sigma_1*sigma_2 sigma_2^2]; % Variance matrix of error terms.

sim = 100; % Number of simulations.

bhat_ols = zeros(sim,3); % OLS estimates.
bhat_iv = zeros(sim,3); % IV estimates.
bhat_liml = zeros(sim,3); % LIML estimates.

%% Estimation

x = ((1:N)'./N).*normrnd(0,1,N,1); % Exogenous regressor.
z = ((1:N)'./N).*normrnd(0,1,N,1); % Instrument.
z_2_weak_instru = ((1:N)'./N).*normrnd(0,1,N,1); % Instrument.

% Simulate

for s = 1:sim
  
    e = mvnrnd(mu,sigma,N); % Variance-covariance matrix of error terms.
    
    u = e(:,1); % Error terms for y.                                                     
    v = e(:,2); % Error terms for w.
    
    w = beta3 + beta4*z + v; % Endogenous regressor.
    w_many_instru = beta3 + beta4*z + beta5_extra_weak_instrument*z_2_weak_instru + v;

    y = beta0 + beta1*w + beta2*x + u;
    
    %% OLS estimates.
    x_ols = [ones(N,1) w x];
    bols = inv(x_ols'*x_ols)*x_ols'*y; % OLS estimate of beta.
    bhat_ols(s,:) = bols';
    
    %% IV estimates.
    x_fs = [ones(N,1) z x];
    what = x_fs*inv(x_fs'*x_fs)*x_fs'*w;
    
    x_ss = [ones(N,1) what x];
    biv = inv(x_ss'*x_ss)*x_ss'*y; % IV estimate of beta.
    bhat_iv(s,:) = biv';
    
    %% LIML estimates.
    
    z_liml = [ones(N,1) z ];
    w_liml = [ones(N,1) z x];
    
    Mz = eye(N) - z_liml*inv(z_liml'*z_liml)*z_liml';
    Mw = eye(N) - w_liml*inv(w_liml'*w_liml)*w_liml';
    
    Wz = [y w]'*Mz*[y w];
    Ww = [y w]'*Mw*[y w];
    kappa = min(eig((Ww^(-1/2))*Wz*(Ww^(-1/2))));
    bliml = inv(x_ols'*(eye(N)-(kappa*Mw))*x_ols)*(x_ols'*(eye(N)-(kappa*Mw))*y);
    
    bhat_liml(s,:) = bliml'; % LIML estimate of beta.
    
end

%% OLS

fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_ols(:,1)), std(bhat_ols(:,1)))

fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_ols(:,2)), std(bhat_ols(:,2)))

fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_ols(:,3)), std(bhat_ols(:,3)))

% Plot endog regressor

subplot(321);
histogram(bhat_ols(:,2), 10)  

% F-test for each simulation for First Stage:
OLS_independent_var_fs = [ones(N,1) z];
OLS_dependent_var_fs = w;
first_stage_OLS_para= inv(OLS_independent_var_fs'*OLS_independent_var_fs)*(OLS_independent_var_fs'*OLS_dependent_var_fs);
beta_3_fs_OLS = first_stage_OLS_para(1);
beta_4_fs_OLS = first_stage_OLS_para(2);
restricted(s) = sum((w-mean(w)).^2);
unrestricted(s) = sum((w-beta_3_fs_OLS-beta_4_fs_OLS*z).^2);

% OLS Extra Instruments

x_ols_extra_weak_instrument = [ones(N,1) w_many_instru x];
bols_extra_weak_instrument= inv(x_ols_extra_weak_instrument'*x_ols_extra_weak_instrument)*x_ols_extra_weak_instrument'*y;                                      % OLS estimate of beta.
bhat_ols_extra_weak_instrument(s,:) = bols_extra_weak_instrument';

% Graph
subplot(322);
histogram(bhat_ols_extra_weak_instrument(:,2), 10)

%% 2SLS

fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_iv(:,1)), std(bhat_iv(:,1)))

fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_iv(:,2)), std(bhat_iv(:,2)))

fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_iv(:,3)), std(bhat_iv(:,3)))

% Endog Reg

subplot(323);
histogram(bhat_iv(:,2), 10)

% IV Xtra Instru

x_fs_ew_instr = [ones(N,1) z z_2_weak_instru x];
what_ew_instr = x_fs_ew_instr*inv(x_fs_ew_instr'*x_fs_ew_instr)*x_fs_ew_instr'*w;
    
x_ss_ew_instr = [ones(N,1) what_ew_instr x];
biv_ew_instr = inv(x_ss_ew_instr'*x_ss_ew_instr)*x_ss_ew_instr'*y;                                          % IV estimate of beta.
bhat_iv_ew_instr(s,:) = biv_ew_instr';
    
subplot(324);
histogram(bhat_iv_ew_instr(:,2), 10)

%% LIML

fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_liml(:,1)), std(bhat_liml(:,1)))

fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_liml(:,2)), std(bhat_liml(:,2)))

fprintf(' %0.4f\n(%0.4f)\n', mean(bhat_liml(:,3)), std(bhat_liml(:,3)))

% Endog Reg

subplot(325);
histogram(bhat_liml(:,2), 10)

% Xtra Instru

z_liml_ew_instr = [ones(N,1) z z_2_weak_instru];
w_liml_ew_instr = [ones(N,1) z z_2_weak_instru x];
    
Mz_ew_instr = eye(N) - z_liml_ew_instr*inv(z_liml_ew_instr'*z_liml_ew_instr)*z_liml_ew_instr';
Mw_ew_instr = eye(N) - w_liml_ew_instr*inv(w_liml_ew_instr'*w_liml_ew_instr)*w_liml_ew_instr';
    
Wz_ew_instr = [y w]'*Mz_ew_instr*[y w];
Ww_ew_instr = [y w]'*Mw_ew_instr*[y w];
kappa_ew_instr = min(eig((Ww_ew_instr^(-1/2))*Wz_ew_instr*(Ww_ew_instr^(-1/2))));
bliml_ew_instr = inv(x_ols'*(eye(N)-(kappa*Mw_ew_instr))*x_ols_extra_weak_instrument)*(x_ols_extra_weak_instrument'*(eye(N)-(kappa*Mw_ew_instr))*y);
    
bhat_liml_extra_weak_instrument(s,:) = bliml_ew_instr'; 

subplot(326);
histogram(bhat_liml_extra_weak_instrument(:,2), 10)

%% F Test

F_test_FS =  mean(restricted./unrestricted);

%% QUESTION 2 %%

clear
clc

N = 100; % Number of observations.

beta0 = 0.2; % Beta 0.
beta1 = 0.5; % Beta 1.
beta2 = 0.3; % Beta 2.

sigma = 2; % Standard deviation.

dist = "Cauchy"; % Distribution of error terms.

bootsim = 200; % Number of "bootstrap repetitions".
sim = 200; % Number of simulations.

b_ols = zeros(sim,3); % OLS estimates.
se_ols = zeros(sim,3); % OLS standard errors.                
se_boot = zeros(sim,3);             

%% Simulate

x1 = 2 + ((1:N)'/N).*normrnd(0,1,N,1);
x2 = 3 + 0.5*x1 + normrnd(0,1,N,1);

x = [ones(N,1) x1 x2];

for s = 1:sim
    
    % Draw u and generate y.
    
    if dist == "Cauchy"    
        u = normrnd(0,sigma,N,1)./normrnd(0,sigma,N,1); % Errors that follows Cauchy dist.
    end
    
    y = beta0 + beta1*x1 + beta2*x2 + u;
    
    % OLS estimates.

    bhat = (x'*x)\x'*y; % OLS estimates.
    b_ols(s,:) = bhat';

    % OLS standard errors.

    uhat = y - x*bhat; % Residuals.
    s2 = (uhat'*uhat)/(N-3); % Sigma-squared.
    
    vmat = inv(x'*x).*s2;
    se_ols(s,:) = (diag(vmat).^(0.5))'; % OLS standard errors.
    
    % Non-parametric bootstrap standard errors.
    
    b_boot = zeros(bootsim,3);
    
    for i = 1:bootsim
        
        yb = x*bhat; % Fitted values.
        e = randsample(uhat,N,true); % Resample residuals.
            
        y_bs = yb+e; % New y.
        
        bhat_bs = (x'*x)\x'*y_bs; % OLS estimates from each bootstrap repetition.
        b_boot(i,:) = bhat_bs';
        
    end

    se_boot(s,:) = std(b_boot,0,1); % Bootstrap standard errors.

end

% Print Results
fprintf(' %0.4f\n(%0.4f)\n(%0.4f)\n(%0.4f)\n(%0.4f)\n(%0.4f)\n', mean(b_ols(:,1)), std(b_ols(:,1)), mean(se_ols(:,1)), std(se_ols(:,1)), mean(se_boot(:,1)), std(se_boot(:,1)))
fprintf(' %0.4f\n(%0.4f)\n(%0.4f)\n(%0.4f)\n(%0.4f)\n(%0.4f)\n', mean(b_ols(:,2)), std(b_ols(:,2)), mean(se_ols(:,2)), std(se_ols(:,2)), mean(se_boot(:,2)), std(se_boot(:,2)))
fprintf(' %0.4f\n(%0.4f)\n(%0.4f)\n(%0.4f)\n(%0.4f)\n(%0.4f)\n', mean(b_ols(:,3)), std(b_ols(:,3)), mean(se_ols(:,3)), std(se_ols(:,3)), mean(se_boot(:,3)), std(se_boot(:,3)))
