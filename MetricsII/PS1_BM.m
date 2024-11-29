%% Probit 

%1. Set the parameters

close all
clear
clc
global x z N

N = 300; % Number of observations.
beta0 = 0.5; % Intercept.
beta1 = 3; % Coefficient on X.
sigma = 1; % Standard deviation.
sim = 100; % Number of simulations.

results_mat = zeros(sim, 2); 

%2. Maximum Likelyhood Estimation 

x = [ones(N, 1), ((1:N)'./N).*normrnd(0,1,N,1)]; % Generate X; 
beta = [beta0; beta1]; % Vector of coefficients
for s = 1:sim

    % Generate the data.
    u = normrnd(0,sigma,N,1); % Generate U.
    y = x*beta + u; % Generate Y.
    z = double((y > 0)); % Generate Z.

    %Likelihood Probit
    logl_probit = @(beta) -sum(z .* log(normcdf(x * beta)) + (1 - z) .* log(1 - normcdf(x * beta)));

    % Estimation using ML.

    b0 = [1; 1]; % Initial values.
    options = optimset('Display','off'); % Turn off display.
    [b_mle, ~, ~, ~, ~, hess] = fminunc(logl_probit, b0, options); % Minimization
    b0 = b_mle;

    % Store estimates.
    results_mat(s, 1:size(results_mat,2)) = b_mle'; % Store results.
    vmat = inv(hess);
end 

%3. Display results of simulation 
fprintf('%0.4f\n(%0.4f)\n',b_mle(1),sqrt(vmat(1,1)))
fprintf('%0.4f\n(%0.4f)\n',b_mle(2),sqrt(vmat(2,2)))

%4. Empirical Results 
fprintf('%0.4f\n(%0.4f)\n', mean(results_mat(:,1)), std(results_mat(:,1)))
fprintf('%0.4f\n(%0.4f)\n', mean(results_mat(:,2)), std(results_mat(:,2)))

%Plots.

%b_0 
figure(1)
hold on
histogram(results_mat(:,1),5)
xlabel('$\beta_{0}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{0}$','interpreter','LaTex')
hold off

%b_1
figure(2)
hold on
histogram(results_mat(:,2),5)
xlabel('$\beta_{1}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{1}$','interpreter','LaTex')
hold off

%% Logit 
% Using the same global variables as last time.

results_mat_logit = zeros(sim, 2);

for s = 1:sim

    % Generate the data.

    u = normrnd(0,sigma,N,1); % Generate U.
    y = x*beta + u; % Generate Y.
    z = double((y > 0)); % Generate Z.

    %Likelihood Logit

    logl_logit = @(beta) -sum(z .* log(1 ./ (1 + exp(-x * beta))) + (1 - z) .* log(1 - 1 ./ (1 + exp(-x * beta))));

    % Estimation using ML.

    b0 = [1; 1]; % Initial values.
    options = optimset('Display','off'); % Turn off display.
    [b_mlegit, ~, ~, ~, ~, hess] = fminunc(logl_logit, b0, options); % Minimization.
    b0 = b_mlegit; % Updating guess

    % Store estimates.

    results_mat_logit(s, 1:size(results_mat_logit,2)) = b_mlegit'; % Store results.
    vmat_logit = inv(hess);
end 

%3. Display results of simulation 
fprintf('%0.4f\n(%0.4f)\n',b_mlegit(1),sqrt(vmat_logit(1,1)))
fprintf('%0.4f\n(%0.4f)\n',b_mlegit(2),sqrt(vmat_logit(2,2)))

%4. Empirical Results 
fprintf('%0.4f\n(%0.4f)\n', mean(results_mat_logit(:,1)), std(results_mat_logit(:,1)))
fprintf('%0.4f\n(%0.4f)\n', mean(results_mat_logit(:,2)), std(results_mat_logit(:,2)))

%Plots

%b_0 
figure(1)
hold on
histogram(results_mat_logit(:,1),5)
xlabel('$\beta_{0}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{0}$','interpreter','LaTex')
hold off

%b_1
figure(2)
hold on
histogram(results_mat_logit(:,2),5)
xlabel('$\beta_{1}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{1}$','interpreter','LaTex')
hold off