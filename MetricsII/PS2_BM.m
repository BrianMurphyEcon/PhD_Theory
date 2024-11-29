%set the parameters
clc
clear

global x w y z

N = 100; % Number of observations.
beta0 = 0.5; % Beta 0.
beta1 = 3; % Beta 1.
gamma0 = 1; % Gamma 0.
gamma1 = 2; % Gamma 1.
sigmau = 2; % Standard deviation for u.
rho = -0.9; % Correlation between u and v.
sim = 50; % Number of simulations.
results_mat = zeros(sim,3);                      

%MLE estimation
x = ((1:N)'/N).*normrnd(0,1,N,1);
w = ((1:N)'/N).*normrnd(0,1,N,1);
    
for s = 1:sim 
    e = mvnrnd([0; 0],[sigmau^2 rho*sigmau; rho*sigmau 1],N); % Correlated error terms. Standard deviation of v is 1.
    u = e(:,1); % Error terms for y.                                                     
    v = e(:,2); % Error terms for q.
    
    y0 = beta0 + beta1*x + u; % Latent variable y.
    z0 = gamma0 + gamma1*w + v; % Latent variable z.
    
    y = y0.*double((z0 > 0)); % Observed y.
    y((z0 <= 0)) = 0; % Truncation based on z0.
    
    z = double((z0 > 0)); % Observed z.
    
    options = optimoptions(@fminunc,'Display','off');
    
    b0 = [0 0 1]; % Initial values.
    [b_mle, ~, ~, ~, ~, hess] = fminunc(@fun_linear, b0, options); % Minimization.

    results_mat(s,:) = b_mle'; % Store results.
    vmat = inv(hess);
end

fprintf(' %0.4f\n(%0.4f)\n',b_mle(1),sqrt(vmat(1,1)))
fprintf(' %0.4f\n(%0.4f)\n',b_mle(2),sqrt(vmat(2,2)))
fprintf(' %0.4f\n(%0.4f)\n',b_mle(3),sqrt(vmat(3,3)))
fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,1)), std(results_mat(:,1)))
fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,2)), std(results_mat(:,2)))
fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,3)), std(results_mat(:,3)))

bins = 10;
close all
figure(1)
hold on
histogram(results_mat(:,1),bins)
xlabel('$\beta_{0}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{0}$','interpreter','LaTex')
hold off

figure(2)
hold on
histogram(results_mat(:,2),bins)
xlabel('$\beta_{1}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{1}$','interpreter','LaTex')
hold off

figure(3)
hold on
histogram(results_mat(:,3),bins)
xlabel('$\sigma$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\sigma$','interpreter','LaTex')
hold off

sim = 50; 
results_mat_selection = zeros(sim,6); 


for s = 1:sim
    e = mvnrnd([0; 0],[sigmau^2 rho*sigmau; rho*sigmau 1],N); % Correlated error terms. Standard deviation of v is 1.
    u = e(:,1); % Error terms for y.                                                     
    v = e(:,2); % Error terms for q.
    
    y0 = beta0 + beta1*x + u; % Latent variable y.
    z0 = gamma0 + gamma1*w + v; % Latent variable z.
    
    y = y0.*double((z0 > 0)); % Observed y.
    y((z0 <= 0)) = 0; % Truncation based on z0.
    
    z = double((z0 > 0)); % Observed z.
    
    options = optimoptions(@fminunc,'Display','off');
    
    b0 = [0.5 0.5 0.5 0.5 0.5 0.01]; % Initial values.
    [b_mle_selection, ~, ~, ~, ~, hess_selection] = fminunc(@fun_selection, b0, options); % Minimization.

    results_mat_selection(s,:) = b_mle_selection'; % Store results.
    vmat_selection = inv(hess_selection);

end

fprintf(' %0.4f\n(%0.4f)\n',b_mle_selection(1),sqrt(vmat_selection(1,1)))
fprintf(' %0.4f\n(%0.4f)\n',b_mle_selection(3),sqrt(vmat_selection(3,3)))
fprintf(' %0.4f\n(%0.4f)\n',b_mle_selection(2),sqrt(vmat_selection(2,2)))
fprintf(' %0.4f\n(%0.4f)\n',b_mle_selection(4),sqrt(vmat_selection(4,4)))
fprintf(' %0.4f\n(%0.4f)\n',b_mle_selection(6),sqrt(vmat_selection(6,6)))
fprintf(' %0.4f\n(%0.4f)\n',b_mle_selection(5),sqrt(vmat_selection(5,5)))
fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,1)), std(results_mat(:,1)))
fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat_selection(:,3)), std(results_mat_selection(:,3)))
fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,2)), std(results_mat(:,2)))
fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat_selection(:,4)), std(results_mat_selection(:,4)))
fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat_selection(:,6)), std(results_mat_selection(:,6)))
fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat_selection(:,5)), std(results_mat_selection(:,5)))

figure(4)
hold on
histogram(results_mat_selection(:,1),bins)
xlabel('$\beta_{0}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{0}$','interpreter','LaTex')
hold off

figure(5)
hold on
histogram(results_mat_selection(:,3),bins)
xlabel('$\gamma_0$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\gamma_0$','interpreter','LaTex')
hold off

figure(6)
hold on
histogram(results_mat_selection(:,2),bins)
xlabel('$\beta_{1}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{1}$','interpreter','LaTex')
hold off

figure(7)
hold on
histogram(results_mat_selection(:,4),bins)
xlabel('$\gamma_1$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\gamma_1$','interpreter','LaTex')
hold off

figure(8)
hold on
histogram(results_mat_selection(:,6),bins)
xlabel('$\rho$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\rho$','interpreter','LaTex')
hold off

figure(9)
hold on
histogram(results_mat_selection(:,5),bins)
xlabel('$\sigma_u$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\sigma_u$','interpreter','LaTex')
hold off

function [L] = fun_linear(b)

% log likelihood function for the linear model-regression of y on x- without accounting for censoring.
global x y

b0 = b(1);
b1 = b(2);
sigma = b(3);

X_B = b0 + b1*x;

L = - 0.5*log(2*pi) - 0.5*log(sigma^2) - 0.5.*(((y - X_B).^2)./(sigma^2));
L = sum(L);
L = -L;

end

function [L] = fun_selection(b)

% loglikelihood function for a model with selection.

global x w y z

b0 = b(1);
b1 = b(2);

g0 = b(3);
g1 = b(4);

sigma = b(5);
rho   = b(6);

XB = b0 + b1*x ;
WG = g0 + g1*w ;

L = log((1/sigma).*normpdf((y - XB)./sigma)) + log(normcdf((1/sqrt(abs(1-rho^2))).*(WG + (rho/sigma).*(y - XB))));
L(z==0) = log(normcdf(-WG(z==0)));
L = sum(L);
L = -L;

end



