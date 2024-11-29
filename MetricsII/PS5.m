%% Question 1

clear
clc

global x z

N = 300; % Number of observations.
beta_1 = 3; % Coefficient on X.
gamma_0 = -1; % First cutoff.
gamma_1 = 1; % Second cutoff.

if gamma_1 <= gamma_0
    error('gamma_0 < gamma.\n')
end

sigma = 1; % Standard deviation.

sim = 50; % Number of simulations.
results_mat = zeros(sim, 3); % Results matrix.

% MLE
x = ((1:N)'./N).*normrnd(0,1,N,1);
    
for s = 1:sim 

    % Generate the data
    u = normrnd(0,sigma,N,1);
    y = beta_1*x + u;
    
    z = double((gamma_0 < y));
    z = z + double((y >= gamma_1));
    
    % ML estimation.
    b0 = [0.5 -0.5 0.5];

    options = optimset('Display','off');
    [b_mle, ~, ~, ~, ~, hess] = fminunc(@logl_oprob, b0, options);

    results_mat(s, 1:size(results_mat,2)) = b_mle';
    vmat = inv(hess);

end

fprintf(' %0.4f\n(%0.4f)\n',b_mle(1),sqrt(vmat(1,1)))

fprintf(' %0.4f\n(%0.4f)\n',b_mle(2),sqrt(vmat(2,2)))

fprintf(' %0.4f\n(%0.4f)\n',b_mle(3),sqrt(vmat(3,3)))

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,1)), std(results_mat(:,1)))

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,2)), std(results_mat(:,2)))

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,3)), std(results_mat(:,3)))

% Plots
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
xlabel('$\gamma_{0}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\gamma_{0}$','interpreter','LaTex')
hold off

figure(3)
hold on
histogram(results_mat(:,3),bins)
xlabel('$\gamma_{1}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\gamma_{1}$','interpreter','LaTex')
hold off

%% Question 2

clear


global x T t N X Y 

N = 300; % Number of observations.
T = 100; % Number of periods in calendar time.
beta_0 = 0.5; % Beta 0.
beta_1 = 3; % Beta 1.

sim = 50; % Number of simulations.
results_mat_2 = zeros(sim, 2); % Results matrix.

x = ((1:N)'./N).*normrnd(0,1,N,1);
    
for s = 1:sim

    % Generate duration data from exponential distribution.
    theta = exp(beta_0 + beta_1*x);
    t = exprnd(1./theta,N,1);
    
    % MLE
    
    b0 = [0 0.5];

    options = optimset('Display','off');
    [b_mle_2, ~, ~, ~, ~, hess_2] = fminunc(@logl_dur, b0, options);

    results_mat_2(s,1:size(results_mat_2,2)) = b_mle_2';
    vmat_2 = inv(hess_2);

end

fprintf(' %0.4f\n(%0.4f)\n',b_mle_2(1),sqrt(vmat_2(1,1)))

fprintf(' %0.4f\n(%0.4f)\n',b_mle_2(2),sqrt(vmat_2(2,2)))

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat_2(:,1)), std(results_mat_2(:,1)))

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat_2(:,2)), std(results_mat_2(:,2)))

% Graphs

bins = 10;

figure(4)
hold on
histogram(results_mat_2(:,1),bins)
xlabel('$\beta_{0}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{0}$','interpreter','LaTex')
hold off

figure(5)
hold on
histogram(results_mat_2(:,2),bins)
xlabel('$\beta_{1}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{1}$','interpreter','LaTex')
hold off

%% Question 3

results_mat_3 = zeros(sim, 2);   % Results matrix.
rng('default')

x = ((1:N)'./N).*normrnd(0,1,N,1);
    
for s = 1:sim

    % Generate duration data from exponential distribution.

    theta = 1 ./ (1 + exp(-(beta_0 + beta_1 * x)));
    t = exprnd(1./theta);
    
    % MLE
    
    b0 = [0 0.5];

    options = optimset('Display','off');
    [b_mle_3, ~, ~, ~, ~, hess_3] = fminunc(@logl_dur, b0, options);

    results_mat_3(s,1:size(results_mat_3,2)) = b_mle_3';
    vmat_3 = inv(hess_3);

end

fprintf(' %0.4f\n(%0.4f)\n',b_mle_3(1),sqrt(vmat_3(1,1)))

fprintf(' %0.4f\n(%0.4f)\n',b_mle_3(2),sqrt(vmat_3(2,2)))

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat_3(:,1)), std(results_mat_3(:,1)))

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat_3(:,2)), std(results_mat_3(:,2)))

%Graphs 

bins = 10;

figure(6)
hold on
histogram(results_mat_3(:,1),bins)
xlabel('$\beta_{0}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{0}$, Q3','interpreter','LaTex')
hold off

figure(7)
hold on
histogram(results_mat_3(:,2),bins)
xlabel('$\beta_{1}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{1}$, Q3','interpreter','LaTex')
hold off

%% Question 4

%setup
results_mat_4 = zeros(sim, 2);    
t_i = t;


for i = 1:length(t_i)
    
    integer_part = floor(t_i(i));
    
    if integer_part >= 2
        y_i(i+ integer_part) = 1;
        y_i(i:i+integer_part-1) = 0;
        x_it(i:integer_part+i) = x(i);
    else 
        y_i(i) = 1;
        x_it(i) = x(i);
    end

end
Y = y_i';
X = x_it';


%estimation
beta_0 = 0.5;                                                     
beta_1 = 3; 

b0 = [1 1];                                                                 
options = optimset('Display','off');                                        
[b_mle_4, ~, ~, ~, ~, hess_4] = fminunc(@logl_logit_panel, b0, options);  


results_mat_4(s, 1:size(results_mat_4,2)) = b_mle_4';
vmat_4 = inv(hess_4);

disp(b_mle_4);

%% Functions

function log_likelihood = logl_oprob(beta)

    global x z

    beta1 = beta(1);
    gamma0 = beta(2);
    gamma1 = beta(3);
    mu = beta1 * x;

    log_likelihood(z == 0) = log(normcdf(gamma0-mu(z == 0)));
    log_likelihood(z == 1) = log(normcdf(gamma1-mu(z == 1)) - normcdf(gamma0-mu(z == 1)));
    log_likelihood(z == 2) = log(1-normcdf(gamma1-mu(z == 2)));
    log_likelihood = -sum(log_likelihood);


end

function Logl = logl_dur(b)

global x t N

mu = [ones(N,1) x] * b';
Logl = (t < 40) .* (log(exp(mu)) - exp(mu) .* t) - ( t >= 40) .* exp(mu) .* 40; 
Logl = -sum(Logl);

end

function L = logl_logit_panel(b)

global X Y

b0 = b(1);

b1 = b(2);

Xb = b0*ones(size(X,1),1)+b1*X;

L = -sum(Y .* log(cdf('Logistic',Xb)) + (1-Y) .* log(cdf('Logistic',-Xb)));

end
