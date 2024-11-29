clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Homework 5 - BM')
disp('Macro 2')
disp('Spring 2024')
disp('Feb 29, 2024')
disp('-------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the Data
clear
data = readtable('Macro 2\PS5\Data_PS5.xlsx');

beta = 0.99;
sigma = 10;

consumption_growth = table2array(data(:,3)); 
equity_returns = table2array(data(:,4)); 
treasury_returns = table2array(data(:,5)); 

mean_cg = mean(consumption_growth); 
std_cg = std(consumption_growth); 
auto_cg = autocorr(consumption_growth,'NumLags',1); 


Pii = (1-auto_cg(2))/2;
Pij = (1+auto_cg(2))/2;

P = [Pii Pij; Pij Pii];

lambda1 = 1 + mean_cg + std_cg;
lambda2 = 1 + mean_cg -std_cg;

lambda = [lambda1 lambda2];

 % Define function for fsolve
    fun = @(p) my_equations(p, beta, sigma, P, lambda);

    % Initial guess
    p0 = [0, 0];

    % Solve equations using fsolve
    p = fsolve(fun, p0);

    % Display results
    p1 = p(1);
    p2 = p(2);

% Return state j from i

rij = [(((1+p(1))*lambda(1))/p(1))-1 (((1+p(2))*lambda(2))/p(1))-1;
    (((1+p(1))*lambda(1))/p(2))-1 (((1+p(2))*lambda(2))/p(2))-1];

% Part a

ri = [P(1,1)*lambda(1)*(p(1)+1)/(p(1)) + P(1,2)*lambda(2)*(p(2)+1)/(p(1))-1;
   P(2,1)*lambda(2)*(p(1)+1)/(p(2)) + P(2,2)*lambda(2)*(p(2)+1)/(p(2))-1];

% Part b

pf = beta.*[P(1,1)*(lambda(1))^(-sigma)+P(1,2)*(lambda(2))^(-sigma);
   P(2,1)*(lambda(1))^(-sigma)+P(2,2)*(lambda(2))^(-sigma)];

rf_i = [1/pf(1)-1;
   1/pf(2)-1];


% Part c
ep = mean(ri) - mean(rf_i);
disp('eq prem')
disp(ep)

% Data equity premium
ep_data = mean(equity_returns) - mean(treasury_returns);
disp('data eq prem')
disp(ep_data)

disp('')
disp('In the above code we have calculated the expected return on equity, conditional on state, expected return of one period bond and the unconditional equity premium.');
disp('Interestingly, the results we get are different from what is equity return from the Mehra and Prescott paper. However, I suspect this is because they used annual data, while we used quarterly.')
function F = my_equations(p, beta, sigma, P, lambda)
    % Extract p1 and p2
    p1 = p(1);
    p2 = p(2);

    % Define equations
    eq1 = beta*(P(1,1)*((lambda(1))^-sigma)*(lambda(1)+p1) + P(1,2)*((lambda(2))^-sigma)*(lambda(2)+p2)) - p1;
    eq2 = beta*(P(2,1)*((lambda(1))^-sigma)*(lambda(1)+p1) + P(2,2)*((lambda(2))^-sigma)*(lambda(2)+p2)) - p2;

    % Return equations as a vector
    F = [eq1; eq2];
end

