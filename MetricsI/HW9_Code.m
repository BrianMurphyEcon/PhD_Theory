clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Homework 9 Brian Murphy')
disp('Econometrics 1')
disp('Spring 2024')
disp('March 27, 2024')
disp('-------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code estimates the coefficients of the following model:
%          c(t) - c(t-1) = b0 + b1r(t) + b2(y(t) - y(t-1)) + u(t)
% using OLS. An application of the Frisch-Waugh Theorem is included.
%
% Income and consumption time series data (per capita; US dollars; US; 
% 1997-2015) are taken from the Bureau of Economic Analysis.
% 
% Real interest rate data (per cent; US; 1997-2015) are taken from the 
% World Bank.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load and format the data. Make sure the data is in the same working
% directory. Each variable is contained in separate excel files and will be 
% loaded as row vectors. The elements of each vector are observations from 
% 1997 to 2015 when read from left to right. 

consumption = xlsread('consumption.xls');
income = xlsread('income.xls');
interest = xlsread('interest.xls');

% Transpose the row vectors into column vectors. 

consumption = consumption';
income = income';
interest = interest';

% Take the log of consumption and income.

logc = log(consumption);
logy = log(income);


clag = logc(1:size(logc,1)-1,1);
ylag = logy(1:size(logy,1)-1,1);

% Take the first difference of each series. The vectors logc and logy will 
% be modified so that the dimensions match here.

logc = logc(2:size(logc,1));
logy = logy(2:size(logy,1));

c = logc - clag;
y = logy - ylag;

% Take the observations of r that correspond to the years that are not
% omitted.

r = interest(2:size(interest,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q3)
% OLS Estimation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate the model using OLS. Create the X matrix. The first column is a
% column of ones for the constant.

X = [ones(size(y,1),1) r y];

b = inv(X'*X)*(X'*c);

% Define the log-likelihood function
log_likelihood = @(params) -log_likelihood_regression(params, c, r, y);

% Estimate parameters using fminunc
initial_guess = [0, 0, 0, .5]; % Initial guess for coefficients and standard deviation
params_mle = fminunc(log_likelihood, initial_guess);

% Display estimated parameters
disp('MLE Coefficients and Standard Deviation:');
disp(params_mle);

%2b)

% Define the log-likelihood function
log_likelihood_het = @(params) -log_likelihood_regression_het(params, c, r, y);

% Estimate parameters using fminunc
initial_guess_het = [0, 0, 0, .5]; % Initial guess for coefficients and standard deviation
params_mle_het = fminunc(log_likelihood_het, initial_guess_het);

% Display estimated parameters
disp('MLE Coefficients and Standard Deviation with heteroskedasticity:');
disp(params_mle_het);

%2c)

[x,fval,exitflag,output,grad,hessian] = fminunc(log_likelihood, initial_guess);

% Estimate the variance-covariance matrix
variance_covariance_matrix = inv(hessian);

disp('Estimated Variance-Covariance Matrix:');
disp(variance_covariance_matrix);

%%% 2d)

% 2 stage GLS

e = c- X*b; % residuals

e_square = e.^2; 
r_square = r.^2; 

gamma = inv(r_square'*r_square)*r_square'*e_square; 

e_square_est = r_square*gamma; 
omega = diag(e_square_est); 

P = omega^(-1/2); 

c_tilde = P*c; % transformed dependent variable
X_tilde = P*X; % transformed regressors
b_fgls = inv(X_tilde'*X_tilde)*(X_tilde'*c_tilde);


disp('2d: Regression Results')

disp(' FGLS  Estimates')
disp(b_fgls)
disp('Note: FGLS estimates for b0, b1 and b2 in that order.')

var_cov_fgls = inv(X'*omega^(-1)*X);

%%
% Define the log-likelihood function for linear regression with normally distributed errors
function log_likelihood = log_likelihood_regression(params, y, x1, x2)
    % Extract parameters
    beta0 = params(1);
    beta1 = params(2);
    beta2 = params(3);
    sigma = params(4);

    % Compute predicted values
    y_pred = beta0 + beta1*x1 + beta2*x2;

    % Compute log-likelihood
    log_likelihood = -0.5 * sum(log(2*pi*sigma^2) + ((y - y_pred).^2) / sigma^2);
end

% Define the log-likelihood function for linear regression with normally distributed errors
function log_likelihood_het = log_likelihood_regression_het(params, y, x1, x2)
    % Extract parameters
    beta0 = params(1);
    beta1 = params(2);
    beta2 = params(3);
    sigma = params(4);

    % Compute predicted values
    y_pred = beta0 + beta1*x1 + beta2*x2;

     % Compute log-likelihood
    log_likelihood_het = -0.5 * sum(log(2*pi*sigma^2*x1.^2) + ((y - y_pred).^2) ./ (sigma^2*x1.^2));
end

