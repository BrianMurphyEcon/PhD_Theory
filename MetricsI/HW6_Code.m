clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Homework 6 Brian Murphy')
disp('Econometrics 1')
disp('Spring 2024')
disp('March 6, 2024')
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

% Create a vector of lags for both series. The command size(zzz,1) counts 
% the number of rows for vector/matrix zzz. The command www(a:b,1) takes 
% the ath until the bth element of the first column of a vector/matrix www.

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
%
% OLS Estimation and Application of Frisch-Waugh.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regress consumption growth on income growth and interest rate:
X = [ones(size(y,1),1) r y];

X = [ones(size(y,1),1) r y];

% OLS estimates. The inv command takes the inverse of a matrix if it
% exists.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Regress c on X. Name the vector of OLS estimates b.

b = inv(X'*X)*X'*c;
res = c-(X*b);
disp('OLS Coefficient Estimates:');
disp(b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2A

%take the absolute value of residuals, & regress them on each of the regressors.

abs_res = abs(res); %the absolute value of residual
abs_r_beta = (r'*r)\(r'*abs_res); %regressing residuals on r
abs_y_beta = (y'*y)\(y'*abs_res); %regressing residuals on y

% Heteroscedasticity test using t-tests
t_test_r = abs_r_beta/sqrt(inv(r'*r)/length(abs_res));
t_test_y = abs_y_beta/sqrt(inv(y'*y)/length(abs_res));

%Level of significance at 5%
alpha = 0.05; 

% Corrected dof
p_value_r = 2*(1 -tcdf(abs(t_test_r), length(c)-4)); 
p_value_y = 2*(1 -tcdf(abs(t_test_y), length(c)-4)); 

% Display results
disp(['T-statistic for r: ', num2str(t_test_r), ', p-value: ', num2str(p_value_r)]);
disp(['T-statistic for y: ', num2str(t_test_y), ', p-value: ', num2str(p_value_y)]);

% Check for heteroscedasticity
if p_value_r < alpha || p_value_y < alpha
    disp('Evidence of heteroscedasticity.');
else
      disp('No evidence of heteroscedasticity.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2B

% Plot residuals against each regressor
figure;

subplot(2,1,1);
scatter(r,res);
title('2B -Residuals vs. r');
xlabel('r');
ylabel('Residuals');

subplot(2,1,2);
scatter(y,res);
title('2B- Residuals vs. y');
xlabel('y');
ylabel('Residuals');


% Display results
disp('The scatterplot of residuals for both variables show no clear pattern. Therefore, the residuals are homoscedastic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2C
% Calculate the weights based on the square of the interest rates
weights = 1 ./(r.^2);

% Create a diagonal matrix with weights
W = diag(weights);

% GLS estimates
b_gls = inv(X'*W*X)*X'*W*c;

% Display GLS coefficient estimates
disp('GLS Coefficient Estimates:');
disp(b_gls);

% Calculate GLS residuals
residuals_gls = c-X*b_gls;

% Check for heteroscedasticity in GLS residuals
figure;
scatter(r, residuals_gls);
title('2C- GLS Residuals vs. r');
xlabel('r');
ylabel('GLS Residuals');

figure;
scatter(y, residuals_gls);
title('2C - GLS Residuals vs. y');
xlabel('r');
ylabel('GLS Residuals');

% Display results
disp('The scatterplot of residuals for both variables show a pattern, so I believe the residuals are heteroscedastic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2D
% OLS estimates
b_ols = inv(X'*X)*X'*c;

% Calculate OLS residuals
residuals_ols = c - X * b_ols;

% White robust standard errors
robust_covariance = inv(X'*X)*(X'*diag(residuals_ols.^2)*X)*inv(X'*X);

% White robust standard errors for each coefficient
robust_standard_errors = sqrt(diag(robust_covariance));

% Display OLS coefficient estimates and White robust standard errors
disp('OLS Coefficient Estimates:');
disp(b_ols);
disp('White Robust Standard Errors:');
disp(robust_standard_errors);


