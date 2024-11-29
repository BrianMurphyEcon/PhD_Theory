clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Homework 8 Brian Murphy')
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

b=inv(X.'*X)*X.'*c;

disp(' ')
disp('Matlab command OLS estimate.')

disp('   Estimate')
disp(b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step b) IV Estimation
% First Stage: Regress current income growth on lagged income growth
iv_X = [ones(size(y,1), 1) ylag];
iv_b = inv(iv_X' * iv_X) * iv_X' * y;

% Get fitted values
fitted_y = iv_X * iv_b;

% Second Stage: Regress consumption growth on fitted values and interest rate
iv_X_second = [ones(size(y,1), 1) r fitted_y];
iv_b_second = inv(iv_X_second'*iv_X_second) * iv_X_second' * c;

% Step c) Try Different Instrument

iv_Xr = [ones(size(y,1), 1) r];
iv_br = inv(iv_Xr' * iv_Xr) * iv_Xr' * y;

% Get fitted values
fitted_y_r = iv_Xr * iv_br;

iv_X_second_r = [ones(size(y,1), 1) r fitted_y_r];
iv_b_second_r = inv(iv_X_second_r' * iv_X_second_r) * iv_X_second_r' * c;


% Step d) Calculate Standard Errors and Test Coefficient Significance
% Calculate standard errors of coefficients from IV estimation

% Note, Bent does not ask for this to be done on the other option of IV
% done in part c

iv_residuals = c - iv_X_second * iv_b_second;
iv_sigma_squared = iv_residuals' * iv_residuals / (size(iv_X_second, 1) - size(iv_X_second, 2));
iv_cov_matrix = iv_sigma_squared * inv(iv_X_second' * iv_X_second);
iv_standard_errors = sqrt(diag(iv_cov_matrix));

% Compare standard errors with those from OLS regression

% Hypothesis testing: Test if coefficient on income growth is zero
t_stat = iv_b_second(3) / iv_standard_errors(3);
p_value = 2 * (1 - tcdf(abs(t_stat), size(iv_X_second, 1) - size(iv_X_second, 2)));

% Print IV estimates, standard errors, and test results
disp('IV Estimates:');
disp(iv_b_second);
disp('IV Standard Errors:');
disp(iv_standard_errors);
disp(['t-statistic: ', num2str(t_stat)]);
disp(['p-value: ', num2str(p_value)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%