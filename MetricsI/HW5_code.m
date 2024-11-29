clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Homework 5')
disp('Econometrics 1')
disp('Spring 2024')
disp('Feb 28, 2024')
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
% OLS Estimation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate the model using OLS. Create the X matrix. The first column is a
% column of ones for the constant.

X = [ones(size(y,1),1) r y];

% OLS estimates. The inv command takes the inverse of a matrix if it
% exists.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Regress c on X. Name the vector of OLS estimates b.

%Regression
b = inv(X'*X)*X'*c;

%Check

beta = regress(c,X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Splitting the data into two periods
% Assuming your data is stored in a MATLAB matrix named 'data' (an 18x1 matrix)

% Splitting the data into two periods
X_Before2008 = X(1:11, :);
X_After2008 = X(12:end, :);

beta_before2008 = inv(X_Before2008' * X_Before2008) * X_Before2008' * c(1:11);
beta_after2008 = inv(X_After2008' * X_After2008) * X_After2008' * c(12:end);

% Compute the sum of squared residuals for the combined model
SSR_combined = sum((c - X * b).^2);

% Compute the sum of squared residuals for the separate models
SSR_before2008 = sum((c(1:11) - X_Before2008 * beta_before2008).^2);
SSR_after2008 = sum((c(12:end) - X_After2008 * beta_after2008).^2);

% Compute degrees of freedom
n = size(c, 1);
k = 3;
%df_combined = n - k;
df_separate = n - 2*k;

% Compute the F-statistic
F = ((SSR_combined - (SSR_before2008 + SSR_after2008)) / 3) / ((SSR_before2008 + SSR_after2008) / df_separate);

% Set significance level
alpha = 0.05;

% Compare F-statistic to critical value
critical_value = finv(1 - alpha, 2, df_separate);
if F > critical_value
    disp('Reject null hypothesis: Coefficients are different before and after 2008.');
else
    disp('Fail to reject null hypothesis: Coefficients are the same before and after 2008.');
end

% Part B

% Extract interest rate data for the entire period and before/after 2008
r_total = X(:, 2);
r_before = X_Before2008(:, 2);
r_after = X_After2008(:, 2);

% Perform regressions
beta_r_total = inv(r_total' * r_total) * r_total' * c;
beta_r_before = inv(r_before' * r_before) * r_before' * c(1:11);
beta_r_after = inv(r_after' * r_after) * r_after' * c(12:end);

% Compute the sum of squared residuals for the combined model
SSR_combined_2 = sum((c - r_total * beta_r_total).^2);

% Calculate t-statistic for testing difference in coefficients
t_test2 = (beta_r_before - beta_r_after) / sqrt(SSR_combined_2 / (length(c) - 2 * 3));

% Calculate the p-value
pvalue_2 = 2 * tcdf(-abs(t_test2), length(c) - 2 * 3);

% Display results
disp([Since the p-value is greater than 0.05, we fail to reject the null that the coefficients are different.);
