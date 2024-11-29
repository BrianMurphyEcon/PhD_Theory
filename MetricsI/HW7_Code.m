clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Homework 7 Brian Murphy')
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


% Computing residuals of the regression

I = eye(length(c));
M = I - X*inv(X'*X)*X';

e = M*c;

% Regressing e_t on e_t-1

e_t = e(2:length(e)); 
e_t_minus1 = e(1:length(e) - 1); 
e_coef = inv(e_t_minus1'*e_t_minus1)*(e_t_minus1'*e_t); 

% T-Test

resid_e = e_t - e_coef*e_t_minus1; 

%unbiased estimator of sigma square

s_2 = (resid_e'*resid_e)/(length(e_t_minus1) - size(e_t_minus1,2)); 

% lets first create a matrix for (X'X)^(-1)
e_tminus1_inv = inv(e_t_minus1'*e_t_minus1); 

% create the t stats for the coefficients
t_stat = e_coef/sqrt(s_2*e_tminus1_inv(1,1));


% Degrees of freedom 
df = length(e_t_minus1) - size(e_t_minus1,2);  

% Calculate p-values for two-tailed tests
p_values = 2 * (1 - tcdf(abs(t_stat), df));

disp('P-Value for coefficient of e_tminus1');
disp(p_values);

%%%%%%% 3b)
%%%% Conduct Prais Winsten estimation

rho = e_coef; 

c_tilde = c(2:size(c)) - rho*c(1:size(c) - 1); 

X_tilde = X(2:end, :) - rho*X(1:end-1, :); % Calculated tilde variables for 
% the second observation onwards

%Now add PW correction for the first observation

c_tilde_PW = vertcat(sqrt(1 - rho^2)*c(1),c_tilde); 

X_tilde_PW = vertcat(sqrt(1 - e_coef^2)*X(1, :), X_tilde) ; 

b_PW = inv(X_tilde_PW'*X_tilde_PW)*(X_tilde_PW'*c_tilde_PW);

disp('   Estimates from PW estimation')
disp(b)
disp('Note: OLS estimates for b0, b1 and b2 with PW estimation in that order.')

e_PW = c_tilde_PW - X_tilde_PW*b_PW; 

%unbiased estimator of sigma square

s_2PW = (e_PW'*e_PW)/(length(c_tilde_PW) - size(X_tilde_PW,2)); 

% t stats for the regressors

% lets first create a matrix for (X'X)^(-1)
X_tilde_PW_inv = inv(X_tilde_PW'*X_tilde_PW); 

% create the t stats for the coefficients
t_b0PW = b_PW(1,1)/sqrt(s_2PW*X_tilde_PW_inv(1,1)); 

t_b1PW = b_PW(2,1)/sqrt(s_2PW*X_tilde_PW_inv(2,2));

t_b2PW = b_PW(3,1)/sqrt(s_2PW*X_tilde_PW_inv(3,3)); 


% Assuming 't_stats' is your vector of t-statistics
t_stats_PW = [t_b0PW,t_b1PW,t_b2PW];  % Replace this with your actual t-stats

% Degrees of freedom 
df_PW = length(c_tilde_PW) - size(X_tilde_PW,2);  

% Calculate p-values for two-tailed tests
p_values = 2 * (1 - tcdf(abs(t_stats_PW), df_PW));

disp('P-Values for b0, b1 and b2 with Prais-Winsten');
disp(p_values);


%Now doing Cochran-Orcutt

c_tilde_CO = c_tilde ;
X_tilde_CO = X_tilde ;

b_CO = inv(X_tilde_CO'*X_tilde_CO)*(X_tilde_CO'*c_tilde_CO);

disp('   Estimates from CO estimation')
disp(b)
disp('Note: OLS estimates for b0, b1 and b2 with CO estimation in that order.')

e_CO = c_tilde_CO - X_tilde_CO*b_CO; 

%unbiased estimator of sigma square

s_2_CO = (e_CO'*e_CO)/(length(c_tilde_CO) - size(X_tilde_CO,2)); 

% t stats for the regressors

% lets first create a matrix for (X'X)^(-1)
X_tilde_CO_inv = inv(X_tilde_CO'*X_tilde_CO); 

% create the t stats for the coefficients
t_b0_CO = b_CO(1,1)/sqrt(s_2_CO*X_tilde_CO_inv(1,1)); 

t_b1_CO = b_CO(2,1)/sqrt(s_2_CO*X_tilde_CO_inv(2,2));

t_b2_CO = b_CO(3,1)/sqrt(s_2_CO*X_tilde_CO_inv(3,3)); 


% Assuming 't_stats' is your vector of t-statistics
t_stats_CO = [t_b0_CO,t_b1_CO,t_b2_CO];  % Replace this with your actual t-stats

% Degrees of freedom 
df_CO = length(c_tilde_CO) - size(X_tilde_CO,2);  

% Calculate p-values for two-tailed tests
p_values = 2 * (1 - tcdf(abs(t_stats_CO), df_CO));

disp('P-Values for b0, b1 and b2 with Cochrane-Orcutt');
disp(p_values);
