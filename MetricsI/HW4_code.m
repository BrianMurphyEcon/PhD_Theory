clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Homework 4 Brian Murphy')
disp('Econometrics 1')
disp('Spring 2024')
disp('Feb 14, 2024')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part A. Calculate the T-Test
disp('Part A')

P = X*inv(X'*X)*X';
I = eye(18,18);
M = I - P;

% OLS Estimation
[b] = regress(c, X);

e = M*c; 

sigma_sq = (e'*e)/(length(c) - size(X,2));  
sigma = sqrt(sigma_sq);

% Calculate t-test statistics
Xinv = inv(X'*X); 

t_b0 = b(1,1)/(sqrt(sigma_sq*Xinv(1,1))); 

t_b1 = b(2,1)/(sqrt(sigma_sq*Xinv(2,2)));

t_b2 = b(3,1)/(sqrt(sigma_sq*Xinv(3,3)));


t_stat = [t_b0,t_b1,t_b2];   
% Degrees of freedom
n = length(c); % Number of observations
p = size(X, 2); % Number of regressors (including intercept)
df = n - p; % Degrees of freedom

% Calculate p-values
p_values = 2 * (1 - tcdf(abs(t_stat), df));

% Display the results
disp('P-Value');
disp(p_values);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part B.
disp('Part B')

ry=r+y;

X1= [ones(size(y,1),1) ry];

g = inv(X1'*X1)*X1'*c;

P1= X1*inv(X1'*X1)*X1';
M1=I-P1;

% Calculate REsiduals

e1=M1*c;

sigma_sq1 = (e1'*e1)/(length(c) - size(X1,2));  
sigma1 = sqrt(sigma_sq1);

% Calculate T-Stats


Xinv1 = inv(X1'*X1);

t_g0 = g(1,1)/(sqrt(sigma_sq1*Xinv1(1,1))); 

t_g1 = g(2,1)/(sqrt(sigma_sq1*Xinv1(2,2)));

% t-statistics
t_stats1 = [t_g0,t_g1];  

% Degrees of freedom 
df = length(c) - size(X1,2);   

p_values1= 2 * (1 - tcdf(abs(t_stats1), df));  

% F-test

R = [0,1,-1];
r = 0; 
j = 1;     

% F-statistic
F_stat_1 = ((R*b -r)'*inv(R*(sigma_sq1*inv(X'*X))*R')*(R*b - r))/(j); 


%Degrees of freedom

df_n = j; % degrees of freedom for numerator
df_d = length(c) - size(X1,2); % degrees of freedom for denominator
% P-Val for 2 tailed 
F_p_value = 2 * min(fcdf(F_stat_1, df_n, df_d), 1 - fcdf(F_stat_1, df_n, df_d));

disp('p-values')
disp(p_values1)

disp('P-value for H0 b1 = b2')
disp(F_p_value)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Part C')

R = [0,1,0;0,0,1];

r = [0,0]'; 

% number of restrictions
j = 2;        

% F-Statistic

F_stat2 = ((R*b -r)'*inv(R*(sigma_sq*inv(X'*X))*R')*(R*b - r))/(j);

% Calculate critical value for the F-test
crit_value_1 = finv(1 - 0.0250, 2, 15);

disp('F-stat for H0: b1 = 0, b2 = 0')
disp(F_stat2)

disp('Critical value')
disp(crit_value_1)

disp('We can reject the null hypothesis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Part D')

R = [1,-1,0;1,0,-1];

r = [0,0]'; 

j = 2; % number of restrictions

% F-statistic
F_stat3 = ((R*b -r)'*inv(R*(sigma_sq*inv(X'*X))*R')*(R*b - r))/(j);

% Calculate critical value for the F-test
crit_value2 = finv(1 - 0.0250, 2, 15);

disp('F-stat for H0: b1 = 0 and b2 = 0')
disp(F_stat3)

disp('Critical value')
disp(crit_value2)

disp('We can reject the null hypothesis.')


