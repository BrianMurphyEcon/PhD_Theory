clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Homework 1 Brian Murphy')
disp('Econometrics 1')
disp('Spring 2024')
disp('Jan 24, 2024')
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

% Estimate the model using OLS. Create the X matrix. The first column is a
% column of ones for the constant.

X = [ones(size(y,1),1) r y];

% OLS estimates. The inv command takes the inverse of a matrix if it
% exists.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Regress c on X. Name the vector of OLS estimates b.
% c is the dependent variable, X is independent

b=inv(X.'*X)*X.'*c;
check1 = regress(c,X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display the results in a table.

disp(' ')
disp('Part A')
disp(' ')

disp('Model: c(t) - c(t-1) = b0 + b1r(t) + b2(y(t) - y(t-1)) + u(t)')
disp(' ')

disp('Regression Results')

disp('   Estimates')
disp(b)
disp(check1)
disp('Note: OLS estimates for b0, b1 and b2 in that order.')

% Matlab also has a built in command for obtaining the OLS coefficients.
% You can use this to double check the results.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Use qq = regress(pp,kk) to find the estimates qq for a regression
% of pp on kk. Take pp as c and kk as X. Label the vector of estimates as 
% beta.

beta = regress(c,X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Matlab command OLS estimates.')

disp('   Estimates')
disp(beta)
disp('Note: OLS estimates for b0, b1 and b2 in that order.')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Regress income on the interest rate. The regression includes the 
% intercept. Compute and record the residuals.

X1 = [ones(size(y,1),1) r];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Regress y on X1. Label the vector of OLS estimates g. Compute the
% residuals and label these Minc.
g = inv(X1.'*X1)*X1.'*y;
check2 = regress(y, X1);
yhat = X1*g;
Minc = y - yhat;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp(' ')
disp('Part B')
disp(' ')

disp('Model: y(t) = g0 + g1r(t) + v(t)')
disp(' ')

disp('Regression Results')

disp('   Estimates')
disp(g)
disp('Note: OLS estimates for g0 and g1 in that order.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Use qq = regress(pp,kk) to find the estimates qq for a regression
% of pp on kk. Take pp as y and kk as X1. Label the vector of estimates as 
% gamma.
gamma = regress(y,X1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Matlab command OLS estimate.')

disp('   Estimate')
disp(gamma)
disp('Note: OLS estimates for g0 and g1 in that order.')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Regress consumption on the interest rate and the residuals from Part B.
% The regression includes the intercept term.

X3 = [ones(size(y,1),1) r Minc];

disp('-------------------------------------------------------------------')
disp(' ')
disp('Part C')
disp(' ')

disp('Model: c(t) - c(t-1) = a0 + a1r(t) + a2Minc(t) + e(t)')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Regress c on X3. Label the vector of OLS estimates A1.
A1 = inv(X3.'*X3)*X3.'*c;
check3 = regress(c,X3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Regression Results')

disp('   Estimates')
disp(A1)
disp('Note: OLS estimates for a0, a1 and a2 in that order.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Use qq = regress(pp,kk) to find the estimates qq for a regression
% of pp on kk. Take pp as c and kk as X3. Label the vector of estimates as 
% alpha1.
alpha1 = regress(c,X3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Matlab command OLS estimates.')

disp('   Estimate')
disp(alpha1)
disp('Note: OLS estimates for a0, a1 and a2 in that order.')
disp(' ')

% Regress consumption on just the residuals from Part B.

disp('-------------------------------------------------------------------')
disp(' ')

X2 = Minc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Regress c on X2. Label the OLS estimate A2.
A2 = inv(X2.'*X2)*X2.'*c;
check4=regress(c,X2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Model: c(t) - c(t-1) = a2Minc(t) + r(t)')
disp(' ')

disp('Regression Results')

disp('   Estimates')
disp(A2)
disp('Note: OLS estimate for a2.')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To do: Use qq = regress(pp,kk) to find the estimates qq for a regression
% of pp on kk. LTake pp as c and kk as X2. abel the vector of estimates as 
% alpha2.

alpha2=regress(c,X2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Matlab command OLS estimates.')

disp('   Estimate')
disp(alpha2)
disp('Note: OLS estimate for a2.')
disp(' ')

disp('-------------------------------------------------------------------')
disp(' ')

disp('Analysis: ')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('I can confirm that the coefficients to income/Minc are the same in all the regressions.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%