clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------------------------------------------')
disp('Homework 3 Brian Murphy')
disp('Econometrics 1')
disp('Spring 2024')
disp('Feb 7, 2024')
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

% Regress income growth on interest rate
X1 = [ones(size(y,1),1) r];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the Residual Maker, M

P = X*inv(X.'*X)*X.';

M = eye(18,18) - P;
eigM = eig(M, "matrix");

[C,A]=eig(M);
disp(A)

M_check = C.*A.*C.';

disp(C)
disp(A)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%