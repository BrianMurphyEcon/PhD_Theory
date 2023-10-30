
clear; % Clear variables from the workspace
clc; % Clear the command window
clear global;

% --------------------------------------------------------------------- %
%                              Example 1                        %
% --------------------------------------------------------------------- %

% (b) Initial guess 
x0 = [1; 1; 1];

% (c) Solving the system of equations with fsolve
solution = fsolve(fun1, x0);

% Extracting the solution
f = solution(1);
g = solution(2);
h = solution(3);

% --------------------------------------------------------------------- %
%                             Example 2.                        %
% --------------------------------------------------------------------- %
% fsolve
clear; % Clear variables from the workspace
clc; % Clear the command window
clear global;

% (a) See files HW2fun.m 

% (b) Initial guess 
y0 = [1; 1];

% (c) Solving the system of equations with fsolve
solu = fsolve(@fun2, y0);

% Extracting the solution
x_1 = solu(1);
x_2 = solu(2);

% --------------------------------------------------------------------- %
%                             Example 3.                        %
% --------------------------------------------------------------------- %

% Clear variables and command window
clear;
clc;
clear global;

% Define global variables
global  A_0 A_1 k_0 k_1 c_0 c_1 beta alpha

% Set the parameters
alpha = 1/2;
A_0 = 10;
A_1 = A_0;
k_0 = 6;
beta = 0.95;

% Initial guess for the model
c_0_g = 1;
c_1_g = 2;
k_1_g = 1;

guess_sol = [c_1_g c_0_g k_1_g];

% Solving the model

solution = fsolve(@fun3, guess_sol);

% Extracting the solution
c_0 = solution(1);
c_1 = solution(2);
k_1 = solution(3);

% --------------------------------------------------------------------- %
%                             Functions                       %
% --------------------------------------------------------------------- %

%% Example 1

function F = fun1(x)
F(1) =  3*x(1) + x(2) - 2*x(3) - 6;
F(2) =  2*x(1) - x(2) + 4*x(3) - 0;
F(3) = - x(1) + 2*x(2) + 3*x(3) -7;
end

%% Example 2
function M = fun2(y)
M(1) = exp(-exp(-(y(1)+y(2)))) - y(2)*(1+y(1)^2);
M(2) = y(1)*cos(y(2)) + y(2)*sin(y(1)) - 0.5;
end

%% Example 3
%x(1) = c_0;
%x(2) = c_1;
%x(3) = k_1;

function G = fun3(x)
global A_0 A_1 k_0 beta alpha

G(1) =   x(1) + x(3) - k_0^alpha * A_0^(1-alpha);
G(2) = x(2) - beta * alpha * x(3)^(alpha-1) * A_1^(1-alpha) * x(1);
G(3) =  x(2) - x(3)^alpha * A_1^(1-alpha);
end
