clear; % Clear variables from the workspace
clc; % Clear the command window
clear global;

% --------------------------------------------------------------------- %
%                             Problem 1a.                               %
% --------------------------------------------------------------------- %

% Clear variables and command window
clear;
clc;
clear global;

% Define global variables
global alpha beta r_star A_0 A_1 L_0 L_1 k_0

% Set the parameters
alpha = 1/3;
beta = 0.96;
r_star = 0.04/0.96;
A_0 = 10;
A_1 = 20;
L_0 = 10;
L_1 = 10;
k_0 = 1;

% Initial guess for the model
c_0_g = 1;
c_1_g = 1;
k_1_g = 1;
b_1_g = 1;

guess_sol = [c_0_g c_1_g k_1_g b_1_g];

% Solving the model

solution = fsolve(@fun3, guess_sol);

% Extracting the solution
c_0 = solution(1);
c_1 = solution(2);
k_1 = solution(3);
b_1 = solution(4);

% Using Solved Values to Solve for R1, W1, NX0 and NX1
R_1 = 1 + r_star;
W_1 = (1-alpha)*A_1^(1-alpha)*k_1^alpha;
NX_0 = b_1 * L_0;
NX_1 = -(1+r_star)*b_1*L_1;

% --------------------------------------------------------------------- %
%                             Functions                                 %
% --------------------------------------------------------------------- %

%% Part A
%x(1) = c_0;
%x(2) = c_1;
%x(3) = k_1;
%x(4) = b_1


function G = fun3(x)
global alpha beta r_star A_0 A_1  k_0

G(1) = x(1) + x(3) + x(4) - (A_0^(1-alpha)*k_0^(alpha));
G(2) = (1 + r_star) - (alpha * x(3)^(alpha-1) * A_1 ^(1-alpha));
G(3) = x(2) - (A_1^(1-alpha)*x(3)^alpha + (1+ r_star)*x(4));
G(4) = x(2) - (beta*(1+r_star)*(A_0^(1-alpha)*k_0-x(3)-x(4)));
end
