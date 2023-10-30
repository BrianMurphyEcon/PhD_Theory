clear; % Clear variables from the workspace
clc; % Clear the command window
clear global;

% --------------------------------------------------------------------- %
%                             Problem 1.                               %
% --------------------------------------------------------------------- %

% Define Global
global alpha beta delta A_0 A_new k_0 kss T

% Set Parameters
alpha = 0.3;
beta = 0.96;
delta = 0.1;
A_0 = 15;
A_new = 35;
T = 150;
k_0 = A_0 * (alpha/(1/beta - (1-delta)))^(1/(1-alpha));
kss = A_new * (alpha/(1/beta - (1-delta)))^(1/(1-alpha));


% Initial Guess
guess_sol= linspace(k_0,kss,T);

% Solve
soln = fsolve(@fun, guess_sol);

k_val = zeros(1,T+2);
k_val(1) = k_0;
k_val(2:T+1) = soln;
k_val(T+2) = kss;

y_val = zeros(1,T+2);
c_val = zeros(1, T+2);
for i=1:T+2
    y_val(i)=(A_new^(1-alpha)*k_val(i)^alpha);
        if i<T+2 
            c_val(i) = y_val(i) - (k_val(i+1)-(1-delta)*k_val(i));
        else
            c_val(i) = y_val(i) - (kss - (1-delta)*k_val(i));
        end
end

figure;

t=0:151;
plot(t,k_val, 'b', 'DisplayName','Capital');
hold on;
plot(t,y_val, 'r', 'DisplayName','Output');
plot(t,c_val, 'k', 'DisplayName','Consumption');
hold off;

xlabel('Time');
ylabel('Capital,Output or Consumption Level');
legend('location', 'best');


p=table(t',k_val',y_val',c_val','VariableNames', {'time', 'capital', 'output', 'consumption'});

writetable(p, 'PS5_BrianMurphy_SSPath.xlsx', 'Sheet', 'ExportedData');


A1 = [A_0 A_new*ones(1,151)];
w_t=(1-alpha)*A1.^(1-alpha).*k_val.^alpha;
r_t=alpha*A1.^(1-alpha).*k_val.^(alpha-1)-delta;

figure;
t=0:151;
plot(t,w_t, 'b', 'DisplayName','Wage');
xlabel('Time');
ylabel('Wage');
legend('location', 'best');


figure;
plot(t,r_t, 'r', 'DisplayName','Interest_Rate');
hold off;

xlabel('Time');
ylabel('Interest Rate');
legend('location', 'best');


% --------------------------------------------------------------------- %
%                             Functions                                 %
% --------------------------------------------------------------------- %

function K = fun(k)

global alpha beta delta A_new k_0 kss T

% k_t
x1 = zeros(1,T);
x1(1) = k_0;
x1(2:T) = k(1:(T-1));

% k_{t+1}
x2 = k;

% k_{t+2}
x3 = zeros(1,T);
x3(1:(T-1)) = k(2:T);
x3(T) = kss;

K=((A_new^(1-alpha)*x2.^alpha-(x3-(1-delta)*x2))./(A_new^(1-alpha)*x1.^alpha-(x2-(1-delta)*x1)))-beta*(A_new^(1-alpha)*alpha*x2.^(alpha-1)+1-delta);

end