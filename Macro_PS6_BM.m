%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% value function iteration for simple non-stochastic growth model
% This is based on Russel Cooper program Chap. 5 grow.m
% Modified on 10/1/2020, Victor Zuluaga
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; % this clears variables from the workspace
clc; % this clears the command window


% Section 1: Functional Forms

% The model is the simple growth model
% preferences: u(c)=(c^(1-theta)-1)/(1-theta)
% discount factor: beta
% rate of capital depreciation: delta
% production function: A*k^(alpha)


% Section 2: Basic Parameters

theta = .75;
beta = .96;
delta = 0.05;
alpha = 1/2;
A = 1.2;


% Section 3: Spaces

% here you need to specify the space that the capital stock lives in
% you must make an educated guess here so that the capital stock space
% includes the steady state

% input the size of the state space
% this will be the number of points in the grid

N = 500;

% Determine the steady state capital stock 
% so here kstar is something we solved for analytically

kstar = (A*alpha/((1/beta)+delta-1))^(1/(1-alpha));

% put in bounds here for the capital stock space
% you may want to adjust these bounds

klo = kstar*.9;
khi = kstar*1.1;
step = (khi-klo)/(N-1);

% Now, construct the state space using the step starting with klo and
% ending with khi. By construction, the state space includes kstar (as
% it should be, right?)

k = (klo:step:khi)';
k1 = k;


% Section 4: Value Function Iteration

% Iterations: using this first guess at V iterate

T = 750;

% here we set our tolerance for convergence; 
% play with this to see how the results depend on tolerance

toler = .0001;

% Guess for the value function

V = zeros(N,1); % We set the initial value function to 0 for all values of K
v = V; % This will become the 'new' value function

matV = NaN(N,T+1); % Matrix to store the value function
matV(:,1) = V;

k_p = V; % This will be the policy function

tic
for j=1:T
    
    for i=1:N % for each value of k(t), find the value of k(t+1) that maximizes value function
		
		C = A*(k(i)^alpha) + (1-delta)*k(i) - k1; %k1 is next period's capital
		
		if theta == 1
            U = log(C);
        else
            U = (C.^(1-theta)-1)/(1-theta);
        end
        
		negc = find(C<0); % this line and the next line is a trick to keep consumption from going negative
		U(negc) = -1e50;
		
		r = U+beta*V;
		[vscal,kind] = max(r); % this is the actual maximization line
		v(i) = vscal; % store the values of the new value function
		k_p(i) = k(kind); % store the new policy function values
        
    end
    
    matV(:,j+1) = v; % Storing the new guess for the value function
    
    % Now, we have a new guess of the value function (v)
    % Now calculate the difference between the last (V) and current (v)
    % guess of the value function
    
    diff=(V-v)./V;
    
    % test that the difference is lower than the tolerance level

    if abs(diff) <= toler
        break
    else
    % so here if the diff exceeds the tolerance, then do a replacement
    % and go back to the top of the loop
        V=v;
    end

end

toc; % list the total time taken

% now that the program has converged.  You will have 
% the vector of values explicity: this is v.
% in addition, you have the policy function: ktp1.


% Section 5: Plots for the value and policy functions
% This displays the main results

% This plots the value function
close all
figure(1)
plot(k,matV(:,50),'k:',...
    k,matV(:,100),'k--',...
    k,matV(:,j),'k-');
xlabel('$$ k_{t} $$','interpreter','LaTex')
ylabel('$$ v(k_{t}) $$','interpreter','LaTex')
legend('50 iterations','100 iterations','Last iteration','Location','southeast')

% This plots the policy function against the 45 degree line
figure(2)
plot(k,k_p,'k-',...
    k,k,'r:');
hold on;
xlabel('$$ k_{t} $$','interpreter','LaTex')
ylabel('$$ k_{t+1} $$','interpreter','LaTex')
legend('Policy function','45-degree line','Location','southeast')


% Section 6: Use the policy function to plot the time series of tilde vars

T = 100; % Time periods to plot over
ts_k = zeros(T,1); % vector holding capital values at each time period
ts_k(1,1) = klo; % set initial value of capital in time period 1
ts_y = zeros(T,1); %vector holding output at each time period
ts_c = zeros(T,1); % vector holding consumption values at each time period
ts_R = zeros(T,1); % vector holding marginal product of capital at each time period
time = 1:1:T; % vector holding time periods

for t=1:T
    kdiff = abs(k - ts_k(t,1)); % find abs. value of differences of ts_k value from possible k values
    [kclose,index]=min(kdiff); % find location of minimum difference
    ts_k(t+1,1) = k_p(index); % Look up that location in policy function, which tells us next value of k to choose
    ts_y(t,1) = A*ts_k(t,1)^alpha; % output
    ts_c(t,1) = A*ts_k(t,1)^alpha + (1-delta)*ts_k(t,1) - ts_k(t+1,1); % consumption 
    ts_R(t,1) = A*alpha*ts_k(t,1)^(alpha-1); % marginal product of capital
end


% Section 7: Plots for the time series

% This plot capital stock over time
figure(3)
plot(time,ts_k(1:100))
xlabel('Time period')
ylabel('$$k_{t}$$','interpreter','LaTex')

% This plots output over time
figure(4)
plot(time,ts_y)
xlabel('Time period')
ylabel('$$y_{t}$$','interpreter','LaTex')

% This plots consumption over time
figure(5)
plot(time,ts_c)
xlabel('Time period')
ylabel('$$c_{t}$$','interpreter','LaTex')

% This plots MPK over time
figure(6)
plot(time,ts_R)
xlabel('Time period')
ylabel('$$R_{t}$$','interpreter','LaTex')
