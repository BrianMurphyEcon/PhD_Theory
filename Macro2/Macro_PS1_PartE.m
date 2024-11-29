% Section 2: Basic Parameters

theta = 1;
beta = .96;
delta = 0.05;
alpha = 0.36;
A = 1;


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

klo = kstar*.3;
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

% Calculate the real value function
A_real = (log(1 - alpha * beta) + ((alpha * beta) / (1 - alpha * beta) * log(alpha * beta))) / (1 - beta);
B_real = alpha / (1 - alpha * beta);
V_real = A_real + B_real * log(k);

% Calculate cstar
cstar = kstar^(alpha) - kstar;

% Calculate the analytical policy function
g_analytical = alpha * beta * k.^alpha;

% Plot results
figure;
plot(step, V, 'b-', 'LineWidth', 1.5);
hold on;
plot(step, V_real, 'r--', 'LineWidth', 1.5);
xlabel('Capital');
ylabel('Value Function');
title('Approximated vs Real Value Function');
legend({'Approximated', 'Real'}, 'Location', 'northwest', 'FontSize', 8);

% Plot the policy functions
figure;
plot(k, k_p, 'k-', 'LineWidth', 1.5);
hold on;
plot(k, g_analytical, 'b--', 'LineWidth', 1.5);
xlabel('Current Capital');
ylabel('Next Period Capital');
title('Real Policy Function vs Analytical Policy Function');
legend({'Real policy function', 'Analytical policy function'}, 'Location', 'southeast', 'FontSize', 8);


