clc;
clear;

%% Parameters
sigma = 2;
beta = 0.9;
r = 0.07;
mu = 10;
rho = 0.85; %0.95
epsilon = 1;
lambda = 3; 
states = 9;

tol = 1e-6; 
max_iter = 1000; 

agrid = 499; 
a_min = -10; %0;
a_max = 10; 
grid = (a_max - a_min) / agrid; 

a_mat = linspace(a_min, a_max, agrid+1)';

N = numel(a_mat); % total number of asset grid points

%% Tauchen Method

% Std of Discritized
sigma_y = epsilon / sqrt(1 - rho^2);

% Bounds for States
y_ub = mu + lambda*sigma_y;
y_lb = mu - lambda*sigma_y;

% State Space
state_space = linspace(y_lb, y_ub, states)';

% Midpoints
m = zeros(states-1, 1); 
for i = 1:(states-1)
    m(i) = (state_space(i+1) + state_space(i)) / 2;
end

transition_probs = zeros(states, states);
for i = 1:states
    for j = 1:states
        if j == 1
            transition_probs(i, j)= normcdf((m(1)-(1-rho)*mu -rho*state_space(i))/epsilon);
        elseif j == states
            transition_probs(i, j)= 1 -normcdf((m(states-1)-(1-rho)*mu-rho*state_space(i))/epsilon);
        else
            transition_probs(i, j)=normcdf((m(j)-(1-rho)*mu-rho*state_space(i))/epsilon)-normcdf((m(j-1) - (1-rho)*mu - rho*state_space(i)) / epsilon);
        end
    end
end

%% Value Function Iteration

% Initial guess for value function
v0 = zeros(N, states);

% Notes on "while" function: 
% https://sites.nd.edu/esims/files/2023/05/val_fun_iter_sp17.pdf

% Iterate on value function
its = 0;
dif = tol + 1;
while dif > tol & its < max_iter
    v1 = zeros(N, states);
    for j = 1:states
        for i = 1:N
            a0= a_mat(i);
            y0= state_space(j);
            c= a0 *(1 +r)+ y0- a_mat;
            c(c < 0) = NaN;
            util=(c.^(1 -sigma)-1)/(1-sigma);
            v1(i, j) =max(util + beta *(transition_probs(j,:)*v0(:,:)')'); 
        end
    end
    dif = norm(v1 - v0); %Difference between updated and previous
    v0 = v1;
    its = its + 1;
end

%% Policy Function

%Policy Matrices
policy_fn =zeros(agrid+1, states);
c_fn = zeros(agrid+1, states);


for j = 1:states
    for i = 1:N
        a0 = a_mat(i);
        y0 = state_space(j);
        c = a0 * (1 + r) + y0 - a_mat;
        c(c < 0) =NaN;
        util =(c.^(1 -sigma) -1) / (1 -sigma); 
        [~, idx] =max(util + beta*(transition_probs(j, :) * v0(:, :)')');
        policy_fn(i, j) = a_mat(idx);
        c_fn(i, j) = c(idx);
    end
end

%% Graph the Functions

figure;
hold on;
for j = 1:states
    plot(a_mat, v0(:, j), 'LineWidth', 1);
end
title('Value Function');
xlabel('Assets');
ylabel('Value');
legend('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8', 'State 9');
grid on;
hold off;

% Assets
figure;
hold on;
for j = 1:states
    plot(a_mat, policy_fn(:, j), 'LineWidth', 1);
end
title('Policy Function - Assets)');
xlabel('Assets');
ylabel('Policy');
legend('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8', 'State 9');
grid on;
hold off;

% Consumption
figure;
hold on;
for j = 1:states
    plot(a_mat, c_fn(:, j), 'LineWidth', 1);
end
title('Policy Function - Consumption');
xlabel('Assets');
ylabel('Consumption');
legend('State 1', 'State 2', 'State 3', 'State 4', 'State 5', 'State 6', 'State 7', 'State 8', 'State 9');
grid on;
hold off;


%% Simulation

obs_ps = 1000;
samples = 1000;

% Create Matrices
income_matrix = zeros(samples,obs_ps); 
asset_matrix = zeros(samples,obs_ps);
cons_matrix = zeros(samples,obs_ps);

% Income for Samples
initial_state = 1;
income_states = zeros(samples, obs_ps);
for sigma = 1:samples
    income_states(sigma, 1) = initial_state;
    for t = 2:obs_ps
        income_states(sigma, t) = randsample(1:9, 1, true, transition_probs(income_states(sigma, t-1), :));
    end
end

% Map state indices to income levels
income_levels=state_space';

for sigma = 1:samples
    for t = 1:obs_ps
        income_matrix(sigma, t)=income_levels(income_states(sigma,t));
    end
end

% Sim Assets
asset_vals_org=a_mat(randi(length(a_mat))); 
for sigma = 1:samples
    asset_matrix(sigma, 1)=asset_vals_org;
    for t = 2:obs_ps
        asset_matrix(sigma, t)=policy_fn(find(a_mat==asset_matrix(sigma,t-1)), income_states(sigma,t-1));
    end
end

% Sim Cons
for sigma = 1:samples
    for t = 1:obs_ps
        cons_matrix(sigma, t)=c_fn(find(a_mat==asset_matrix(sigma,t)), income_states(sigma,t));
    end
end

%% Statistics
num_samples = 1000;
num_periods = 1000;

income_mean = zeros(num_samples, 1);
income_std = zeros(num_samples, 1);
income_autocorrelation = zeros(num_samples, 1);

consumption_mean = zeros(num_samples, 1);
consumption_std = zeros(num_samples, 1);
consumption_autocorrelation = zeros(num_samples, 1);

for sample = 1:num_samples
    income_data = income_matrix(sample, :);
    cons_data = cons_matrix(sample, :);
    
    % mean, std, and autocorrelation for income
    income_mean(sample) = mean(income_data);
    income_std(sample) = std(income_data);
    i_autocorrelation_store = autocorr(income_data, 1);
    income_autocorrelation(sample) = i_autocorrelation_store(2);
    
    % mean, std, and autocorrelation for cons
    consumption_mean(sample) = mean(cons_data);
    consumption_std(sample) = std(cons_data);
    c_autocorrelation_store = autocorr(cons_data, 1);
    consumption_autocorrelation(sample) = c_autocorrelation_store(2);
end

% Calculate means and standard errors of the stats above
mean_income_mean = mean(income_mean);
mean_income_std_dev = mean(income_std);
mean_income_autocorr = mean(income_autocorrelation);

mean_consumption_mean =mean(consumption_mean);
mean_consumption_std_dev =mean(consumption_std);
mean_consumption_autocorr =mean(consumption_autocorrelation);

std_error_income_mean = std(income_mean) /sqrt(num_samples);
std_error_income_std_dev = std(income_std) /sqrt(num_samples);
std_error_income_autocorr = std(income_autocorrelation) /sqrt(num_samples);

std_error_consumption_mean = std(consumption_mean) /sqrt(num_samples);
std_error_consumption_std_dev = std(consumption_std) /sqrt(num_samples);
std_error_consumption_autocorr = std(consumption_autocorrelation) /sqrt(num_samples);

%% Create Table

% Define the data for the table
data = [mean_income_mean, mean_income_std_dev, mean_income_autocorr; 
        mean_consumption_mean, mean_consumption_std_dev, mean_consumption_autocorr;
        std_error_income_mean, std_error_income_std_dev, std_error_income_autocorr;
        std_error_consumption_mean, std_error_consumption_std_dev, std_error_consumption_autocorr];

% Define row names and column names for the table
rowNames = {'Mean Income', 'Mean Consumption', 'Std Error Income', 'Std Error Consumption'};
colNames = {'Mean', 'Std Dev', 'Autocorr'};

% Create a figure
fig = figure;
fig.Position = [100, 100, 500, 250]; % Set figure position and smaller size

% Create a uitable within the figure
uitable('Data', data, 'RowName', rowNames, 'ColumnName', colNames, ...
        'Position', [50, 50, 400, 150], ...  % Set uitable position and smaller size
        'ColumnWidth', {120}, 'FontSize', 10, 'FontWeight', 'bold', ...
        'Parent', fig);  % Specify the parent figure for the uitable

% Set title for the figure
fig.Name = 'Statistics Table';

% Graph has a scroll bar. To view entire graph, make sure to scroll.

%% Plot 150 Periods of Income, Assets and Consumption

% Sample from 1, first 150 periods

income_data_plot = income_matrix(1, 1:150);
consumption_data_plot = cons_matrix(1, 1:150);
asset_data_plot = asset_matrix(1, 1:150);

% Plot income
figure;
plot(1:150, income_data_plot, 'b', 'LineWidth', 1.5);
xlabel('Period');
ylabel('Income');
title('Income Over 150 Periods');

% Plot asset holdings
figure;
plot(1:150, asset_data_plot, 'r', 'LineWidth', 1.5);
xlabel('Period');
ylabel('Asset Holdings');
title('Asset Holdings Over 150 Periods');

% Plot consumption
figure;
plot(1:150, consumption_data_plot, 'g', 'LineWidth', 1.5);
xlabel('Period');
ylabel('Consumption');
title('Consumption Over 150 Periods');


%% Plot Stationary Dist of Wealth

obs_h = 100;
num_samples_h = 30000;
asset_matrix_h = zeros(num_samples_h,obs_h);

% Generate income observations for all samples
initial_state = 1;
income_state_h = zeros(num_samples_h, obs_h);
for sigma = 1:num_samples_h
    income_state_h(sigma, 1) = initial_state;
    for t = 2:obs_h
        income_state_h(sigma, t) = randsample(1:9, 1, true, transition_probs(income_state_h(sigma, t-1), :));
    end
end

asset_values_h = a_mat(randi(length(a_mat))); 
for sigma = 1:num_samples_h
    asset_matrix_h(sigma, 1) = asset_values_h;
    for t = 2:obs_h
        asset_matrix_h(sigma, t) = policy_fn(find(a_mat == asset_matrix_h(sigma,t-1)), income_state_h(sigma,t-1));
    end
end

for q= 3:width(asset_matrix_h)

m_t = [mean(asset_matrix_h(:,q)) var(asset_matrix_h(:,q)) ];
m_t_1 = [mean(asset_matrix_h(:,q-1)) var(asset_matrix_h(:,q-1))];

diff_m= max(abs(m_t-m_t_1));

if diff_m <= 0.001
    break
end

end

figure
histogram(asset_matrix_h(:,q));
title('Stationary Wealth Distribution');
