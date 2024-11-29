% Import the Data
data = readtable('Macro 2\PS5\Data_PS5.xlsx');

% Seperate the data
consumption_growth = data{:, 3};
equity_returns = data{:, 4};
treasury_returns = data{:, 5};

% Calculate mean
mean_cg = mean(consumption_growth);

% Calculate standard deviation
std_cg = std(consumption_growth);

% Calculate autocorrelation
autocorr_cg = autocorr(consumption_growth, 'NumLags', 1);

Pii = (1-autocorr_cg(2))/2;
Pij = (1+autocorr_cg(2))/2;

% Create Transition Matrix
P = [Pii, Pij; Pij, Pii];

mu1 = 1+mean_cg+std_cg;
mu2 = 1+mean_cg-std_cg;

mu = [mu1, mu2];


% Part A
beta = 0.99;
sigma=10;

syms p1 p2

eqn1 = beta*(P(1,1)*((mu1)^(-sigma))*(mu1+p1)+P(1,2)*((mu2)^(-sigma))*(mu2+p2))-p1;
eqn2 = beta*(P(2,1)*((mu1)^(-sigma))*(mu1+p1)+P(2,2)*((mu2)^(-sigma))*(mu2+p2))-p2;

sols= solve ([eqn1, eqn2], [p1,p2]);
p1sol = sols.p1;
p2sol = sols.p2;
p1 = double(p1sol);
p2 = double(p2sol);





% Part B

% Part C


% Question 2

%Set Up Parameters

n = 50;
beta_2 = linspace(0.85, 0.99, n);
sigma_2 = linspace(2, 4, n);
delta_2 = linspace(0.5, 0.9, n);

P = [Pii, Pij; Pij, Pii];

for i = 1:n
    for j = 1:n
        for k = 1:n


real_beta = beta_2(i);
real_sigma = sigma_2(j);
real_delta = delta_2(k);

syms p1_2 p2_2

eqn1 = beta*(P(1,1)*((mu1-delta_2)/(1-);
eqn2 = beta*(P(2,1)*((mu1);

sols_2= solve ([eqn1, eqn2], [p1_2,p2_2]);
p1sol_2 = sols_2.p1_2;
p2sol_2 = sols_2.p2_2;
Rp1_2 = double(p1sol_2);
Rp2_2 = double(p2sol_2);


        end
    end
end









