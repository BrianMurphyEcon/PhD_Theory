%% PS10

clc
clear

global y T 

T_values = [50,100,500]; 
beta1 = 1; 
beta1_2 = 0.5; %comment out for p1
sigma = 1; % std dev

sim = 100;                        
ols_mat = zeros(sim,3);                                      

average_gamma_hats = [T_values
    zeros(1,3)];

t_stat = [T_values
    zeros(1,3)];

for i = 1:length(T_values)
    T = T_values(i);

for s = 1:sim
    
    e = normrnd(0,sigma,T,1);
    u = normrnd(0,sigma,T,1);
    
    y = zeros(T,1);
    y(1) =  1;   

    x = zeros(T,1);
    x(1) =1;   

    
    for j = 2:T            
        y(j) = beta1*y(j-1) + e(j);   
        x(j) = beta1_2*y(j-1) + u(j);   
    end

    % OLS.
    
    X = [ones(T, 1), x]; 
    b_ols = (X'*X)\X'*y;
        
    ols_mat(s,1:2) = b_ols';

end

average_gamma_hats(2,i) = mean(ols_mat(:,2));

t_stat(2,i) = (mean(ols_mat(:,2)) - 1) / (std(ols_mat(:,2)) / sqrt(sim));

end