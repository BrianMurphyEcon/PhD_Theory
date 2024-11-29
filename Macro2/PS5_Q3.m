% Import the Data
clear
data = readtable('Macro 2\PS5\Data_PS5.xlsx');

%Question 3

consumption_growth = table2array(data(:,3)); 
equity_return = table2array(data(:,4));
treasury_return = table2array(data(:,5)); 

mean_cg = mean(consumption_growth); %mean consumption growth
std_cg = std(consumption_growth); %std dev consumption growth
auto_cg = autocorr(consumption_growth,'NumLags',1); %1-lag autocorr of cons growth

Pii = (1-auto_cg(2))/2; 
Pij = (1+auto_cg(2))/2;

P = [Pii Pij; Pij Pii]; %transition matrix

lambda1 = 1 + mean_cg + std_cg; 
lambda2 = 1 + mean_cg -std_cg;
lambda = [lambda1 lambda2];

n=20;

betas = linspace(0.9, 0.99,n);
sigmas = linspace(2,5,n);
deltas = linspace(0.5,1.5,n);

for i = 1:n
   for j =1:n
       for k =1:n

        beta =  betas(i);
        sigma = sigmas(j);
        delta = deltas(k);
        
        % Define function to be solved
        fun = @(p)[beta*(P(1,1)*((((lambda1-delta))/(1-(delta/lambda1)))^(-sigma))*(lambda1+p(1)) + P(1,2)*((((lambda2-delta))/(1-(delta/lambda1)))^(-sigma))*(lambda2+p(2)))-p(1);
        beta*(P(2,1)*((((lambda1-delta))/(1-(delta/lambda2)))^(-sigma))*(lambda1+p(1)) + P(2,2)*((((lambda2-delta))/(1-(delta/lambda2)))^(-sigma))*(lambda2+p(2)))-p(2)];

        % Initial guess
        x0 = [0, 0];

        % Solve using fsolve
        p = fsolve(fun, x0);

        % Display results
        p1 = p(1);
        p2 = p(2);

        %exp return
        
        ri = [P(1,1)*(lambda1)*(p(1)+1)/(p(1)) + P(1,2)*(lambda2)*(p(1)+1)/(p(1));
           P(2,1)*(lambda1)*(p(2)+1)/(p(2)) + P(2,2)*(lambda2)*(p(2)+1)/(p(2))];
        
        %risk free
        pf = beta.*[P(1,1)*(((lambda1-delta))/(1-(delta/lambda1)))^(-sigma)+P(1,2)*(((lambda2-delta))/(1-(delta/lambda1)))^(-sigma);
           P(2,1)*(((lambda1-delta))/(1-(delta/lambda2)))^(-sigma)+P(2,2)*(((lambda2-delta))/(1-(delta/lambda2)))^(-sigma)];
        
        rfi = [1/pf(1);
           1/pf(2)];
        
        %equity premium
        
        E(i,j,k) = mean(ri) - mean(rfi);
        R(i,j,k) = mean(ri);           

       end
   end
end

epd = mean(equity_return) - mean(treasury_return);
MSE = (E-(epd)).^2 + (R - mean(treasury_return)).^2;
mindiff = min(min(min(MSE)));
[b,s,d] = ind2sub(size(MSE), find(MSE == mindiff));


beta1= betas(b);

sigma1=sigmas(s);

delta1=deltas(d);