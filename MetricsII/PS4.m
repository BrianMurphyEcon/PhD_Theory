%% Question 1

clear
clc

p = 1; % pth difference.
k = 1; % kth lag.

%% Data

load data

% Label data.

pop = [pops6396(1:8,:); pops6396(10:51,:)];                             
pop = pop(:,1:33); % Population  
cpi_vector = cpi6396(1:33);                                             
cpi     = kron(ones(size(pop,1),1),cpi_vector'); % CPI
inc = [dpi6396(1:8,:); dpi6396(10:51,:)];                               
inc = inc(:,1:33); % Income
agg_inc = sum(inc)';                                                      
inc   = inc./cpi;                  
agg_inc = agg_inc./cpi_vector; % Aggregate Income
ndur = ndur6095(:,4:36);                                               
non_dur     = [ndur(1:8,:); ndur(10:51,:)]; % Non-Durable Goods
sale_pc = non_dur./cpi;                                                           
perc = sum(non_dur)';                                                          
perc = perc./cpi_vector;                                                      

clear non_dur cpi cpi_vector pops6396 cpi6396 dpi6396 ndur6095 ndur

% Make everything per capita.

dpi     = inc./pop;
sale_pc    = sale_pc./pop;
perc    = perc./sum(pop)';
agg_inc = agg_inc./sum(pop)';
agg_inc = kron(ones(size(pop,1),1),agg_inc');  
perc    = kron(ones(size(pop,1),1),perc');  

clear pop

% Take logs

log_dpi     = log(dpi);
log_sale    = log(sale_pc);
log_perc    = log(perc);
log_dpi_agg = log(agg_inc);

clear dpi sale_pc perc agg_inc

% Differences

D = size(log_dpi,2);
diff_dpi     = my_diff(log_dpi,p,D);
diff_sale    = my_diff(log_sale,p,D);
diff_perc    = my_diff(log_perc,p,D);
diff_dpi_agg = my_diff(log_dpi_agg,p,D);

clear D p log_dpi log_sale log_perc log_dpi_agg

dpi_tp = diff_dpi';
dpi_1ag1 = lagmatrix(diff_dpi',1);
n1 = size(dpi_tp,1);
n2 = size(dpi_1ag1,1);
dpi_tp = dpi_tp(k+1:n1,:)';
dpi_1ag1 = dpi_1ag1(k+1:n2,:)';

clear n1 n2 diff_dpi

sale_tp = diff_sale' ;
n3 = size(sale_tp,1);
sale_tp = sale_tp(k+1:n3,:)';

clear n3 diff_sale

dpi_agg_tp = diff_dpi_agg';
dpi_agg_1ag1 = lagmatrix(diff_dpi_agg',1);
n4 = size(dpi_agg_tp,1);
n5 = size(dpi_agg_1ag1,1);
dpi_agg_tp = dpi_agg_tp(k+1:n4,:)';
dpi_agg_1ag1 = dpi_agg_1ag1(k+1:n5,:)';

clear n4 n5 diff_dpi_agg

perc_tp = diff_perc';

n6 = size(perc_tp,1);

perc_tp = perc_tp(k+1:n6,:)';

clear n6 diff_perc

%% Estimate by FGLS %%

N = size(sale_tp,1);
T = size(sale_tp,2);
W = ones(N,T);

% Time fixed-effects

[dpi_t_ft,dpi_1_ft] = my_fe(dpi_tp,dpi_1ag1,1,0,W);
[dpi_atft,dpi_a1ft] = my_fe(dpi_tp-dpi_agg_tp,dpi_1ag1-dpi_agg_1ag1,1,0,W);
[saleatft,~] = my_fe(sale_tp-perc_tp,perc_tp,1,0,W);
[sale_tft,perc_tft] = my_fe(sale_tp,perc_tp,1,0,W);

% Cross section fixed-effects

[dpi_t_fx,dpi_1_fx] = my_fe(dpi_tp,dpi_1ag1,0,1,W);
[dpi_atfx,dpi_a1fx] = my_fe(dpi_tp-dpi_agg_tp,dpi_1ag1-dpi_agg_1ag1,0,1,W);
[sale_tfx,~] = my_fe(sale_tp,perc_tp,0,1,W);
[saleatfx,perc_tfx] = my_fe(sale_tp-perc_tp,perc_tp,0,1,W);

% Cross section and time fixed effects

[dpi_t_fxt,dpi_1_fxt] = my_fe(dpi_t_fx,dpi_1_fx,1,0,W);
[sale_tfxt,~] = my_fe(sale_tfx,perc_tfx,1,0,W);
[saleatfxt,perc_tfxt] = my_fe(saleatfx,perc_tfx,1,0,W);

% FGLS

O = ones(N,T);

[gls1,glsstdev1] = my_gls(sale_tp, [O dpi_tp]); % Benchmark case
[gls2,glsstdev2] = my_gls(sale_tft,dpi_t_ft); % Benchmark case including time fes
[gls3,glsstdev3] = my_gls((sale_tp-perc_tp),[O dpi_tp-dpi_agg_tp]); % Risk share          
[gls4,glsstdev4] = my_gls(sale_tft, dpi_1_ft); % Excess sensitive with time fes
[gls5,glsstdev5] = my_gls((sale_tp-perc_tp), [O (dpi_1ag1-dpi_agg_1ag1)]); % Excess sensitivity with flux
[gls6,glsstdev6] = my_gls(saleatfx,dpi_atfx); % Excess smoothnessw state fes
[gls7,glsstdev7] = my_gls(sale_tfxt,dpi_t_fxt); % Excess smoothnessw state and time fes

% Benchmark
fprintf('%0.4f \t %0.4f    \n', gls1(1,1), gls1(2,1))
fprintf('(%0.4f)  (%0.4f)    \n', glsstdev1(1,1), glsstdev1(2,1))

% Benchmark with time fes
fprintf('%0.4f    \n', gls2(1,1))
fprintf('(%0.4f)    \n', glsstdev2(1,1))

% Risk Share
fprintf('%0.4f \t %0.4f    \n', gls3(1,1), gls3(2,1))
fprintf('(%0.4f) (%0.4f)    \n', glsstdev3(1,1), glsstdev3(2,1))

% Excess Sensitive w time fes
fprintf('%0.4f    \n', gls4(1,1))
fprintf('(%0.4f)    \n', glsstdev4(1,1))

% Excess Sensitive w agg flux
fprintf('%0.4f \t %0.4f    \n', gls5(1,1), gls5(2,1))
fprintf('(%0.4f)  (%0.4f)    \n', glsstdev5(1,1), glsstdev5(2,1))

% Excess Smooth w State Fes
fprintf('%0.4f    \n', gls6(1,1))
fprintf('(%0.4f)    \n', glsstdev6(1,1))

% Excess Smooth w state and times fes
fprintf('%0.4f    \n', gls7(1,1))
fprintf('(%0.4f)    \n', glsstdev7(1,1))

%% Question 4

global x T MA

T = 300;  % Number of time periods.
MA = 1;% MA order.
beta0 = 0; % Beta 0.
beta1 = 0.5; % Beta 1.
if MA == 2
    beta2 = 0.5; % Beta 2.
end
sigma = 2; % Std

sim = 50;% Number of simulations
results_mat = zeros(sim,3); % Results matrix
if MA == 2
    results_mat = [results_mat zeros(sim,1)];
end

for s = 1:sim
   
    % Draw the error terms.
    u = normrnd(0,sigma,T,1);
    
    % Construct x recursively.
    x = zeros(T,1);
    
    if MA == 1 % For MA(1) process.        
        x(1) = u(1) + beta1*normrnd(0,sigma,1,1);
        for j = 2:T            
            x(j) = u(j) + beta1*u(j-1);
        end
        
        b0 = [0.1 0.5 0.75];
        
    elseif MA == 2 % For MA(2) process.        
        c = normrnd(0,sigma,1,1);
        x(1) = u(1) + beta1*c + beta2*normrnd(0,sigma,1,1);
        x(2) = u(2) + beta1*u(1) + beta2*c;
        
        for j = 3:T            
            x(j) = u(j) + beta1*u(j-1) + beta2*u(j-2);
        end
        
        b0 = [0.5 0.5 0.5 0.5];
        
    end
    
    % MLE.
    options = optimset('Display','off'); 
    [b_mle,~,~,~,~,hess] = fminunc(@logl_MA,b0,options);
    hess = inv(hess);
    results_mat(s,:) = b_mle';
end


fprintf(' %0.4f\n(%0.4f)\n',b_mle(1),sqrt(hess(1,1)))

fprintf(' %0.4f\n(%0.4f)\n',b_mle(2),sqrt(hess(2,2)))

if MA == 2
    fprintf(' %0.4f\n(%0.4f)\n',b_mle(4),sqrt(hess(4,4)))
else
    fprintf('None.\n')
end

fprintf(' %0.4f\n(%0.4f)\n',b_mle(3),sqrt(hess(3,3)))

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,1)), std(results_mat(:,1)))

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,2)), std(results_mat(:,2)))

if MA == 2
    fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,4)), std(results_mat(:,4)))
else
    fprintf('None.\n')
end

fprintf(' %0.4f\n(%0.4f)\n', mean(results_mat(:,3)), std(results_mat(:,3)))

bins = 10;

close all
figure(1)
hold on
histogram(results_mat(:,1),bins)
xlabel('$\beta_{0}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{0}$','interpreter','LaTex')
hold off



figure(2)
hold on
histogram(results_mat(:,2),bins)
xlabel('$\beta_{1}$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\beta_{1}$','interpreter','LaTex')
hold off


if MA == 2    
    figure(3)
    hold on
    histogram(results_mat(:,4),bins)
    xlabel('$\sigma$','interpreter','LaTex'); ylabel('Frequency')
    title('Plot of $\beta_2$','interpreter','LaTex')
    hold off    
else    
    fprintf('Not applicable.\n')
end


figure(4)
hold on
histogram(results_mat(:,3),bins)
xlabel('$\sigma$','interpreter','LaTex'); ylabel('Frequency')
title('Plot of $\sigma$','interpreter','LaTex')
hold off

%% Functions

function L = logl_AR(b)

% The following constructs the loglikelihood function for an MA(1).

global y T AR

b1 = b(1);                                                                  % Intercept.
b2 = b(2);                                                                  % AR coefficient.
s  = b(3);                                                                  % Standard deviation.

L = - 0.5*log(s^2/(1-b2^2)) - 0.5*((y(1) - b1/(1-b2))^2/(s^2/(1-b2^2)));

for t = 2:T
    L = L - 0.5*log(s^2) - 0.5*(((y(t)-b1-b2*y(t-1))^2/s^2));
end

L = L - 0.5*T*log(2*pi);
L = -L;                                                                     % Negative of loglikelihood function (for minimization).

end

function [ L ] = logl_MA(b)

% Loglikelihood for MA(1) and MA(2).

global x T MA

b0 = b(1);                                                                  % Intercept.
b1 = b(2);                                                                  % MA(1) coefficient.
s2 = b(3)^2;                                                                % Standard deviation.

if MA == 1                                                                  % For MA(1) process.

    omega = eye(T)*s2*(1+b1^2);                                             % Variances.
    cov1  = ones(T-1,1)*s2*b1;                                              % First order covariances.
    omega = omega + diag(cov1,1) + diag(cov1,-1);
    
elseif MA == 2                                                              % For MA(2) process.

    b2 = b(4);                                                              % MA(2) coefficient.

    omega = eye(T)*s2*(1+(b1^2)+(b2^2));                                    % Variances.
    cov1  = ones(T-1,1)*s2*(b1+(b1*b2));                                    % First order covariances.
    cov2  = ones(T-2,1)*s2*b2;                                              % Second order covariances.

    omega = omega + diag(cov1,1) + diag(cov1,-1);
    omega = omega + diag(cov2,2) + diag(cov2,-2);

end 
    
L = -0.5*log(det(omega)) - 0.5*((x-b0)'/omega)*(x-b0);                      % Loglikelihood function.
L = L - 0.5*T*log(2*pi);
L = -L;                                                                     % Negative of loglikelihood function (for minimization).

end

%% Funcations are copied in from Bent's provided code

function [diff_data] = my_diff(data,k,T)

% This code takes non-overlapping kth differences. The data set has T periods (columns).

if k == 1
    diff_data = data(:,2:end)-data(:,1:end-1);
else
    T_diff = floor((T-1)/k);
    diff_data = zeros(size(data,1),T_diff);
    
    for i = 1:T_diff
        diff_data(:,T_diff-i+1) = data(:,T-(i-1)*(k+1)) - data(:,T-(i-1)*(k+1)-k);
    end
end

end



function[Y,X] = my_fe(Y,X,tfe,cfe,W)

% This code removes the fixed effects (time or cross section; tfe and cfe, respectively) with wighting matrix, W, where W'W is the covariance matrix.

N = size(Y,1);
T = size(Y,2);

if tfe == 1
    % Time fixed effects.
    
    fe_Y = (sum(Y.*W))./(sum(W.^2));
    Y = Y - (W.*(ones(N,1)*fe_Y));
    
    fe_X = (sum(X.*W))./(sum(W.^2));
    X = X - (W.*(ones(N,1)*fe_X));
end

if cfe == 1 
    % Cross section fixed effects.
    
    fe_Y = (sum((Y.*W),2))./(sum((W).^2,2));
    Y = Y - (W.*(fe_Y*ones(1,T)));   
   
    fe_X = (sum((X.*W),2))./(sum((W).^2,2));
    X = X - (W.*(fe_X*ones(1,T)));
end

end

function [bfgls,se] = my_gls(Y,X)

% Panel estimation.

N = size(Y,1);                                                              % Number of states.
T = size(Y,2);                                                              % Number of years.
nobs = N*T;                                                                 % Sample size.

numx = size(X,2)/T;                                                         % Number of regressors.

% Step 1: OLS.

yvec = reshape(Y,nobs,1);
xvec = reshape(X,nobs,[]);

bols = (xvec'*xvec)\xvec'*yvec;
umat = reshape(yvec - xvec*bols,N,T);

% Estimate of variance matrix (with ridge).

ridge = 0.1;

varn = mean(umat.^2,2);
gmat = (1./T).*((umat*umat')./((varn.^(1/2))*(varn.^(1/2))'));                                       
gmat = ridge*eye(N) + (1-ridge).*gmat;
gmat = inv(gmat);

% Step 2: FGLS.

weights = varn.^(-1/2);
weights = kron(weights,ones(1,T));
weights = reshape(weights,[],1);

yhet = ((weights.*yvec)'*kron(eye(T),gmat))';                               % Normalize so residuals have variance 1 and weight by variance matrix.
xhet = ((repmat(weights,1,size(xvec,2)).*xvec)'*kron(eye(T),gmat))';

bfgls = (xhet'*xhet)\xhet'*yhet;

resid_het = yhet - xhet*bfgls;
s2 = (resid_het'*resid_het)/(nobs-numx);

vcmat = inv(xhet'*xhet)./s2;

se = (diag(vcmat)).^(1/2);

end

