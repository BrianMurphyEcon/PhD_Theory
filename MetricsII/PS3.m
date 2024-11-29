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




