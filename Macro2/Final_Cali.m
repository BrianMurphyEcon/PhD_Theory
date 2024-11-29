%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Question 3, Part A.   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear;

% Import and Name Data

% My data cuts off after Q4 of 2019. I could only find updated working pop
% data until 2019, nothing was updated after. I ran into problems
% calibrating the data with NaN's so I removed data until all variables 
% were completely filled out. The additional data is available on the
% "unused_data" tab on the excel file I turned in

date=PS6Data{:,1};
t=PS6Data{:,2};
GDP=PS6Data{:,3};
GNP=PS6Data{:,4};
NICUR=PS6Data{:,5};
COFC=PS6Data{:,6};
CProfit=PS6Data{:,7};
PropInc=PS6Data{:,8};
RentIn=PS6Data{:,9};
GDP_Def=PS6Data{:,10};
GNP_Def=PS6Data{:,11};
NNP=PS6Data{:,12};
Net_Interest=PS6Data{:,13};
RGNP=PS6Data{:,14};
Total_RGNP=PS6Data{:,15};
Population=PS6Data{:,16};
Total_Pop=PS6Data{:,17};
Working_Pop=PS6Data{:,18};
Total_Working_Pop=PS6Data{:,19};
WP_Growth=PS6Data{:,20};
RGNP_WP=PS6Data{:,21};
UI=PS6Data{:,22};
Tot_Alpha=PS6Data{:,23};
AI=PS6Data{:,24};
Fixed_Assets=PS6Data{:,25};
Depreciated_Assets=PS6Data{:,26};
Cons_Durable_Goods=PS6Data{:,27};
Cons_Dur_Goods_Depr=PS6Data{:,28};
PCEC=PS6Data{:,29};
Avg_Hrs_Worked= PS6Data{:,30};
Wages=PS6Data{:,32};
K=PS6Data{:,32};
K_Cons_Dur=PS6Data{:,33};
RK=PS6Data{:,34};
r=PS6Data{:,35};
I=PS6Data{:,36};
avg_I = PS6Data{:,37};
alpha_cd = PS6Data{:,38};
r_cd = PS6Data{:,39};
R = PS6Data{:,40};
Y = PS6Data{:,41};
X= PS6Data{:,42};
Total_K = PS6Data{:,43};
gr_RGNP = PS6Data{:,44};
K_div_Y= PS6Data{:,45};
X_div_Y=PS6Data{:,46};
C_div_Y=PS6Data{:,47};
R_div_Y=PS6Data{:,48};
Y_div_C = PS6Data{:,49};
g_Avg_Hrs_Worked = PS6Data{:,50};
Ldiv_n=PS6Data{:,51};
mean_K_div_Y = PS6Data{1,52};
mean_X_div_Y = PS6Data{1,53};
mean_C_div_Y =PS6Data{1,54};
mean_R_div_Y = PS6Data{1,55};
mean_gr_RGNP = PS6Data{1,56};
eta_all = PS6Data{:,57};
eta= PS6Data{1,58};
agg_RGNP = PS6Data{:,59};
lambda = PS6Data{1,60};
alpha = PS6Data{1,61};
X_div_K = PS6Data{:,62};
tot_delta = PS6Data{:,63};
delta = PS6Data{1,64};
h = PS6Data{1,65};
theta=(RK+r_cd)./(GNP+r_cd);
theta=mean(theta);
beta=(1+lambda)*(theta*(mean(K_div_Y))^(-1)-delta+1)^(-1);
pi=(1-theta)*mean((C_div_Y).^-1)*((1-h)/h);

n=size(PS6Data,1);

%agg_RGNP=zeros(floor(n/4)-1,1);

%for i=1:floor(n/4)-1
 %   agg_RGNP(i) = ((RGNP_WP(4*i+1)+RGNP_WP(4*i+2)+RGNP_WP(4*i+3)+RGNP_WP(4*i+4))-(RGNP_WP(4*i-3)+RGNP_WP(4*i-2)+RGNP_WP(4*i-1)+RGNP_WP(4*i)))./(RGNP_WP(4*i-3)+RGNP_WP(4*i-2)+RGNP_WP(4*i-1)+RGNP_WP(4*i));
%end

% I copied and pasted this back into excel to continue working there. I
% kept having problems with the above in excel.

% Plot Results
plot(date,K_div_Y,date,X_div_Y,date,C_div_Y,date,R_div_Y)
legend('K/Y','X/Y','C/Y','R/Y')


% Add in Shocks

% Calculate the Tech Shocks
K_diff=zeros(n,1);
Y_diff=zeros(n,1);
h_diff=zeros(n,1);

for i=2:n
    K_diff(i)=log(K(i))-log(K(i-1));
end
K_diff(1)=K_diff(2);

for i=2:n
    Y_diff(i)=log(Y(i))-log(Y(i-1));
end
Y_diff(1)=Y_diff(2);

for i=2:n
    h_diff(i)=log(g_Avg_Hrs_Worked(i))-log(g_Avg_Hrs_Worked(i-1));
end
h_diff(1)=h_diff(2);

plot (date,Y_diff,date,K_diff)

% Solow residual
z = log(Y)-theta*log(K)-(1-theta)*log(g_Avg_Hrs_Worked);

% Normalize
norm_z = z(1); 

for i = 1:length(z)
    z(i) = z(i)/norm_z;
end

% Getting the autocorrelation and eps
ln_z = log(z);

% AR(1) Process
autocorr = autocorr(ln_z, 1);
rho = autocorr(2);
sigma_z = var(ln_z);
eps = ones(length(z-1),1);

for i = 2:length(z)
    eps(i) = z(i) - rho*z(i-1);
end
% The estimator of the variance would be:
eps = (eps'*eps)/(n-1);
%% Parameters from the calibrations
parameters=[beta, eta, pi, lambda, theta, delta, rho, eps];















