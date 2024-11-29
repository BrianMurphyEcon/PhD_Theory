%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Question 3, Part A.   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Import and Name Data

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
Alpha=PS6Data{:,23};
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
theta = PS6Data{:,52};



n=size(PS6Data,1);

% Continue Creating Variables Needed (was getting hard on Excel)

% Growth Rates

agg_RGNP=zeros(floor(n/4)-1,1);

for i=1:floor(n/4)-1
    agg_RGNP(i) = ((RGNP_WP(4*i+1)+RGNP_WP(4*i+2)+RGNP_WP(4*i+3)+RGNP_WP(4*i+4))-(RGNP_WP(4*i-3)+RGNP_WP(4*i-2)+RGNP_WP(4*i-1)+RGNP_WP(4*i)))./(RGNP_WP(4*i-3)+RGNP_WP(4*i-2)+RGNP_WP(4*i-1)+RGNP_WP(4*i));
end


mean(agg_RGNP)
mean(gr_RGNP)


% Part iii)


% Plot Results
plot(date,K_div_Y,date,X_div_Y,date,C_div_Y,date,R_div_Y)
legend('K/Y','X/Y','C/Y','R/Y')

% Calculate the Mean Var
mean_var=[mean(K_div_Y), mean(X_div_Y), mean(C_div_Y), mean(R_div_Y)];


% Part iv)
eta=zeros(n-1,1);
for i=1:n-1
    eta(i)=4*(Working_Pop(i+1)-Working_Pop(i))./Working_Pop(i);
end
eta=mean(eta,'all');
lambda=mean(agg_RGNP);
theta=(RK+r_cd)./(GNP+r_cd);
theta=mean(theta);
X_div_K=X./K;
delta=X_div_K+(1-(1+lambda)*(1+eta))*ones(n,1);
delta=mean(delta);
beta=(1+lambda)*(theta*(mean(K_div_Y))^(-1)-delta+1)^(-1);
h=0.31;
pi=(1-theta)*mean((C_div_Y).^-1)*((1-h)/h);

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
% Parameters from the calibrations
parameters=[beta, eta, pi, lambda, theta, delta, rho, eps];





