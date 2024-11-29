% Simulation - PS6
% Following PS4: 
%clear;

% After calibrating, we have our parameters.  They are reported below:
% Step 1: Input the parameter values

% Note, run the calibration code first to have the data imported into
% matlab.
%alpha is theta

alpha = 0.32324884926337;
beta = 0.961593698946701;	
eta = 0.00295488450935293;	
phi = 2.04523016146548;	
lambda = 0.0180805407126184;	
theta = 0.323335095140840;	
delta = 0.0251147610497855;	
rho = 0.988678542753938;	
eps = 0.00453508542486091;

t=60;

% -- Set the no. of periods of the impulse response to 60.

% Step 2: Compute the steady state 
% -- Steady-state values: y, k, c, i, kdivn, n
% Note we can solve without f-solve, since some variables are only
% functions of parameters

kdivn_ss=(1/theta*((1+lambda)/beta+delta-1))^(1/(theta-1));
n_ss= (1-theta)/(phi*(1-delta-lambda-eta+lambda*eta)*kdivn_ss^(1-theta)+1-theta);
k_ss= kdivn_ss*n_ss;
y_ss =(k_ss^theta)*n_ss^(1-theta);
c_ss =k_ss^theta*n_ss^(1-theta)+(1-delta)*k_ss-(1+lambda)*(1+eta)*k_ss;
i_ss= (delta+eta)*k_ss;

% Step 3: Linearizing the first-order conditions
% -- calculation of first and second derivatives of u & f (fill in the derivatives)

u1 =1/c_ss;
u11= -1/c_ss^2;
u2 =-phi/(1-n_ss);
u22=-phi/(1-n_ss)^2;
f1 =theta*k_ss^(theta-1)*n_ss^(1-theta);
f11=theta*(theta-1)*k_ss^(theta-2)*n_ss^(1-theta);
f13=theta*k_ss^(theta-1)*n_ss^(1-theta);
f2 =(1-theta)*k_ss^theta*n_ss^-theta;
f22=(1-theta)*(-theta)*k_ss^theta*n_ss^(-theta-1);
f23=(1-theta)*k_ss^theta*n_ss^-theta;
f12=theta*(1-theta)*k_ss^(theta-1)*n_ss^(-theta);
f21=f12;
f3 =(k_ss^theta)*(n_ss^(1-theta));

% -- fill in the coefficients in the linearized labor equation (4.14)

om =-u22-(f2^2)*u11-u1*f22; % omega in the notes
b1 =1/om*((f2*u11*(1+lambda+eta))/beta+u1*f21);
b2 =1/om*(f2*f3*u11+u1*f23);
b3 =-1/om*f2*u11*(1+lambda+eta);

% -- fill in the coefficients in the linearized Euler equation (4.15)

d1 =k_ss*((1+lambda+eta)/beta+f2*b1);
d2 =k_ss*(f2*b3-(beta*u1*f11/u11+(1+lambda+eta)*((1+lambda+eta)/beta)-(f2*(1+lambda+eta)+beta*u1*f11/u11)*b1));
d3 =k_ss*(-(f2*(1+lambda+eta)+beta*(u1/u11)*f12)*b3+1);
d4 =f2*b2+f3-rho*(f2*(1+lambda+eta)+beta*(u1/u11)*f12)*b2-rho*(f3*(1+lambda+eta)+beta*(u1/u11)*f13);

% Step 4: Compute the undetermined coefficients in (4.17)

% -- Solve the quadratic equation in pi_1.
p = [d3 d2 d1];
r = roots(p);

% -- Using pi_1, compute pi_2.
pi_1 = r(abs(r)<1);
pi_2 = -d4/(d2+pi_1*d3+d3*rho);

% Step 5: Generate a sequence of shocks for the impulse response

e= normrnd(0,eps,[t,1]);

z_hat= zeros(t,1);
k_hat= zeros(t,1);
n_hat =zeros(t,1);
y_hat=zeros(t,1);
i_hat =zeros(t,1);
c_hat= zeros(t,1);

z_hat(1) = e(1);

for t=2:t
    z_hat(t)= rho*z_hat(t-1)+e(t);
end

k_hat(1)= 0;
for t=2:t
    k_hat(t)= pi_1*k_hat(t-1)+pi_2*z_hat(t-1);
end

n_hat(1)=0;
for t=1:t
    n_hat(t)= kdivn_ss*(b1+b3*pi_1)*k_hat(t)+(b2/n_ss+kdivn_ss*b3*pi_2)*z_hat(t);
end

y_hat(1)=0;
for t=1:t
    y_hat(t)= z_hat(t)+theta*k_hat(t)+(1-theta)*n_hat(t);
end

i_hat(1)=0;
for t=1:t
    i_hat(t) =(pi_1/delta-((1-delta)/delta))*k_hat(t)+pi_2/delta*z_hat(t);
end

c_hat(1)=0;
for t=1:t
    c_hat(t) =y_ss/c_ss*y_hat(t)-i_ss/c_ss*i_hat(t);
end


% Step 7: Plot the impulse respones.

% Create a figure with subplots
figure;

% Plot Capital Response
subplot(3,2,1);
plot(k_hat);
title('Capital Response');

% Plot Output Response
subplot(3,2,2);
plot(y_hat);
title('Output Response');

% Plot Consumption Response
subplot(3,2,3);
plot(c_hat);
title('Consumption Response');

% Plot Labor Response
subplot(3,2,4);
plot(n_hat);
title('Labor Response');

% Plot Investment Response
subplot(3,2,5);
plot(i_hat);
title('Investment Response');



% HP Filter to Detrend the Data. HP gives a trend and a cyclical variable.
% Name them as such

[trend_K,cyclical_K]= hpfilter(K); % Capital
[trend_i,cyclical_i]=hpfilter(X); % Invest
[trend_N,cyclical_N]= hpfilter(g_Avg_Hrs_Worked); % Labor
[trend_C,cyclical_C]= hpfilter(PCEC); % Note; PCEC is consumption
[trend_Y,cyclical_Y]= hpfilter(Y); % Output

% Simulating data for 225 periods. Create new time variable by naming T

T = 225; 

e=normrnd(0,eps,[T,1]);
z_hat = zeros(T,1);
k_hat = zeros(T,1);
c_hat = zeros(T,1);
n_hat = zeros(T,1);
i_hat = zeros(T,1);
y_hat = zeros(T,1);

for t=2:T
    z_hat(t)=rho*z_hat(t-1)+e(t);
end

% Estimate the 'hat' variables

for i = 2:T
    k_hat(i+1)= pi_1*k_hat(i)+pi_2*z_hat(i);
end 

for i = 2:T
    n_hat(i)= kdivn_ss*(b1+b3*pi_1)*k_hat(i)+(b2/n_ss+kdivn_ss*b3*pi_2)*z_hat(i);
end

for i = 2:T
    y_hat(i)= z_hat(i)+theta*k_hat(i)+(1-theta)*n_hat(i);
end

for i = 2:T
    i_hat(i) =(pi_1/delta-((1-delta)/delta))*k_hat(i)+pi_2/delta*z_hat(i);
end

for i = 2:T
    c_hat(i)= y_ss/c_ss*y_hat(i)-i_ss/c_ss*i_hat(i);
end

% Next, drop the first 25 observations
k_hat= k_hat(25:end);
y_hat = y_hat(25:end);
c_hat = c_hat(25:end);
i_hat = i_hat(25:end);
n_hat = n_hat(25:end);
z_hat = z_hat(25:end);

% Now, apply the HP Filter
[trend_k_hat,cyclical_k_hat]= hpfilter(k_hat);
[trend_y_hat,cyclical_y_hat] =hpfilter(y_hat);
[trend_c_hat,cyclical_c_hat]= hpfilter(c_hat);
[trend_n_hat,cyclical_n_hat]= hpfilter(n_hat);
[trend_i_hat,cyclical_i_hat] =hpfilter(i_hat);

% Then calculate the standard deviation divided by y of both hat and
% non-hat

standard_dev_K= std(trend_K)/std(trend_Y);
standard_dev_i =std(trend_i)/std(trend_Y);
standard_dev_C= std(trend_C)/std(trend_Y);
standard_dev_n= std(trend_N)/std(trend_Y);
sd_Y=std(trend_Y);

standard_dev_K_hat= std(trend_k_hat)/std(trend_y_hat);
standard_dev_c_hat=std(trend_c_hat)/std(trend_y_hat);
standard_dev_i_hat= std(trend_i_hat)/std(trend_y_hat);
standard_dev_y_hat= std(trend_y_hat)/std(trend_y_hat);
standard_dev_n_hat= std(trend_n_hat)/std(trend_y_hat);


% Create final time variable used to show time periods not yet removed.

Time = 200;

% Finally, calculate the auto correlation for hat and non-hat

% For k
k_hat_1 =k_hat(2:Time,:);
k_hat_2= k_hat(1:Time-1,:);

corr_k_hat= corrcoef(k_hat_1, k_hat_2);
ac_k_hatv =corr_k_hat(1,2);

% For c
c_hat_1=c_hat(2:Time,:);
c_hat_2= c_hat(1:Time-1,:);

corr_c_hat =corrcoef(c_hat_1, c_hat_2);
ac_c_hat=corr_c_hat(1,2);

% For i
i_hat_1=i_hat(2:Time,:);
i_hat_2= i_hat(1:Time-1,:);

corr_i_hat =corrcoef(i_hat_1, i_hat_2);
ac_i_hat=corr_i_hat(1,2);

% For y
y_hat_1= y_hat(2:Time,:);
y_hat_2= y_hat(1:Time-1,:);

corr_y_hat= corrcoef(y_hat_1, y_hat_2);
ac_y_hat=corr_y_hat(1,2);

% For n
n_hat_1=n_hat(2:Time,:);
n_hat_2=n_hat(1:Time-1,:);

corr_n_hat=corrcoef(n_hat_1, n_hat_2);
ac_n_hat=corr_n_hat(1,2);

% For k
cyclical_K_1=cyclical_K(2:n,:);
cyclical_K_2=cyclical_K(1:n-1,:);

corr_cyclical_K=corrcoef(cyclical_K_1, cyclical_K_2);
ac_cyclical_K=corr_cyclical_K(1,2);

% For y
cyclical_Y_1=cyclical_Y(2:n,:);
cyclical_Y_2=cyclical_Y(1:n-1,:);

corr_cyclical_Y =corrcoef(cyclical_Y_1, cyclical_Y_2);
ac_cyclical_Y= corr_cyclical_Y(1,2);

% For c
cyclical_C_1= cyclical_C(2:n,:);
cyclical_C_2=cyclical_C(1:n-1,:);

corr_cyclical_C =corrcoef(cyclical_C_1, cyclical_C_2);
ac_cyclical_C=corr_cyclical_C(1,2);

% For i
cyclical_i_1 =cyclical_i(2:n,:);
cyclical_i_2=cyclical_i(1:n-1,:);

corr_cyclical_i= corrcoef(cyclical_i_1, cyclical_i_2);
ac_cyclical_i=corr_cyclical_i(1,2);

