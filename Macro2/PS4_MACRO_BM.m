%=======================================
%
%   Guidelines to PS #4
%
%=======================================

clear all 
close all 

% Step 1: Input the parameter values
global alpha delta beta phi sigma rho T

%Param Values

alpha = 0.3;
delta = 0.02;
beta = 0.9901;
phi = 2.0417;
sigma = 0;
rho = 0.8;
T=40;

% Step 2: Compute the steady state 
x0 = [0.5, 0.75, 0.5]; % Initial guess for [k,n,c]

steady_state = fsolve(@funs, x0);

% Display results
k_ss = steady_state(1);
n_ss = steady_state(2);
c_ss = steady_state(3);

y_ss = c_ss-(delta*k_ss);
i_ss= delta*k_ss;
% Step 3: Linearizing the first-order conditions
% -- calculation of first and second derivatives of u & f (fill in the derivatives)

u1 = 1/c_ss;
u11 = -1/c_ss^2;
u2 = -phi/1-n_ss;
u22 = -phi/(1-n_ss)^2;
f1 = alpha*(k_ss^(alpha-1))*n_ss^(1-alpha);
f11= (alpha)*(alpha-1)*(k_ss^(alpha-2))*(n_ss^(1-alpha));
f12= alpha*(1-alpha)*(k_ss^(alpha-1))*(n_ss^(-alpha));
f13= alpha*(k_ss^(alpha-1))*(n_ss^(1-alpha));
f2 = (1-alpha)*(k_ss^(alpha))*(n_ss^(-alpha));
f21 = (1-alpha)*(alpha)*(k_ss^(alpha-1))*(n_ss^(-alpha));
f22= (1-alpha)*(-alpha)*(k_ss^(alpha))*(n_ss^(-alpha-1));
f23= (1-alpha)*(k_ss^(alpha))*(n_ss^(-alpha));
f3 = k_ss^(alpha)*n_ss^(1-alpha);

% -- fill in the coefficients in the linearized labor equation (4.14)
om = -u22-(f2^2)*u11-u1*f22 ; 
b1 = (1/om)*((f2*u11)*(1/beta)+u1*f21);
b2 = (1/om)*(f2*f3*u11 + u1*f23);
b3 = -(1/om)*f2*u11;

% -- fill in the coefficients in the linearized Euler equation (4.15)
d1 = k_ss*((1/beta)+f2*b1);
d2 = k_ss*(f2*b3-(1+(1/beta)+beta*(u1/u11)*f11)-(f2+beta*(u1/u11)*f12)*b1);
d3 = k_ss*(-(f2+beta*(u1/u11)*f12)*b3+1);
d4 = f2*b2+f3-rho*(f2+beta*(u1/u11)*f12)*b2-rho*(f3+beta*(u1/u11)*f13);


% Step 4: Compute the undetermined coefficients in (4.17)

% -- Solve the quadratic equation in pi_1.
p = [d3 d2 d1];
r = roots(p);

% -- Using pi_1, compute pi_2.
pi_1 = r(abs(r)<1);
pi_2 = -d4/(d2+pi_1*d3+d3*rho);

% Step 5: Generate a sequence of shocks for the impulse response
% -- Define a sequence z_hat.
z_hat = NaN*ones(T,1);
z_hat(1) = .01;
for t = 1:T-1 % Iterate up to T-1 to make zhat 39x1
    z_hat(t+1) = rho*(z_hat(t)) + sigma*randn(1) ;
end

% Step 6: Compute the impluse responses
% -- For k_hat, use (4.17).
k_hat=NaN*ones(T,1);
k_hat(1) = 0;
for T = 1:T-1
    k_hat(T+1) = pi_1*k_hat(T) + pi_2*z_hat(T);
end
% -- For other variables, see p.35.
n_hat=NaN*ones(T,1);
n_hat = (k_ss/n_ss)*(b1+b3*pi_1)*k_hat+((1/n_ss)*b2+(k_ss/n_ss)*b3*pi_2)*z_hat;
y_hat=NaN*ones(T,1);
y_hat = z_hat+alpha*k_hat+(1-alpha)*n_hat;
i_hat=NaN*ones(T,1);
i_hat = (pi_1/delta)*k_hat+(pi_2/delta)*z_hat-((1-delta)/delta)*k_hat;
c_hat=NaN*ones(T,1);
c_hat = (y_ss/c_ss)*y_hat-(i_ss/c_ss)*i_hat;

what = z_hat+alpha.*k_hat-alpha.*n_hat;
rhat = (1-beta*(1-delta))*(y_hat-[0;k_hat(1:end-1)]);

% Step 7: Plot the impulse respones.
figure
plot(k_hat)
title('Capital Response')

figure
plot(y_hat)
title('Output Response')

figure
plot(c_hat)
title('Consumption Response')

figure
plot(n_hat)
title('Labor Response')

figure
plot(i_hat)
title('Investment Response')

figure
plot(rhat)
title('R Response')

figure
plot(what)
title('W Response')


function F = funs(x)
global alpha delta beta phi
%x(1) = k_ss
%x(2) = n_ss
%x(3) = c_ss
    
    F(1) = (beta * (alpha * x(1)^(alpha - 1)*x(2)^(1-alpha)+(1-delta))) - 1;
    F(2) = (1/x(3))*x(1)^(alpha)*x(2)^(-alpha)*(1-alpha)+(-phi/(1-x(2)));
    F(3) =  x(1)^alpha * x(2)^(1-alpha)-x(3)-delta*x(1);
end
