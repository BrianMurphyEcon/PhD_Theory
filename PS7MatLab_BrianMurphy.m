clear; % this clears variables from the workspace
clc; % this clears the command window

% Section 1: Parameters
beta = 0.96;
b = 75;
N = 1000;
T = 100;
toler = .0001;

% Create Random Wage
w_l = 50;
w_h = 150;
wt = linspace(w_l, w_h, N);


%Define Value Function

V = zeros(N,1); % Initial Value Function 
v = V; % This will become the 'new' value function

matV = NaN(N,T+1); % Matrix to store the value function
matV(:,1) = V;

w_p = []; 
r=[];

tic
for j=1:T	
    for i = 1: N
% Assuming log utility
    if (log(wt(i)) / (1 - beta)) >= (log(b) + beta * sum((1/N)*V))
         r = log(wt(i)) / (1 - beta);
         w_p(i)=1;
    else
         r = log(b) + beta * sum((1/N)*V);
         w_p(i)=0;
    end  
    [vscal] = max(r); % this is the actual maximization line
    v(i) = vscal; % store the values of the new value function
    end 
matV(:,j+1) = v; % Storing the new guess for the value function
    
    diff=(V-v)./V;
    
    if abs(diff) <= toler
        break
    else
    % so here if the diff exceeds the tolerance, then do a replacement
    % and go back to the top of the loop
        V=v;
    end

    %Calculate the reservation wage
    
    wstar = exp((log(b) + beta * sum((1/N)*V))*(1-beta));

end

toc;


% This plots the value function
close all
figure(1)
plot(wt,matV(:,j),'k-');
xlabel('$$ w_{t} $$','interpreter','LaTex')
ylabel('$$ v(w_{t}) $$','interpreter','LaTex')

% This plots the policy function
figure(2)
plot(wt, w_p','k-');
xlabel('$$ w_{t} $$','interpreter','LaTex')
ylabel('$$ choice $$','interpreter','LaTex')



    








 