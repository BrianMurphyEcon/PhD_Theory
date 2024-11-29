function W = weight(b,mat)

%{
    weight.m

    This code estimates the optimal weighting matrix. It implements
    pre-whitening using a VAR(1) model, calculates the optimal bandwidth,
    estimates the variance matrix according to the user-selected scheme
    then recolors to obtain the optimal weighting matrix.
%}

global lag n T

    if mat == 0                                                             % Identity weighting matrix.

        W = eye(n);

    elseif mat == 1 || mat == 2                                             % 1 = Newey-West, 2 = Quadratic Spectral.  
        %% Pre-whiten with VAR(1).

        f = orth(b);                                                        % Errors.
        
        X = f(1:T-lag-1,:);
        Y = f(2:T-lag,:);
        
        Ahat = (X'*X)\X'*Y;                                                 % VAR(1) coefficients.
        fstar = Y - X*Ahat;                                                 % Residuals.
        
        %% Calculate optimal bandwidth.
        
        X = f(1:T-lag-2,:);                                                 
        Y = f(2:T-lag-1,:);                                                 % Note we remove the most recent observation. See Bent's GMM3 notes.
        
        rho = diag(diag(X'*X))\diag(diag(X'*Y));                            % AR(1) coefficients.
        res = Y - X*rho;
        sigma2 = diag(res'*res)./(T-lag-2);
        
        rho = diag(rho);
        denom = sum((sigma2.^2)./((1-rho).^4));
        
        a1 = 4*((rho.^2).*(sigma2.^2))./(((1-rho).^6).*((1+rho).^2));        
        a1 = sum(a1)/denom;                                                 % alpha(1) as in Andrews.
        
        a2 = 4*((rho.^2).*(sigma2.^2))./((1-rho).^8);
        a2 = sum(a2)/denom;                                                 % alpha(2) as in Andrews.
        
        maxnw = round(10*(T-lag)^(1/3));
        
        K1 = min(round(1.1447*(a1*(T-lag-1))^(1/3)),maxnw);                 % Optimal bandwidth for Newey-West.
        K2 = min(round(1.3221*(a2*(T-lag-1))^(1/5)),25);                    % Optimal bandwidth for Quadratic Spectral.

        %% Estimate variance matrix.
        
        %omega = fstar'*fstar;                                              % For pre-whitening.
        omega = f'*f;                                                       % No pre-whitening.
        
        if mat == 1                                                         % Newey-West.
            
            for i = 1:K1                
                omegai = f(i+1:T-lag-1,:)'*f(1:T-lag-1-i,:);
                omega  = omega + (1-(i/K1))*(omegai+omegai');
            end
            
            omega = omega./(T-lag-1);
            
        elseif mat == 2                                                     % Quadratic Spectral.
            
            for i = 1:3*K2
                x = i/K2;
                kx = 25*((sin(6*pi*x/5)/(6*pi*x/5))...
                                            - cos(6*pi*x/5))/(12*pi^2*x^2);
                omegai = f(i+1:T-lag,:)'*f(1:T-lag-i,:);
                omega = omega + kx*(omegai+omegai');
            end
            
            omega = omega/(T-lag);
            
        end
        
        %% Recolor and obtain optimal weighting matrix.
        
        %omega = (eye(n)-Ahat)'\omega/(eye(n)-Ahat);                         % Recolor.
        
        W = inv(omega);                                                     % Optimal weighting matrix.

    end

end

