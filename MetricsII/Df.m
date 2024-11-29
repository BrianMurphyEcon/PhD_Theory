function grad = Df(b)

%{
    Df.m

    This computes the gradient of the moment conditions. It take the
    estimates of beta and gamma, respectively, as inputs.
%}

global lag T

Dgp1 = sum(orth([b(1)+eps b(2)    ]),1)./(T-lag);
Dgm1 = sum(orth([b(1)-eps b(2)    ]),1)./(T-lag);
Dgp2 = sum(orth([b(1)     b(2)+eps]),1)./(T-lag);
Dgm2 = sum(orth([b(1)     b(2)-eps]),1)./(T-lag);

grad(:,1) = (Dgp1-Dgm1)'./(2*eps);                                          % Derviative wrt beta.
grad(:,2) = (Dgp2-Dgm2)'./(2*eps);                                          % Derivative wrt gamma.

end

