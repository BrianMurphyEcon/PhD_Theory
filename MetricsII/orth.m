function ZXb = orth(guess)

%{
    orth.m

    This is the orthogonality condition from Hansen and Singleton (1982)
    E(z(t)*((beta*((C(t)/C(t-1))^(-gamma))*r(t)) - 1)) = 0.
%}

global c lag n rf T Z

beta  = guess(1);
gamma = guess(2);

C = repmat(c(1+lag:T),1,n);
R = repmat(rf(1+lag:T),1,n);

ZXb = Z.*((beta.*(C.^(-gamma)).*R)-1);

end

