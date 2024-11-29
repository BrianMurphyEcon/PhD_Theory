function crit = gmm_obj(guess,W)

%{
    gmm_obj.m

    This is the quadratic GMM objective function for Hansen and Singleton
    (1982).
%}

global lag T

mom = ((sum(orth(guess),1))./(T-lag))';
crit = mom'*W*mom;

end

