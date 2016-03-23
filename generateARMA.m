function [xt] = generateARMA(obs, phi, theta)
%GENERATEARMA generate ARMA(p,q) time series
%
%   [xt] = generateARMA(obs, phi, theta)
%
%   takes three arguments:
%       obs   - number of time points (observations) to generate
%       phi   - array of autoregressive (AR) coefficients
%       theta - array of moving average (MA) coefficients
%
%   Erik Schluter - 8/2/2015
    
    np = length(phi);
    nq = length(theta);
    xt = zeros(1, obs);
    wt = randn(1, obs);
    
    start = max(np,nq)+1;
    for (i = start:obs)
        % compute AR term
        ar = 0;
        for (p = 1:np)
            ar = ar + phi(p)*xt(i-p);
        end
        % compute MA term
        ma = 0;
        for (q = 1:nq)
            ma = ma + theta(q)*wt(i-q);
        end
        % put it all together
        xt(i) = ar + wt(i) + ma;
    end
    
return