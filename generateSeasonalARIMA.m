function [xt] = generateSeasonalARIMA(obs, phi, theta, s_phi, s_theta, s)
%GENERATESEASONALARIMA generate ARIMA(p,q)x(P,Q)s time series
%
%   [xt] = generateARMA(obs, phi, theta)
%
%   takes six arguments:
%       obs     - number of time points (observations) to generate
%       phi     - array of autoregressive (AR) coefficients
%       theta   - array of moving average (MA) coefficients
%       s_phi   - array of seasonal AR coefficients
%       s_theta - array of seasonal MA coefficients
%       s       - seasonal parameter
%
%   Erik Schluter - 8/2/2015
    
    np = length(phi);
    nq = length(theta);
    nP = length(s_phi);
    nQ = length(s_theta);
    xt = zeros(1, obs);
    wt = randn(1, obs);
    
    start = max([np nq s*nP s*nQ])+1;
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
        % compute sAR term
        s_ar = 0;
        for (P = 1:nP)
            s_ar = s_ar + s_phi(P)*xt(i-P*s);
        end
        % compute sMA term
        s_ma = 0;
        for (Q = 1:nQ)
            s_ma = s_ma + s_theta(Q)*wt(i-Q*s);
        end
        % put it all together
        xt(i) = ar + wt(i) + ma + s_ar + s_ma;
    end
    
return