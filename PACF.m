function [pacf] = PACF(data)
%PACF compute partial autocorrelation using Durbin-Levinson
%
%   [pacf] = PACF(data)
%
%   takes one argument
%       data - array containing raw data
%
%   Erik Schluter - 8/5/2015

    global phi;
    dn    = length(data);
    npacf = floor(dn/2);
    pacf  = zeros(1,npacf);
    phi   = zeros(npacf);
    
    % Get autocorrelation
    acf  = ACF(data);
    
    % set initial conditions
    pacf(1)  = acf(1);
    phi(1,1) = pacf(1);
    
    % compute PACF iteratively using Durbin-Levinson
    for (i = 2:npacf)
        d_numer = 0;
        d_denom = 0;
        for (j = 1:(i-1))
            d_numer = d_numer + phi(i-1,j)*acf(i-j);
            d_denom = d_denom + phi(i-1,j)*acf(j);
        end
        %calculate pacf(lag i)
        pacf(i)  = (acf(i) - d_numer) / (1 - d_denom);
        phi(i,i) = pacf(i);
        %calculate the remaining phi terms
        for (j = 1:(i-1))
            PHI(i,j);
        end
    end
    
return

function PHI(n,k)
    global phi;
    phi(n,k) = phi(n-1,k) - phi(n,n)*phi(n-1,n-k);
return