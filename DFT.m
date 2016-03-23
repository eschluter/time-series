function [ft] = DFT(data)
%DFT - compute the discrete fourier transform using covariance
%
%   [ft] = DFT(data)
%
%   takes one argument
%       data - raw data
%
%   returns DFT calculated over first n-1 fundamental frequencies
%    this is a symmetric two-sided spectrum
%
%   Erik Schluter - 8/20/2015

    n = length(data);
    ft = zeros(1,n);
    
    % get the sample autocovariance
    acov = autocov(data);
    
    % Set up the lag vector
    lags = linspace(-(n-1), n-1, 2*n-1);
    
    % calculate the power spectrum at the zeroth fund frequency
    ft(1) = n*(sum(data)/n)^2;
    
    % Calculate the periodogram
    for (j = 1:(n-1))
        for (h = 1:(2*n-1))
            ft(j+1) = ft(j+1) + acov(abs(lags(h))+1)*exp(-2*pi*1i*(j/n)*lags(h));
        end
    end
    
return