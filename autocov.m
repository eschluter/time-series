function [acov] = autocov(data)
%AUTOCOV compute the sample autocovariance
%
%   [acov] = autocov(data)
%
%   takes 1 input argument
%       data - raw data
%   
%   returns n-1 lags
%   The first value in the return array is at lag 0
%   So lag 0 zero is at index 1
%
%   Erik Schluter - 8/17/2015

    n = length(data);
    
    %compute sample mean of data
    m = sum(data)/n;
    
    acov = zeros(1, n);
    for i = 0:(n-1)
        s = 0;
        for j = 1:(n-i)
            s = s + (data(j) - m)*(data(j+i) - m);
        end
        acov(i+1) = s/n;
    end
    
return