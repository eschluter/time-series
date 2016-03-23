function [acf] = ACF(data)
%ACF compute the sample ACF
%
%   [acf] = ACF(data)
%
%   takes 1 input argument
%       data - raw data
%   
%   returns n-1 lags
%   The first value in the return array is at lag 1.
%   The zeroth lag is omitted due to one indexing
%
%   Erik Schluter - 8/5/2015

    n = length(data);
    
    %compute sample mean and variance of data
    m = sum(data)/n;
    v = sum((data - m).^2)/(n - 1);
    
    acf = zeros(1, n-1);
    for i = 1:(n-1)
        s = 0;
        for j = 1:(n-i)
            s = s + (data(j) - m)*(data(j+i) - m);
        end
        acf(i) = s/(v*(n-1));
    end
    
return