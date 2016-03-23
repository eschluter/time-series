% Erik Schluter
% 7/15/2015
% Shumway and Stoffer - 2.8

function glacialVarves()

    xt          = load('datasets\varve.txt');
    n           = length(xt);
    
    mid         = floor(n/2);
    xt_h1       = xt(1:mid);
    xt_h2       = xt((mid+1):end);
    
    [xtVar1,~]  = sampleVar(xt_h1, length(xt_h1));
    [xtVar2,~]  = sampleVar(xt_h2, length(xt_h2));
    xratio      = xtVar1/xtVar2;
    
    disp(['First half data variance = ' num2str(xtVar1)])
    disp(['Second half data variance = ' num2str(xtVar2)])
    disp(['Xt variance ratio = ' num2str(xratio)])
    
    figure(1)
    hx = hist(xt,50);
    bar(hx);
    title('raw data distribution')
    
    % log tranformation
    yt          = log(xt);   
    yt_h1       = yt(1:mid);
    yt_h2       = yt((mid+1):end);
    
    [ytVar1,~]  = sampleVar(yt_h1, length(yt_h1));
    [ytVar2,~]  = sampleVar(yt_h2, length(yt_h2));
    [ytVar, ym] = sampleVar(yt, n);
    yratio      = ytVar1/ytVar2;
    
    disp(['First half data variance = ' num2str(ytVar1)])
    disp(['Second half data variance = ' num2str(ytVar2)])
    disp(['Yt variance ratio = ' num2str(yratio)])
    
    figure(2)
    hy = hist(yt,50);
    bar(hy);
    title('Log transformed data distribution')
    
    figure(3)
    hold on
    plot(xt, 'k')
    plot(yt, 'r')
    hold off
    title('overlayed Log transform and raw varve data')
    
    figure(4)
    plot(yt, 'k')
    title('Log transform of varve data')
    
    % compute the sample ACF
    acf_lag = zeros(1, n-1);
    for i = 1:(n-1)
        s = 0;
        for j = 1:(n-i)
            s = s + (yt(j) - ym)*(yt(j+i) - ym);
        end
        acf_lag(i) = s/(ytVar*(n-1));
    end
    
    % compare to built in matlab function just for kicks
    builtInACF = autocorr(yt, n-1);
    
    figure(5)
    hold on
    plot(2:n, acf_lag, 'o-k', 'markers', 2)
    plot(builtInACF, 'r')
    hold off
    title('Sample ACF')
    
    % plot differenced log data and sample acf
    ut      = diff(yt);
    utACF   = autocorr(ut, n-2);
    figure(6)
    plot(ut, 'k')
    title('Differenced log data')
    figure(7)
    plot(utACF, 'b')
    title('Differenced log data ACF')
    
return

function [v, m] = sampleVar(data, n)
    m = sum(data)/n;
    v = sum((data - m).^2)/(n - 1);
return