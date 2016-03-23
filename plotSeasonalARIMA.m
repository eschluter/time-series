% Erik Schluter
% 8/9/2015
% Shumway and Stoffer - 3.36

function plotSeasonalARIMA(obs)

    if (nargin < 1)
        obs = 200;
    end
    
    xt = generateSeasonalARIMA(obs, 0, 0.5, 0.8, 0, 12);
    
    figure(1)
    subplot(2,1,1); plot(ACF(xt), 'k')
    title('ACF of seasonal ARIMA model')   
    subplot(2,1,2); plot(PACF(xt), 'k')
    title('PACF of seasonal ARIMA model')
    
    figure(3)
    plot(xt, 'k')
    title('seasonal ARIMA model s = 12')
    
return