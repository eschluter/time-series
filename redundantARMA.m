% Erik Schluter
% 8/7/2015
% Shumway and Stoffer - 3.20

function redundantARMA(obs)

    if (nargin < 1)
        obs = 500;
    end
    
    phi     = 0.9;
    theta   = -0.9;
    
    xt = generateARMA(obs, phi, theta);
    
    figure(1)
    plot(xt, 'k')
    title('Redundant ARMA(1,1) process')
    
    % plot ACF and PACF
    figure(2)
    subplot(2,1,1); plot(ACF(xt), 'k'); title('ACF')
    subplot(2,1,2); plot(PACF(xt), 'k'); title('PACF')
    
    b = fitARMAgn(xt, phi, theta);
    
    disp(['estimated coeffs: ' mat2str(b)]);
    
return