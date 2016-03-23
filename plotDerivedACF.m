% Erik Schluter
% 8/28/2015
%
% This example compares my derivation for the autocovariance
%	function of an ARMA(2,2) process with the sample autocov

function plotDerivedACF(obs)

    if (nargin < 1)
        obs = 1000;
    end

    p1_1 = 0.5;     p1_2 = -0.5;
    p2_1 = 0.29;    p2_2 = 0.29;
    t1_1 = 0.21;    t1_2 = 0.21;
    t2_1 = 0.1;     t2_2 = 0.1;
    sw2 = 1;
    
    xt1 = generateARMA(obs, [p1_1 p2_1], [t1_1 t2_1]);
    xt2 = generateARMA(obs, [p1_2 p2_2], [t1_2 t2_2]);

    figure(1)
    subplot(2,1,1); plot(xt1, 'k')
    subplot(2,1,2); plot(xt2, 'k')

    figure(2)
    subplot(2,1,1)
    hold on
    plot(theoreticalACOV(obs, p1_1, p2_1, t1_1, t2_1, sw2), 'k')
    plot(autocov(xt1), 'o-r')
    title('Sample (red dots) and theoretical (black) autocovariance')
    hold off
    subplot(2,1,2)
    hold on
    plot(theoreticalACOV(obs, p1_2, p2_2, t1_2, t2_2, sw2), 'k')
    plot(autocov(xt2), 'o-r')
    hold off
    
    
return

function acov = theoreticalACOV(obs, p1, p2, t1, t2, sw2)

    % derived expression for zeroth lag initial condition
    numer = sw2*((p1*t1)/(1-p2) + (p1*p2*t1)/(1-p2) + p2*t2 + p1*t1 + t1^2 + p2*t2 + t2^2 + 1);
    denom = 1 - (p1^2)/(1-p2) - ((p1^2)*p2)/(1-p2) - p2^2;
    
    % initial conditions
    gamma0 = numer/denom;
    gamma1 = (p1/(1-p2))*gamma0 + (t1/(1-p2))*sw2;
    gamma2 = p1*gamma1 + p2*gamma0 + t2*sw2;
    
    acov    = zeros(1,obs);
    acov(1) = gamma0;
    acov(2) = gamma1;
    acov(3) = gamma2;
    for i = 4:obs
        acov(i) = p1*acov(i-1) + p2*acov(i-2);
    end
    
return