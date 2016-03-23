% Erik Schluter
% 7/15/2015
% Shumway and Stoffer - 2.3

function randomWalkDrift()

    n       = 100;
    t       = 1:n;
    delta   = 0.01;
    
    for i = 1:6
        x       = randn(n, 1);
        X       = delta*(t') + cumsum(x);

        Z  = t';

        b  = (Z'*Z)\(Z'*X);

        xt = b(1)*t;

        subplot(3, 2, i)
        hold on
        plot(t, xt, 'r')
        plot(t, 0.01*t, 'b')
        plot(t, X, 'k')
        hold off
    end
    
return