% Erik Schluter
% 8/9/2015
% Shumway and Stoffer - 3.10

function [phi] = cmortAR(order, horizon, suppress)

    if (nargin < 1)
        order = 2;
    end
    
    if (nargin < 2)
        horizon = 4;
    end
    
    if (nargin < 3)
        suppress = false;
    end

    M       = load('datasets\cmort_series.txt');
    Mtrunc  = M((order+1):end);
    n       = length(M);
    X       = zeros(order, n-order);
    
    % Build lagged data set
    for i = 1:order
        X(i,:) = M((order-i+1):(end-i))';
    end
    
    % Regression
    Z       = [ones(1,n-order); X]';
    
    phi     = (Z'*Z)\(Z'*Mtrunc);
    
    %calculate fit
    xt      = phi(1);
    for (j = 1:order)
        xt  = xt + X(j,:)*phi(j+1);
    end
    
    if (~suppress)
        figure(1)
        hold on
        plot(M, 'k');
        plot(xt, 'r');
        title(['AR(' num2str(order) ') fit using regression'])

        % Forecast
        for (p = 1:horizon)
            f_est = phi(1);
            for (j = 1:order)
                f_est = f_est + M(end-j)*phi(j+1);
            end
            M  = [M; f_est];
        end
        plot(n:(n+horizon), M(n:end), 'g')
        hold off
    end
    
return