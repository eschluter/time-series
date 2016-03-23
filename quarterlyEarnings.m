% Erik Schluter
% 7/15/2015
% Shumway and Stoffer - 2.1

function quarterlyEarnings(shouldAddIntercept)

    if (nargin < 1)
        shouldAddIntercept = false;
    end
    
    Y  = load('datasets\jj_series.txt');
    X  = log(Y);
    n  = length(Y);
    t  = 1:n;
    
    if (shouldAddIntercept)
        I   = ones(1, n);
        si  = 2;
    else
        I   = [];
        si  = 1;
    end
    
    Q1 = zeros(1, n); Q1(1:4:n) = 1;
    Q2 = zeros(1, n); Q2(2:4:n) = 1;
    Q3 = zeros(1, n); Q3(3:4:n) = 1;
    Q4 = zeros(1, n); Q4(4:4:n) = 1;
    
    Z  = [ I; t; Q1; Q2; Q3; Q4 ]';
    
    a  = (Z'*Z)\(Z'*X);
    
    intcpt  = (shouldAddIntercept)*a(1);
    beta    = a(si);
    alpha   = a((si+1):end);
    
    xt      = intcpt + beta*t + alpha(1)*Q1 + ...
                                alpha(2)*Q2 + ...
                                alpha(3)*Q3 + ...
                                alpha(4)*Q4;
                            
    % residuals
    figure(1)
    hold on
    plot(t./4, X - xt')
    hold off
    title('residuals')
    
    % model fit and data
    figure(2)
    hold on
    plot(t./4, xt, 'r')
    plot(t./4, X, 'k')
    hold off
    title('model fit and data')
    xlabel('year')
    ylabel('earnings')
    
return;