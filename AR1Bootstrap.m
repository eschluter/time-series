% Erik Schluter
% 8/10/2015
% Shumway and Stoffer - 3.22

function AR1Bootstrap(obs, bootstraps)

    if (nargin < 1)
        obs = 50;
    end
    
    if (nargin < 2)
        bootstraps = 200;
    end
    
    phi     = 0.99;
    phi_0   = 0.98; % GN initial guess
    theta   = [];
    
    X       = generateARMA(obs, phi, theta);
    
    %estimate the coefficient
    phi_h   = fitARMAgn(X, phi_0, [], 15);
    % variance estimate
    var_h   = (1 - phi_h^2)/obs;
    
    % build innovations
    eps     = zeros(1,obs);
    eps(1)  = X(1)*sqrt(1-phi_h^2);
    for (i = 2:obs)
        eps(i) = X(i) - phi_h*X(i-1);
    end
    
    % run the bootstrap experiment
    phi_hb  = zeros(1,bootstraps);
    var_hb  = zeros(1,bootstraps);
    x_b     = zeros(1,obs);
    for (b = 1:bootstraps)
        index = ceil(rand*obs);
        x_b(1) = eps(index)/sqrt(1-phi_h^2);
        for (j = 2:obs)
            index = ceil(rand*obs);
            x_b(j) = phi_h*x_b(j-1) + eps(index);
        end
        phi_hb(b) = fitARMAgn(x_b, phi_0, [], 15);
    end

    %plot
    [n,x] = hist(phi_hb, 50);
    bar(x,n/trapz(x,n),'BarWidth',1) %convert hist to pdf (divide by area)
    hold on
    plot(x,normpdf(x,phi_h,sqrt(var_h)),'r')
    hold off
    
    disp(['Asymptotic distribution ~ N(' num2str(phi_h) ',' num2str(var_h) ')'])
    
return