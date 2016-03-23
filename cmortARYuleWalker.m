% Erik Schluter
% 8/9/2015
% Shumway and Stoffer - 3.18

function cmortARYuleWalker(order)

    if (nargin < 1)
        order = 2;
    end
    
    M = load('datasets\cmort_series.txt');
    M = M-mean(M);
    n = length(M);
    
    m_acf = [1 ACF(M)]';
    
    YWM = zeros(order);
    g   = zeros(order,1);
    for (k = 0:order-1)
        g(k+1) = m_acf(k+2);
        for (j = 0:order-1)
            YWM(j+1,k+1) = m_acf(abs(k-j)+1);
        end
    end
    
    % calculate phi hat
    phi_h = YWM \ g;
    
    % calculate asymptotic variances (diagonal of inverse information matrix)
    if (order == 2)
        v_phi_h = 1 - (phi_h(2))^2;
    end
    
    %calculate fit
    xt = zeros(1, n);
    w  = randn(1, n);
    for (i = order+1:n)
        xt(i) = w(i);
        for (j = 1:order)
            xt(i)  = xt(i) + xt(i-j)*phi_h(j);
        end
    end
    
    figure(1)
    hold on
    plot(M, 'k')
    plot(xt, 'r')
    
    % get the linear regression coeffs
    phi_r = cmortAR(2, 1, true);
    
    disp(['YW coefficient estimates:           = ' mat2str(phi_h)])
    disp(['Asymptotic variance of YW estimates = ' num2str(v_phi_h)])
    
    disp(['Regression coefficient estimates    = ' mat2str(phi_r(2:end))])
    
return