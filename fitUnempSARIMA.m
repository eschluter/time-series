% Erik Schluter
% 8/10/2015
% Shumway and Stoffer - 3.37

function fitUnempSARIMA(horizon)

    if (nargin < 1)
        horizon = 12;
    end

    X = load('datasets\unemp_series.txt');
    n = length(X);
    
    % plot ACF and PACF of raw data
    % compare my PACF func to built in parcorr func for kicks
    my_pacf = PACF(X);
    figure(1)
    subplot(2,1,1); plot(ACF(X), 'k');
    title('ACF of raw unemployment data')
    subplot(2,1,2)
    hold on
    plot(parcorr(X, n/2-5), 'k')
    plot(2:length(my_pacf)+1, my_pacf, 'r-o')
    title('PACF of raw unemployment data')
    hold off
    figure(2)
    plot(X)
    title('Unemployment data')
    
    %take first difference
    Xd = diff(X);
    figure(3)
    subplot(2,1,1); plot(ACF(Xd), 'k');
    title('first differenced ACF')
    subplot(2,1,2); plot(PACF(Xd), 'k')
    title('first differenced PACF')
    figure(4)
    plot(Xd, 'k')
    title('First differenced unemployment data')
    
    %take seasonal difference
    SXd = seasonalDiff(Xd, 12);
    nsd = length(SXd);
    figure(5)
    subplot(2,1,1); plot(ACF(SXd), 'k'); 
    title('first+seasonal differenced ACF')
    subplot(2,1,2); plot(PACF(SXd), 'k');
    title('first+seasonal differenced PACF')
    
    % build data structure to fit with regression
    % Fit ARIMA (2,1,0) x (3,1,1)_12
    numParams = 6;
    SXd_trunc = SXd(37:end)';
    SXd_trunc = SXd_trunc - mean(SXd_trunc);
    D         = zeros(numParams, nsd-36);
    D(6,:)    = randn(1,nsd-36); % THETA_1
    D(5,:)    = SXd(1:end-36);  % PHI_3
    D(4,:)    = SXd(13:end-24); % PHI_2
    D(3,:)    = SXd(25:end-12); % PHI_1
    D(2,:)    = SXd(35:end-2);  % phi_2
    D(1,:)    = SXd(36:end-1);  % phi_1
    
    % solve normal equations
    Z         = D';
    b         = (Z'*Z)\(Z'*SXd_trunc);
    
    % calculate fit
    fit =    b(1)*D(1,:) + b(2)*D(2,:) + ...                  phi's
             b(3)*D(3,:) + b(4)*D(4,:) + b(5)*D(5,:) + ...    PHI's
             b(6)*D(6,:);                                   % THETA
             
    figure(6)
    hold on
    plot(SXd, 'k')
    plot(37:length(SXd), fit, 'r')
    title('first+seasonal differenced data (black), SARIMA fit (red), forecast (green)')
    
    disp('ARIMA (2,1,0) x (3,1,1)_12')
    disp(['[phi_1 phi_2 PHI_1 PHI_2 PHI_3 THETA_1] = ' mat2str(b,4)])
    
    % forecast
    wt = [D(6,:) randn(1,horizon)];
    for i = 1:horizon
        fit_next = b(1)*SXd(end-1) + b(2)*SXd(end-2) + ...
                   b(3)*SXd(end-12) + b(4)*SXd(end-24) + b(5)*SXd(end-36) + ...
                   b(6)*wt((nsd-36+horizon)-12);
        SXd = [SXd fit_next];
    end
    start           = length(SXd)-horizon;
    plot(start:start+horizon, SXd(start:end), 'g')
    hold off
    
    % invert differences
    complete_SXd    = [Xd(1:12)' SXd];
    iSfit           = seasonalCumSum(complete_SXd, 12);
    complete_iSfit  = [X(1) iSfit];
    iSDfit          = cumsum(complete_iSfit);
    transN          = length(iSDfit);
    figure(7)
    hold on
    plot(iSDfit(1:(end-horizon)), 'k')
    plot((transN-horizon):transN, iSDfit((end-horizon):end), 'r')
    title('Cumsum inverse transformed data to original data set, Forecast (red)')
    hold off
    
return

function [out] = seasonalDiff(data, order)
    n = length(data);
    out = zeros(1, n-12);
    for i = 1:(n-order)
        out(i) = data(order+i) - data(i);
    end
return

function [out] = seasonalCumSum(data, order)
    % data input must contain the first (order) elements from the original
    %  differenced data.
    n = length(data);
    out = zeros(1, n);
    out(1:order) = data(1:order);
    out((order+1):(2*order)) = out(1:order)+data((order+1):(2*order));
    for i = (2*order+1):n
        out(i) = out(i-order) + data(i);
    end
return