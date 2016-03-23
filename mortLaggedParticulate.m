% Erik Schluter
% 7/15/2015
% Shumway and Stoffer - 2.2

function mortLaggedParticulate()

    M = load('datasets\cmort_series.txt');
    T = load('datasets\tempr_series.txt');
    P = load('datasets\part_series.txt');

    T = T - mean(T);
    
    % start series at week 5 so data is aligned with lag
    pLag    = P(1:end-4);
    aM      = M(5:end);
    aT      = T(5:end);
    aP      = P(5:end);
    
    n = length(aM);
    t = 1:n;
    
    Z = [ ones(1,n); t; aT'; (aT.^2)'; aP'; pLag' ]';
    
    beta    = (Z'*Z)\(Z'*aM);
    
    mt      = beta(1) + beta(2)*(t')    + ...
                        beta(3)*aT      + ...
                        beta(4)*(aT.^2) + ...
                        beta(5)*aP      + ...
                        beta(6)*pLag;
                    
    % model fit and data
    figure(1)
    hold on
    plot(t, mt, 'r')
    plot(t, aM, 'k')
    hold off
    title('model fit and data')
    xlabel('week of the year')
    ylabel('mortality')
    
    figure(2)
    hold on
    plot(t, aM - mt)
    hold off
    title('residuals')
    
    % scatterplots
    figure(3)
    a = 20;
    hold on
    subplot(4,4,1); plot(0,0); text(-0.6,0,'Mortality')
    subplot(4,4,2); scatter(aT,aM,a); b = corrcoef(aT,aM); text(max(aT),max(aM), num2str(b(2,1)))
    subplot(4,4,3); scatter(aP,aM,a); b = corrcoef(aP,aM); text(max(aP),max(aM), num2str(b(2,1)))
    subplot(4,4,4); scatter(pLag,aM,a); b = corrcoef(pLag,aM); text(max(pLag),max(aM), num2str(b(2,1)))
    subplot(4,4,5); scatter(aM,aT,a); b = corrcoef(aM,aT); text(max(aM),max(aT), num2str(b(2,1)))
    subplot(4,4,6); plot(0,0); text(-0.85,0,'Temperature')
    subplot(4,4,7); scatter(aP,aT,a); b = corrcoef(aP,aT); text(max(aP),max(aT), num2str(b(2,1)))
    subplot(4,4,8); scatter(pLag,aT,a); b = corrcoef(pLag,aT); text(max(pLag),max(aT), num2str(b(2,1)))
    subplot(4,4,9); scatter(aM,aP,a); b = corrcoef(aM,aP); text(max(aM),max(aP), num2str(b(2,1)))
    subplot(4,4,10); scatter(aT,aP,a); b = corrcoef(aT,aP); text(max(aT),max(aP), num2str(b(2,1)))
    subplot(4,4,11); plot(0,0); text(-0.7,0,'Particulate')
    subplot(4,4,12); scatter(pLag,aP,a); b = corrcoef(pLag,aP); text(max(pLag),max(aP), num2str(b(2,1)))
    subplot(4,4,13); scatter(aM,pLag,a); b = corrcoef(aM,pLag); text(max(aM),max(pLag), num2str(b(2,1)))
    subplot(4,4,14); scatter(aT,pLag,a); b = corrcoef(aT,pLag); text(max(aT),max(pLag), num2str(b(2,1)))
    subplot(4,4,15); scatter(aP,pLag,a); b = corrcoef(aP,pLag); text(max(aP),max(pLag), num2str(b(2,1)))
    subplot(4,4,16); plot(0,0); text(-0.6,0,'P-Lagged')
    hold off
    
return