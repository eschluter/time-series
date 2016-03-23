% Erik Schluter
% 8/28/2015
% Shumway and Stoffer - 4.8

function sunspotzSpectrum(alpha)

    if (nargin < 1)
        alpha = 0.05;
    end

    xt = load('datasets\sunspotz_series.txt');
    n = length(xt);
    
    % Find chi2 values at the specified alpha level
    aLower = chi2inv(1-alpha/2, 2);
    aUpper = chi2inv(alpha/2, 2);
    
    % Use homecooked routine (Shumway equation 4.22)
    % must divide number of data points to account for scaling
    per = DFT(xt)/n;
    % Use fft package
    fT  = fft(xt)/n;
    
    % Express spectrum in terms of absolute values
    % Discard zeroth frequency to get rid of the mean value spike
    %  also so the frequencies are aligned correctly from 1 indexing
    Iwp = abs(per(2:floor(n/2)));
    Iwf = abs(fT(2:floor(n/2))).^2;
    
    % Calculate (1-alpha)% confidence intervals for 3 highest peaks
    [peakMag, pInd] = sort(Iwp, 'descend');
    ci1 = [ (2*peakMag(1))/aLower (2*peakMag(1))/aUpper ];
    ci2 = [ (2*peakMag(2))/aLower (2*peakMag(2))/aUpper ];
    ci3 = [ (2*peakMag(5))/aLower (2*peakMag(5))/aUpper ];
    disp([num2str(100*(1-alpha)) '% confidence interval for peak at ' num2str(pInd(1))])
    disp(mat2str(ci1))
    disp([num2str(100*(1-alpha)) '% confidence interval for peak at ' num2str(pInd(2))])
    disp(mat2str(ci2))
    disp([num2str(100*(1-alpha)) '% confidence interval for peak at ' num2str(pInd(5))])
    disp(mat2str(ci3))
    
    % plot starting at index 2 since the first index is the zeroth frequency
    figure(1)
    hold on
    plot(Iwp, 'or')
    plot(Iwf, 'k')
    xlabel('Fundamental frequencies')
    ylabel('|d(w_j)|^2')
    title('DFT using FFT (black) and covariance (red circles overlayed)')
    hold off
    
    figure(2)
    hold on
    subplot(2,1,1); plot(xt, 'k'); title('raw sunspot data')
    subplot(2,1,2); plot(ACF(xt), 'k'); title('ACF of raw sunspot data')
    hold off
    
return