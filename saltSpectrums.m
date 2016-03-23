% Erik Schluter
% 8/28/2015
% Shumway and Stoffer - 4.9

function saltSpectrums(alpha)

    if (nargin < 1)
        alpha = 0.05;
    end

    xsC = load('datasets\salt.txt');
    xsT = load('datasets\saltemp.txt');
    
    nsC = length(xsC);
    nsT = length(xsT);
    
    % Find chi2 values at the specified alpha level
    aLower = chi2inv(1-alpha/2, 2);
    aUpper = chi2inv(alpha/2, 2);
    
    % must divide number of data points to account for scaling
    % salt concentration spectrum
    per_C = DFT(xsC)/nsC;
    fT_C  = fft(xsC)/nsC;
    
    % Express spectrum in terms of absolute values
    Iwp_C = abs(per_C(2:floor(nsC/2)));
    Iwf_C = abs(fT_C(2:floor(nsC/2))).^2;
    
    % Calculate (1-alpha)% confidence intervals for highest peak
    [peakMagC, pIndC] = sort(Iwp_C, 'descend');
    ciC = [ (2*peakMagC(1))/aLower (2*peakMagC(1))/aUpper ];
    disp([num2str(100*(1-alpha)) '% confidence interval for peak at ' num2str(pIndC(1))])
    disp(mat2str(ciC))
    
    % salt temp spectrum
    per_T = DFT(xsT)/nsT;
    fT_T  = fft(xsT)/nsT;
    
    % Express spectrum in terms of absolute values
    Iwp_T = abs(per_T(2:floor(nsT/2)));
    Iwf_T = abs(fT_T(2:floor(nsT/2))).^2;
    
    % Calculate (1-alpha)% confidence intervals for highest peak
    [peakMagT, pIndT] = sort(Iwp_T, 'descend');
    ciT = [ (2*peakMagT(1))/aLower (2*peakMagT(1))/aUpper ];
    disp([num2str(100*(1-alpha)) '% confidence interval for peak at ' num2str(pIndT(1))])
    disp(mat2str(ciT))
    
    figure(1)
    hold on
    plot(Iwp_C, 'or')
    plot(Iwf_C, 'k')
    xlabel('Fundamental frequencies')
    ylabel('|d(w_j)|^2')
    title('DFT using FFT (black) and covariance (red circles overlayed) of soil concentration data')
    hold off
    
    figure(2)
    hold on
    subplot(2,1,1); plot(xsC, 'k'); title('raw soil concentration data')
    subplot(2,1,2); plot(ACF(xsC), 'k'); title('ACF of raw soil concentration data')
    hold off
    
    figure(3)
    hold on
    plot(Iwp_T, 'or')
    plot(Iwf_T, 'k')
    xlabel('Fundamental frequencies')
    ylabel('|d(w_j)|^2')
    title('DFT using FFT (black) and covariance (red circles overlayed) of soil temperature data')
    hold off
    
    figure(4)
    hold on
    subplot(2,1,1); plot(xsT, 'k'); title('raw soil temperature data')
    subplot(2,1,2); plot(ACF(xsT), 'k'); title('ACF of raw soil temperature data')
    hold off
    
return