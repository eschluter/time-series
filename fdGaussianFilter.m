% Erik Schluter
% 8/28/2015
%
% This example generates an ARMA(2,2) process and filters with
%	a filter defined as -2*t*e^(-(t^2)/2). The filter is displayed
%	in physical and fourier space.

function fdGaussianFilter(obs, filterWidth)

    if (nargin < 1)
        obs = 10000;
    end

    if (nargin < 2)
        filterWidth = 100;
    end
    
    %Constants
    FilterRangeLim = 4;
    
    % Build first derivitive gaussian filter
    t    = linspace(-FilterRangeLim,FilterRangeLim,filterWidth);
    dg   = -2*t.*exp(-(t.^2)./2);
    figure(1)
    subplot(2,1,1); plot(dg, 'k'); title('Filter in physical space')
    
    % generate ARMA(2,2) process with 10000 observations
    xt = generateARMA(obs, [0.5 0.1], [0.21 0.9]);
    figure(2)
    hold on
    subplot(2,1,1); plot(xt(1:500), 'k'); title('raw ARMA(2,2) process')
    
    % perform fourier space filter
    g       = [dg zeros(1,obs-filterWidth)];
    fh      = fft(xt);
    gh      = fft(g);
    temp    = fh.*conj(gh);
    filtXT  = ifft(temp);
    
    subplot(2,1,2); plot(filtXT(1:500), 'r'); title('fdg filtered ARMA(2,2) process')
    hold off
    
    figure(1)
    hold on
    subplot(2,1,2); plot(abs(gh(1:floor(obs/2+1))),'r'); title('Filter in Fourier space')
    hold off
    
return