% Erik Schluter
% 8/28/2015
%
% This example compares the speed of physical space
%	and Fourier space filters for ARMA(2,1) series
%	of differing lengths

function gaussianFilter(analysisLengths, filterWidth)

    if (nargin < 1)
        analysisLengths = [1e3 1e4 1e5 1e6];
    end
    
    if (nargin < 2)
        filterWidth = 100;
    end
    
    nAL      = length(analysisLengths);
    armaCell = cell(nAL,1);
    filtSigP = cell(nAL,1);
    filtSigF = cell(nAL,1);
    pTime    = zeros(nAL,1);
    fTime    = zeros(nAL,1);
    
    % create Gaussian
    t_g = linspace(-3,3,filterWidth);
    gt  = exp(-t_g.^2);
    
    % create ARMA(2,1) processes to convolve
    for (i = 1:nAL)
        armaCell{i} = generateARMA(analysisLengths(i),[0.5 0.1], [0.5]);
    end
    
    % perform physical space filter
    for (i = 1:nAL)
        tic;
        filtSigP{i} = zeros(1,analysisLengths(i));
        tempSig     = zeros(1,analysisLengths(i)+filterWidth);
        tempSig     = [armaCell{i} armaCell{i}(1:100)];
        for (k = 1:analysisLengths(i))
            sum = 0;
            for (ki = 1:filterWidth)
                sum = sum + tempSig(ki+k-1)*gt(ki);
            end
            filtSigP{i}(k) = sum;
        end
        pTime(i) = toc;
    end
    
    % perform fourier space filter
    for (i = 1:nAL)
        g    = [gt zeros(1,analysisLengths(i)-filterWidth)];
        tic;
        fh   = fft(armaCell{i});
        gh   = fft(g);
        temp = fh.*conj(gh);
        filtSigF{i} = ifft(temp);
        fTime(i) = toc;
    end
    
    figure(1)
    hold on
    subplot(nAL,1,1); plot(armaCell{1},'k')
    title('Unfiltered ARMA(2,1) processes')
    for (p = 2:nAL)
        subplot(nAL,1,p); plot(armaCell{p},'k')
    end
    hold off
    
    for (p = 1:nAL)
        figure(p+1)
        hold on
    	subplot(2,1,1); plot(filtSigP{p}, 'b');
        title(['Physical space: ' num2str(pTime(p)) ' seconds'])
        subplot(2,1,2); plot(filtSigF{p}, 'r');
        title(['Fourier space: ' num2str(fTime(p)) ' seconds'])
        hold off
    end
    
return