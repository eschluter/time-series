% Erik Schluter
% 7/15/2015
% Rejection sampling using g(x) = 1

function rejectionSamplingG(numSamps)

    if (nargin < 1)
        numSamps = 10000;
    end
    
    samps       = zeros(1, numSamps);
    time        = zeros(1, numSamps);
    iterations  = zeros(1, numSamps);
    
    for i = 1:numSamps
        tic
        its = 0;
        while (true)
            its = its + 1;
            X       = rand;
            alpha   = X^2;
            if (rand <= alpha)
                break;
            end
        end
        samps(i)        = X;
        time(i)         = toc;
        iterations(i)   = its;
    end
    
    timeMean  = sum(time)/numSamps;
    timeVar   = sum((time - timeMean).^2)/(numSamps - 1);
    
    itsMean  = sum(iterations)/numSamps;
    itsVar   = sum((iterations - itsMean).^2)/(numSamps - 1);
    
    disp(['Elapsed time mean = ' num2str(timeMean)])
    disp(['Elapsed time variance = ' num2str(timeVar)])
    
    disp(['Iterations mean = ' num2str(itsMean)])
    disp(['Iterations variance = ' num2str(itsVar)])
    
    [N,x]   = hist(samps,50);
    dx      = x(2) - x(1);
    plot(x,N/(numSamps*dx), 'r');
    hold on;
    plot(x,3*x.^2, 'k');
    hold off;
    
return