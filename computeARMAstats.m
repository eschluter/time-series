% Erik Schluter
% 8/7/2015
% Shumway and Stoffer - 3.21

function computeARMAstats(phi, theta, numRuns, obs, gn_iterations)

    if (nargin < 1)
        phi = 0.9;
    end
    
    if (nargin < 2)
        theta = 0.5;
    end

    if (nargin < 3)
        numRuns = 10;
    end

    if (nargin < 4)
        obs = 200;
    end
    
    if (nargin < 5)
        gn_iterations = 8;
    end
    
    b = zeros(2,numRuns);
    
    for (i = 1:numRuns)
        xt      = generateARMA(obs, phi, theta);
        b(:,i)  = fitARMAgn(xt, phi, theta, gn_iterations);
    end
    
    estimatorMean  = sum(b,2)/numRuns;
    estimatorVar   = sum(((b - repmat(estimatorMean,1,numRuns)).^2),2)/(numRuns - 1);
    
    disp(['The estimator mean = ' mat2str(estimatorMean)])
    disp(['The estimator variance = ' mat2str(estimatorVar)])
    
return