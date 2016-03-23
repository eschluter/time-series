function [beta] = fitARMAgn(xt, e_phi, e_theta, iterations)
%FITARMAGN fit an ARMA(p,q) model to a time series
%
%   [beta] = fitARMAgn(xt, e_phi, e_theta, iterations)
%
%   takes four arguments:
%       xt         - array of data to fit
%       e_phi      - array of AR paramameter estimates
%       e_theta    - array of MA parameter estimates
%       iterations - number of Gauss-Newton iterations (default = 12)
%
%   Erik Schluter - 8/3/2015
    
    if (nargin < 4)
        iterations = 12;
    end

    % fit model using Gauss-Newton
    np      = length(e_phi);
    nq      = length(e_theta);
    start   = max(np,nq) + 1;
    n       = length(xt);
    beta    = [e_phi e_theta]';
    z       = zeros(n,np+nq);
    w       = zeros(n,1);
    
    for (j = 1:iterations)
        % calculate error
        for (t = start:n)
            %AR part
            ar_term = 0;
            for (i = 1:np)
                ar_term = ar_term + beta(i)*xt(t-i);
            end
            %MA part
            ma_term = 0;
            for (i = (np+1):(np+nq))
                ma_term = ma_term + beta(i)*w(t-(i-np));
            end
            w(t) = xt(t) - ar_term - ma_term;
        end
        %calculate AR z terms
        for (t_outer = 1:np)
            for (t = start:n)
                z(t,t_outer) = xt(t-1) + beta(t_outer)*z(t-t_outer,t_outer);
            end
        end
        %calculate MA z term
        for (t_outer = (np+1):(np+nq))
            for (t = start:n)
                z(t,t_outer) = w(t-1) + beta(t_outer)*z(t-(t_outer-np),t_outer);
            end
        end
        
        beta = beta + (z'*z)\(z'*w);
    end
    
return