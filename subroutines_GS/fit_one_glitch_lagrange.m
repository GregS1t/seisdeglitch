function [theta, yhat, drift] = fit_one_glitch_lagrange(data2fit, weight, index, Greens, ND2, NDS, typefit, lagrangemult)
    %
    % Function to fit a single glitch on a given axe at a given index.
    % The fit is done using the Greens function and it's derivatives
    %        + a linear trend.
    %
    % This version of the fit integrate Lagrange multipliers : first and last 
    % points of the fit are constrainted to be equal to the first and last
    % points of the input time series.
    % Explanation : https://people.duke.edu/~hpgavin/cee201/constrained-least-squares.pdf
    %
    % INPUTS:
    % -------
    %      @data2fit : 1xN array which contains data vector                      
    %      @weight    : Weight vector to be applyed on the fit (same size as the Green functions : (1,  ND2+NDS)     
    %      @index: Integer with the index number of the glitch to fit   
    %      @Greens     : Array of Greens function defined in the function
    %                              get_green_func ( size (4, NB2+NDS + 1))
    %      @ND2        : Size before max
    %      @NDS        : Size after max
    %      @typefit     : if "g" only Green's and derivative
    %                           "gp" fit with Greens and Precursor
    %      @lagrangemult : 0/1 Flag to constrain the limit of the fit
    % OUTPUTS
    % -------- 
    %       @theta     : vector of fitted parameters
    %       @yhat       : fitted signal (same size as data2fit)
    
    % Window around the extremum
    y = data2fit(index-ND2:index+NDS);
    
    % Fitting  : Weighted Normal equation : $\theta=(X^T.W.X)^{-1}.(X^T.W.y)$
    % -----------------------------------------------------------------

    x1 = Greens(1,:).';
    x2 = Greens(2,:).';
    x3 = Greens(3,:).';
    x4 = Greens(4,:).';
    x5 = ones(ND2+NDS+1, 1);
    x6 = [0:ND2+NDS]'*1.0;
    if typefit == "g"                          % only the green functions
        X = [x1, x2, x5, x6];
        A = [[Greens(1,1), Greens(2,1), x5(1), x6(1)];...
                [Greens(1, end), Greens(2, end), x5(end), x6(end)]];
    elseif typefit == "gp"                  % green and precursor
        X = [x1, x2, x3, x4, x5, x6];
        A = [[Greens(1,1), Greens(2,1), Greens(3,1), Greens(4,1), x5(1), x6(1)];...
                [Greens(1, end), Greens(2, end), Greens(3, end), Greens(4, end), x5(end), x6(end)]];
    end
    b = [y(1); y(end)];                   % force the fisrt and end points to be 
    %b = [0; 0];                           % equal to those of the input signal

    % Linear system to solve in case of Lagrange multipliers
    % [ 2X'X    A' ] [ a      ]  = [ 2X'y ]
    % [   A       0 ] [  lbda ]     [ b      ]
    if lagrangemult == 1
        M = [[2*X'.*weight*X, A'];[A, zeros(2)]];
        R = [2*X'.*weight*y ; b];

        theta = M\R;        % theta = inv(M)*R
        theta = theta(1:size(X,2));
    else
        theta = ((X'.*weight*X)\(X'.*weight*(y)));   
    end
    yhat = X*theta; % Fitted model
    drift = [x5,x6]*theta(5:6);    
end