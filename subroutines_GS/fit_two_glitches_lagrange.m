function [theta, yhat] = fit_two_glitches_lagrange(data2fit, weight, index1, index2, Greens, ND2, NDS, typefit, lagrangemult)
    %
    % Function to fit two close glitches on a given axe.
    % The fit is done using the Greens function and it's derivatives
    %        + a linear trend + lagrange multipliers
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
    %      @index1    : Integer with the index number of the first glitch to fit   
    %      @index2    : Integer with the index number of the second glitch to fit   
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
    
    % Window around the both extrema    
    y = data2fit(index1-ND2:index2+NDS);
    
    % Fitting  : Weighted Normal equation : $\theta=(X^T.W.X)^{-1}.(X^T.W.y)$
    % -----------------------------------------------------------------
    
    % Green functions
    x1 = Greens(1,:).';
    x2 = Greens(2,:).';
    x3 = Greens(3,:).';
    x4 = Greens(4,:).';
    
    
    x1_l = [ones(abs(index2-index1),1)*x1(1); x1];
    x2_l = [ones(abs(index2-index1),1)*x2(1); x2];
    x3_l = [ones(abs(index2-index1),1)*x3(1); x3];
    x4_l = [ones(abs(index2-index1),1)*x4(1); x4];
    
    x1_r = [x1;  ones(abs(index2-index1),1)*x1(end)];            % Mettre 0 à la place
    x2_r = [x2;  ones(abs(index2-index1),1)*x2(end)];
    x3_r = [x3;  ones(abs(index2-index1),1)*x3(end)];
    x4_r = [x4;  ones(abs(index2-index1),1)*x4(end)];
    % Slope under signal
    x5 = ones(ND2+NDS+(index2-index1)+1, 1);
    x6 = [0:ND2+NDS+ (index2-index1)]'*1.0;
    
    if typefit == "g"                                                                % only the green functions
        X = [x1_l, x2_l, x1_r, x2_r, x5, x6];
        A = [[x1_l(1), x2_l(1), x1_r(1), x2_r(1), x5(1), x6(1)];...
                [x1_l(end), x2_l(end), x1_r(end), x2_r(end), x5(end), x6(end)]];
    
    elseif typefit == "gp"                                                        % green and precursor
        X = [x1_l, x2_l, x3_l, x4_l, x1_r, x2_r, x3_r, x4_r, x5, x6];
        A = [[x1_l(1), x2_l(1), x3_l(1), x4_l(1), x1_r(1), x2_r(1), x3_r(1), x4_r(1), x5(1), x6(1)];...
                [x1_l(end), x2_l(end), x3_l(end), x4_l(end), x1_r(end), x2_r(end), x3_r(end), x4_r(end), x5(end), x6(end)]];
    end
    
    b = [y(1); y(end)];                                                            % force the fisrt and end points to be
                                                                                          % equal to those of the input signal

    % Linear system to solve in case of Lagrange multipliers
    % [ 2X'X    A' ] [ a      ]  = [ 2X'y ]
    % [   A       0 ] [  lbda ]     [ b      ]
    if lagrangemult == 1
        M = [[2*X'.*weight'*X, A'];[A, zeros(2)]];
        R = [2*X'.*weight'*y ; b];

        theta = M\R;        % theta = inv(M)*R<
        theta = theta(1:size(X,2));
    else
        theta = ((X'.*weight*X)\(X'.*weight*(y)));                         
    end
    yhat =X*theta;                                                                 % Fitted model with fitted params
end