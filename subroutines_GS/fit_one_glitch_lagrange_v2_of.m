function [theta, yhat, drift] = fit_one_glitch_lagrange_of(data2fit, weight, index, Greens, ND2, NDS, typefit, lagrangemult, iax, sub_glitch_cat)
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
    %      @iax                : specify the axis number (1, 2, 3 for VBBU, VBBV, VBBW)
    % OUTPUTS
    % -------- 
    %       @theta     : vector of fitted parameters
    %       @yhat      : fitted signal (same size as data2fit)
    %       @drift      : fitted slope under the glitch 
    %disp("In fit one glitch Lagrange v2")
    % Window around the extremum
    y = data2fit(index-ND2:index+NDS);
    
    % Fitting  : Weighted Normal equation : $\theta=(X^T.W.X)^{-1}.(X^T.W.y)$
    % -----------------------------------------------------------------

    x1 = Greens(1,:).';                 % G
    x2 = Greens(2,:).';                 % G'
    x3 = Greens(3,:).';                 % P
    x4 = Greens(4,:).';                 % P'
    x5 = ones(ND2+NDS+1, 1);  
    x6 = [0:ND2+NDS]'*1.0;
  
    %disp(strcat("This glitch (", num2str(index),") has been first fitted on axes ", num2str(sub_glitch_cat.first_fit_axs)))
    ax = ["u", "v", "w"];
    a = sub_glitch_cat.(strcat("thetaG_", ax(sub_glitch_cat.first_fit_axs)));
    b = sub_glitch_cat.(strcat("thetaGp_", ax(sub_glitch_cat.first_fit_axs)));
    b_sur_a =b/a;

    X = [x1+ b_sur_a*x2, x3+ b_sur_a*x4, x5, x6];
    A2 = [[Greens(1,1) + b_sur_a*Greens(2,1), Greens(3,1) + b_sur_a*Greens(4,1), x5(1), x6(1)];...
            [Greens(1, end) + b_sur_a*Greens(2, end), Greens(3, end) + b_sur_a*Greens(4, end), x5(end), x6(end)]];
    b = [y(1); y(end)];   

    if lagrangemult == 1
        M = [[2*X'.*weight*X, A2'];[A2, zeros(2)]];
        R = [2*X'.*weight*y ; b];

        theta2 = M\R;        % theta = inv(M)*R
        theta2 = theta2(1:size(X,2));
    else
        theta2 = ((X'.*weight*X)\(X'.*weight*(y)));   
    end
    y_hat_tmp =  X*theta2;
    drift = [x5, x6]*theta2(3:4);

    residual = y - y_hat_tmp;
    y_no_drift    = y - drift;
    yhat = X*theta2; % Fitted model
    
    theta = [theta2(1), theta2(1)*b_sur_a, theta2(2), theta2(2)*b_sur_a, theta2(3), theta2(4)];
end