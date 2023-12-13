function [theta, yhat, drift] = fit_one_glitch_lagrange_ff(data2fit, weight, index, Greens, ND2, NDS, typefit, lagrangemult, iax)
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
    x5 = ones(ND2+NDS+1, 1);   % step  under the glitch
    x6 = [0:ND2+NDS]'*1.0;        % slope under the glitch
    
    % New version of fit to constrain the relation between the precursor 
    % and the green function. 
    % 1> fit the Green and its derivative -> a, b
    % 2> fit all the function with a/b = c+d
    
    A_param = [];
    B_param = [];
    C_param = [];
    D_param = [];
    Var_vec = [];
    n_iter = 3;       % # of iteration  
    
    %1> fit the Green and its derivative
    if typefit == "g"                      % only the glitch functions
        X = [x1, x2, x5, x6];
        A = [[Greens(1,1), Greens(2,1), x5(1), x6(1)];...
                [Greens(1, end), Greens(2, end), x5(end), x6(end)]];
    elseif typefit == "gp"               % glitch and precursor
        X = [x1, x2, x3, x4, x5, x6];
        A = [[Greens(1,1), Greens(2,1), Greens(3,1), Greens(4,1), x5(1), x6(1)];...
                [Greens(1, end), Greens(2, end), Greens(3, end), Greens(4, end), x5(end), x6(end)]];
    end
    b = [y(1); y(end)];                   % force the fisrt and end points to be 
    %b = [0; 0];                            % equal to those of the input signal

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
    yhat1 = X*theta; % Fitted model
    drift1 = [x5, x6]*theta(5:6);    
    
    residual = y - yhat1;
    y_no_drift    = y - drift1;
        
    wvar_res = residual'*diag(weight)*residual;
    wvar_y_no_drift = y_no_drift'*diag(weight)*y_no_drift;
    wvar_reduction  = wvar_res / wvar_y_no_drift; 
    Var_vec = [wvar_res];
    
    %2> Fit c and d parameters

    A_param = [theta(1)];
    B_param = [theta(2)];
    C_param = [0];
    D_param = [0];
    
    for idx = 2:n_iter
        b_sur_a = B_param(idx-1)/A_param(idx-1);
        
        X_2 = [x1, x2, x3+ b_sur_a*x4, x5, x6];
        A2 = [[Greens(1,1), Greens(2,1), Greens(3,1)+ b_sur_a*Greens(4,1), x5(1), x6(1)];...
                [Greens(1, end), Greens(2, end), Greens(3, end) + b_sur_a*Greens(4, end), x5(end), x6(end)]];
            
        b = [y(1); y(end)];   

        if lagrangemult == 1
            M = [[2*X_2'.*weight*X_2, A2'];[A2, zeros(2)]];
            R = [2*X_2'.*weight*y ; b];

            theta2 = M\R;        % theta = inv(M)*R
            theta2 = theta2(1:size(X_2,2));
        else
            theta2 = ((X_2'.*weight*X_2)\(X_2'.*weight*(y)));   
        end
        A_param = [A_param, theta2(1)];
        B_param = [B_param, theta2(2)];
        C_param = [C_param, theta2(3)];
        D_param = [D_param, C_param*b_sur_a];
        y_hat_tmp =  X_2*theta2;
        drift = [x5, x6]*theta2(4:5);
        
        residual = y - y_hat_tmp;
        y_no_drift    = y - drift;
        
%         wvar_res = residual'*diag(weight)*residual;
%         wvar_y_no_drift = y_no_drift'*diag(weight)*y_no_drift;
%         wvar_reduction  = wvar_res / wvar_y_no_drift;    % Weighted norm ratio
%         Var_vec = [Var_vec, wvar_res];
        
%          if index == 176469 %|| index == 20657 || index == 168734 || index == 21407
% 
%             figure(200000+index)
%             subplot(3,3, iax)
%             plot(Var_vec, 'DisplayName', "Var residual")
%             hold on 
%             xlabel("Iteration")
%             ylabel("Variance residual")
%             title(strcat("Var residual on axe ", num2str(iax)))
%             grid on
%             subplot(3,3, iax+3)
%             plot(A_param, 'DisplayName', "A param")
%             title(strcat("A on axe ", num2str(iax)))
%             hold on 
%             grid on
%             subplot(3,3, +6+iax)
%             plot(B_param, 'DisplayName', "B param")
%             title(strcat("B on axe ", num2str(iax)))
%             hold on 
%             grid on
%          end

        % a'; b', c', d'
    end % end itération to fit a, b, c, d params
    
    yhat = X_2*theta2; % Fitted model
    drift = [x5, x6]*theta2(4:5);
    
    theta = [theta2(1), theta2(2), theta2(3), theta2(3)*b_sur_a, theta2(4), theta2(5)];
    
%     if index == 176469
%         figure(index+iax)
%         plot(yhat, 'DisplayName', "Fit")
%         hold on 
%         plot(y+1000, 'DisplayName', "Real")
%         plot(y-yhat+2000, 'DisplayName', "Residual")
%         
%         xlabel("samples")
%         ylabel("DU")
%         legend
%         grid on
%     end
    
end