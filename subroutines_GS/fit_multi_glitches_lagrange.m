function [theta, yhat] = fit_multi_glitches_lagrange(data2fit, idx_list, Greens, dt, ND2, NDS, typefit, lagrangemult)
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
    
    idx_list = sort(idx_list);
    % Window around the both extrema    
    y = data2fit(min(idx_list)-ND2:max(idx_list)+NDS);
    
    % For convenience, Weight vector is calculated here.
    dtweight=4.;
    ntweight=floor(dtweight/dt);
    weight = ones(length(Greens(1,:))+abs(max(idx_list)-min(idx_list)),1);
    weight(ntweight:abs(max(idx_list)-min(idx_list))+5*ntweight)=100.;
    
    % Fitting  : Weighted Normal equation : $\theta=(X^T.W.X)^{-1}.(X^T.W.y)$
    % -----------------------------------------------------------------
    
    % Green functions
    x1 = Greens(1,:).';
    x2 = Greens(2,:).';
    x3 = Greens(3,:).';
    x4 = Greens(4,:).';
    
    X = [];
    A_l = [];
    A_r = [];
    for i=1:length(idx_list)
        
        x1_new = [ones(abs(idx_list(i)-min(idx_list)),1)*x1(1); x1 ; ones(abs(idx_list(i)-max(idx_list)),1)*x1(end)];
        x2_new = [ones(abs(idx_list(i)-min(idx_list)),1)*x2(1); x2 ; ones(abs(idx_list(i)-max(idx_list)),1)*x2(end)];
        X = [X, x1_new, x2_new];
        A_l = [A_l, [x1_new(1), x2_new(1)]];
        A_r = [A_r, [x1_new(end), x2_new(end)]];
        
                
        if typefit == "gp"     
            
            x3_new = [ones(abs(idx_list(i)-min(idx_list)),1)*x3(1); x3; ones(abs(idx_list(i)-max(idx_list)),1)*x3(end)];
            x4_new = [ones(abs(idx_list(i)-min(idx_list)),1)*x4(1); x4; ones(abs(idx_list(i)-max(idx_list)),1)*x4(end)];
            X = [X, x3_new, x4_new];
            A_l = [A_l, x3_new(1), x4_new(1)];
            A_r = [A_r, x3_new(end), x4_new(end)];
        end
    end
   
    % Slope under whole signal
    x5 = ones(ND2+NDS+(max(idx_list)-min(idx_list))+1, 1);
    x6 = [0:ND2+NDS+ (max(idx_list)-min(idx_list))]'*1.0;
    A_l = [A_l, x5(1), x6(1)];
    A_r = [A_r, x5(end), x6(end)];
    
    A = [A_l; A_r];
    X = [X, x5, x6];
   
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
    yhat =X*theta;                                                                % Fitted model with fitted params
end
