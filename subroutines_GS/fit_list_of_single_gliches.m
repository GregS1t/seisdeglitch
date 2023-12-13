function [struct_data, Glitch_cat] =  fit_list_of_single_gliches(struct_data, index_list_2_deglitch, ...
                                                        max_glitch_per_axis, green_fn_list, typefit, lagrangemult, weight, ...
                                                        var_thres, ND2, NDS, verbose)
    % Function to fit a full list of single extrema
    % Roughly, it's a double loop over the axes and over the list of
    % glitches
    %
    % -----
    % INPUT:
    %       @struct_data : struct array with all the input data
    %       @index_list_2_deglitch : array with the indexes of the maxima of the glitches
    %       @max_glitch_per_axis  : array with the max number of glitches found per axis 
    %                                           (used to initialize the Variance array)
    %       @green_fn_list             : Green's functions and derivatives used in the fit_one_glitch_lagrange function 
    %       @typefit                      : can be either "g" (Just Green and derivative) 
    %                                                          or "gp" (Green and first, second and third derivatives)
    %       @lagrangemult            : 0/1 to constraint the first and the last point of the fit to be identical to the signal ones
    %       @weight                      : weight vector to overwight the maximum part of the Greens func'
    %       @ND2, @NDS              : respectivaly left and right offsets of the green functions
    %       @results_fit                  : Array with the results of the parameters
    %
    % -------
    % OUPUT:
    %       @struct_data : enlarged struct array with the result of all
    %       fits and the deglitched signal for each axes.
    %       @Variance: Array to store the reduction of variance at each
    %       step of the fit (after removing the fitted glitch).
    %       @results_fit : Array with the results of the parameters,
    %       updated after the fits

     % Initialize Variance array                                           
     if verbose == 1
        disp("In fit_list_of_single_glitches function.")
     
     end
     
    Variance = ones(max_glitch_per_axis, 3)*2.;
    iter_idx = 1; 
    
    
    %% New section to sort glitches by amplitudes over the 3 axes 
 
    
    %%
    
    
    for i_ax=1:3
        % define an array to stack and store the glitches
        struct_data(i_ax).synth_glitch = zeros(length(struct_data(i_ax).data), 1);
        struct_data(i_ax).drift = zeros(length(struct_data(i_ax).data), 1);
        struct_data(i_ax).data_dg = struct_data(i_ax).data; 
       
        % One removes the biggest glitches first so we
        % have to sort the peaks according their amplitude. 
        [amp_Peak, idx_Peak] = sort(struct_data(i_ax).data(index_list_2_deglitch), 'descend', 'ComparisonMethod', 'abs');

        %Parameters of the loop
        idxmin  = 1;
        idxmax = length(amp_Peak);
       
        for idx_glitch=idxmin:idxmax
            index = index_list_2_deglitch(idx_Peak(idx_glitch));
            %disp(strcat("Index", num2str(index)))
            if index + NDS <= length(struct_data(i_ax).data_dg) && index - ND2 > 0
            
                [theta, yhat, drift] = fit_one_glitch_lagrange_v2_ff(struct_data(i_ax).data_dg, weight, ...
                                                  index, green_fn_list, ND2, NDS, typefit, lagrangemult, i_ax);

                y = struct_data(i_ax).data_dg(index-ND2:index+NDS);
                residual  = y - yhat;
                y_no_drift    = y - drift;
                 
                % Var estimation. 
                % var_reduction = var(y-yhat)/var(y);                 % no more used 
                wvar_res = residual'*diag(weight)*residual;
                wvar_y_no_drift = y_no_drift'*diag(weight)*y_no_drift;
                wvar_reduction  = wvar_res / wvar_y_no_drift;    % Weighted norm ratio
                

                %struct_data(i_ax).data_dg(index-ND2:index+NDS) = residual;
                %struct_data(i_ax).drift(index-ND2:index+NDS) = drift; 
                [~, ampEx] = find_extremum(yhat);
                
                if exist('Glitch_cat', 'var') == 0
                    idx = iter_idx;
                    %iter_idx = iter_idx +1;
                else
                    gl_in_cat = find([Glitch_cat.index] == index);
                    if ~isempty(gl_in_cat)
                        idx = gl_in_cat;
                    else
                        idx = iter_idx+1;
                        iter_idx = iter_idx +1;
                    end
                end     
                Glitch_cat(idx).index            = index;
                Glitch_cat(idx).time              = struct_data(i_ax).time(index);
                Glitch_cat(idx).lmst              = struct_data(i_ax).lmst(index);
                Glitch_cat(idx).(strcat("data_",struct_data(i_ax).label))      = struct_data(i_ax).data(index);
                Glitch_cat(idx).(strcat("maxAmpFit_",struct_data(i_ax).label))    = ampEx;
                Glitch_cat(idx).(strcat("thetaG_",struct_data(i_ax).label))    = theta(1);
                Glitch_cat(idx).(strcat("thetaGp_",struct_data(i_ax).label))  = theta(2); 
                Glitch_cat(idx).(strcat("delayG_",struct_data(i_ax).label))  = theta(2)/theta(1); 
                Glitch_cat(idx).(strcat("var_tot_ratio_",struct_data(i_ax).label)) = wvar_reduction;
                if typefit == "gp"
                   Glitch_cat(idx).(strcat("thetaP_",struct_data(i_ax).label))   = theta(3);   
                   Glitch_cat(idx).(strcat("thetaPp_",struct_data(i_ax).label)) = theta(4);  
                   Glitch_cat(idx).(strcat("delayP_",struct_data(i_ax).label))  = theta(4)/theta(3);

                end

                if wvar_reduction<= var_thres
                    Glitch_cat(idx).(strcat("isglitch_",struct_data(i_ax).label)) = 1;
                    struct_data(i_ax).drif(index-ND2:index+NDS) = drift;
                    struct_data(i_ax).data_dg(index-ND2:index+NDS) = residual;
                    struct_data(i_ax).synth_glitch(index-ND2:index+NDS) = yhat;      
                else
                    Glitch_cat(idx).(strcat("isglitch_",struct_data(i_ax).label)) = 0;
                    struct_data(i_ax).drif(index-ND2:index+NDS) = drift;
                    %struct_data(i_ax).data_dg(index-ND2:index+NDS) = struct_data(i_ax).data(index);
                    %struct_data(i_ax).synth_glitch(index-ND2:index+NDS) = ...
                    %                                    zeros(length(struct_data(i_ax).data(index)), 1);   

                end % end if sur la variance
                
                Glitch_cat(idx).comment = "SingleGlitch";
                % Stack of the synthetic glitches
                % Create a synthetic signal of glitches
                                                                                                          % One put 0.0 between two glitches.  
                %iter_idx = iter_idx+1; 
                plot_glitch = 0;
                if plot_glitch == 1
                    if index == 2291 || index == 20657 || index == 168734 
                        disp(strcat("Glitch # ", num2str(index), " - a= ", num2str(theta(1)),   " - b/a= ", num2str(theta(2)/theta(1))))
                        figure(index)
                        tmp_vec = index-ND2:index+NDS;
                        hold on 
                        subplot(3, 1, i_ax)
                        %xlim([min_idx-ND2, max_idx+NDS])
                        plot(tmp_vec, y, "DisplayName", "y")
                        hold on
                        plot(tmp_vec, yhat, "DisplayName", "yhat")
                        plot(tmp_vec, residual, "DisplayName", "residual")
                        grid on
                        xlabel("Samples")
                        ylabel("Amplitudes")
                        legend
                    end
                end

            else
                disp(["Glitch # ",index, " rejected."])
            end
        end %end loop over extrema
    end
    
end