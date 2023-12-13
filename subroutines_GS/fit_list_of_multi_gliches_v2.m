function [struct_data, Variance, Glitch_cat] =  fit_list_of_multi_gliches(struct_data, clusters_of_glitches, ...
                                                        max_glitch_per_axis, green_fn_list, Variance, dt, typefit, lagrangemult, weight, ND2, NDS, ...
                                                        Glitch_cat, var_thres, verbose)
% Function to fit a list of cluster of gliches 
% This version fit the glitches one by one,  in the decreasing of amplitude
%
% -----
% INPUT:
%       @struct_data : struct array with all the input data
%       @clusters_of_glitches : array with the indexes of the maxima of the glitches
%       @max_glitch_per_axis  : array with the max number of glitches found per axis 
%                                           (used to initialize the Variance array)
%       @green_fn_list             : Green's functions and derivatives used in the fit_one_glitch_lagrange function 
%       @Variance                   : Array with the evolution of the variance from the previous fits
%       @typefit                      : can be either "g" (Just Green and derivative) 
%                                                          or "gp" (Green and first, second and third derivatives)
%       @lagrangemult            : 0/1 to constraint the first and the last point of the fit to be identical to the signal ones
%       @ND2, @NDS              : respectivaly left and right offsets of the green functions
%       @Glitch_cat                 = catalog of glitches
% -------
% OUPUT:
%       @struct_data : enlarged struct array with the result of all
%               fits and the deglitched signal for each axes.
%       @Variance: Array to store the reduction of variance at each
%               step of the fit (after removing the fitted glitch).
%       @Glitch_cat                 = catalog of glitches updated
%       @verbose                   = 0/1 Display informations

    length_results_fit_init = length(Glitch_cat);
   
    for i_ax=1:3 
        %disp(strcat("i_ax = ", num2str(i_ax)))
        iter_idx = length_results_fit_init;

        for clust_idx=2:length(clusters_of_glitches) % Loop over all the clusters

            % Sort the glitches by descending amplitude
            absampS =  sort(abs(struct_data(i_ax).data_dg(clusters_of_glitches{clust_idx})), 'descend');
            sort_idx = [];
            for i=1:length(absampS)
                id = find(abs(struct_data(i_ax).data_dg(clusters_of_glitches{clust_idx})) == absampS(i));
                sort_idx = [sort_idx, id];
            end

            % Fit each glitches as a single glitch
            for gl=1:length(sort_idx)
                idx = sort_idx(gl);                                        % Need to decompose each cluster fit
                %disp([idx,  clusters_of_glitches{clust_idx}(idx), struct_data(i_ax).data(clusters_of_glitches{clust_idx}(idx))])
                index =  clusters_of_glitches{clust_idx}(idx);

                if index + NDS <= length(struct_data(i_ax).data_dg)
                    [theta, yhat] = fit_one_glitch_lagrange(struct_data(i_ax).data_dg, weight, ...
                                              index, green_fn_list, ND2, NDS, typefit, lagrangemult);                                                                                                    
                    y                = struct_data(i_ax).data_dg(index-ND2:index+NDS);
                    residual      = y - yhat;
                    
                    % If var(y-yhat)/var(y) < var_thres
                    struct_data(i_ax).data_dg(index-ND2:index+NDS) = residual;
                    
                    %else
                    % struct_data(i_ax).data_dg(index-ND2:index+NDS) = y;
                    % end
                    [~, ampEx] = find_extremum(yhat);
                    
                    Glitch_cat(iter_idx).index                                                               = index;
                    Glitch_cat(iter_idx).time                                                                = struct_data(i_ax).time(index);
                    Glitch_cat(iter_idx).lmst                                                                = struct_data(i_ax).lmst(index);
                    Glitch_cat(iter_idx).(strcat("data_",struct_data(i_ax).label))               = struct_data(i_ax).data(index);
                    Glitch_cat(iter_idx).(strcat("maxAmpFit_",struct_data(i_ax).label))     = ampEx;
                    Glitch_cat(iter_idx).(strcat("thetaG_",struct_data(i_ax).label))            = theta(1);
                    Glitch_cat(iter_idx).(strcat("delayG_",struct_data(i_ax).label))            = theta(2)/theta(1); 
                    Glitch_cat(iter_idx).(strcat("var_",struct_data(i_ax).label))                  = var(y);  
                    Glitch_cat(iter_idx).(strcat("var_res_",struct_data(i_ax).label))           = var(y-yhat);
                    Glitch_cat(iter_idx).(strcat("var_ratio_",struct_data(i_ax).label))        = var(y-yhat)/var(y);
                    Glitch_cat(iter_idx).(strcat("var_tot_",struct_data(i_ax).label))           = var(y);  
                    Glitch_cat(iter_idx).(strcat("var_tot_res_",struct_data(i_ax).label))    = var(y-yhat); 
                    Glitch_cat(iter_idx).(strcat("var_tot_ratio_",struct_data(i_ax).label)) = var(y-yhat)/var(y);
                    if typefit == "gp"
                       Glitch_cat(iter_idx).(strcat("thetaP_",struct_data(i_ax).label))       = theta(3);   
                       Glitch_cat(iter_idx).(strcat("delayP_",struct_data(i_ax).label))       = theta(4)/theta(3);

                    end

                    Glitch_cat(iter_idx).comment = strcat("Polyglitch_", num2str(clust_idx),"_",num2str(gl));
                    % Stack of the synthetic glitches
                    struct_data(i_ax).synth_glitch(index-ND2:index+NDS) = yhat;  % Create a synthetic signal of glitches
                                                                                                               % One put 0.0 between two glitches.  
                    iter_idx = iter_idx+1;                                             
                else
                    disp(["Glitch # ",index, " rejected."])
                end
            end
        end
    end
end