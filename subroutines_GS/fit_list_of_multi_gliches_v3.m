function [struct_data, Glitch_cat] =  fit_list_of_multi_gliches(struct_data, clusters_of_glitches, ...
                                                        max_glitch_per_axis, green_fn_list, dt, typefit, lagrangemult, weight, ND2, NDS, ...
                                                        Glitch_cat, var_thres, verbose,  plot_fit_clusters)
% Function to fit a list of cluster of gliches 
% This version fit the glitches one by one,  in the decreasing of amplitude
%
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

    disp("In fit_list_of_multi_glitches function.")

    length_results_fit_init = length(Glitch_cat);
    iter_idx = length_results_fit_init;
    %disp(strcat("Index at start", num2str(iter_idx)))
    for clust_idx=2:length(clusters_of_glitches) % Loop over all the clusters
        % Sort the glitches by descending amplitude
        
        min_idx = clusters_of_glitches{clust_idx}(1);
        max_idx = clusters_of_glitches{clust_idx}(length(clusters_of_glitches{clust_idx}));
        tempvec = min_idx-ND2:max_idx+NDS;
        
        %disp(clusters_of_glitches{clust_idx})
         if plot_fit_clusters == 1
            figure (1000+clust_idx)
             cpt = 1;
         end
       
        merge_glitches_2_vec = [];
        keep_axes = [];
        keep_idx   = [];
        
        for i_ax=1:3
            keep_axes = [keep_axes, ones(1, length(clusters_of_glitches{clust_idx}))*i_ax];
            keep_idx   = [keep_idx, clusters_of_glitches{clust_idx}];
            merge_glitches_2_vec = [merge_glitches_2_vec,  struct_data(i_ax).data_dg(clusters_of_glitches{clust_idx})'];
        end
         
         disp(keep_axes);
         disp(merge_glitches_2_vec)
         disp(keep_idx)
         
         [absamp_all_Axes, sort_idx_all_axes] =  sort(abs(merge_glitches_2_vec), 'descend');
         disp(absamp_all_Axes);
         disp(sort_idx_all_axes);
         
          for gl=1:length(sort_idx_all_axes)
            disp([keep_idx(gl), merge_glitches_2_vec(gl), keep_axes(gl)])
          end
         
         
        for i_ax=1:3
            %disp(["i_ax=" i_ax])
            [absampS, sort_idx] =  sort(abs(struct_data(i_ax).data_dg(clusters_of_glitches{clust_idx})), 'descend');
            
            % Fit each glitches as a single glitch
            
            for gl=1:length(sort_idx)

                idx = sort_idx(gl);                                                 % Need to decompose each cluster fit
                %disp([idx,  clusters_of_glitches{clust_idx}(idx), struct_data(i_ax).data(clusters_of_glitches{clust_idx}(idx))])
                index =  clusters_of_glitches{clust_idx}(idx);
                %disp(strcat("Index: ", num2str(index)))
                
                gl_in_cat = find([Glitch_cat.index] == index);
                if ~isempty(gl_in_cat)
                    %disp("Glitch already in catalog.")
                    idx_glitch = gl_in_cat;
                    
                else
                     %disp("New glitch")
                    idx_glitch = iter_idx+gl;
                end
                
                %disp(strcat("i_ax = ", num2str(i_ax)))
                if index + NDS <= length(struct_data(i_ax).data_dg)  && index - ND2 > 0
                    [theta, yhat, drift] = fit_one_glitch_lagrange_v2(struct_data(i_ax).data_dg, weight, ...
                                          index, green_fn_list, ND2, NDS, typefit, lagrangemult, i_ax);                                                                                                    
                    
                    y                = struct_data(i_ax).data_dg(index-ND2:index+NDS);
                    residual      = y - yhat;
                    y_no_drift    = y - drift;
                    
                    %%% PLOT TO ANALYSE THE STEPS OF DEGLITCHING %%%
                   
                    if plot_fit_clusters == 1
                        tmp_vec = index-ND2:index+NDS;
                        hold on 
                        subplot(3, 1, cpt)
                        xlim([min_idx-ND2, max_idx+NDS])
                        plot(tmp_vec, y, "DisplayName", "y")
                        hold on
                        plot(tmp_vec, yhat, "DisplayName", "yhat")
                        plot(tmp_vec, residual, "DisplayName", "residual")
                        grid on
                        xlabel("Samples")
                        ylabel("Amplitudes")
                        legend
                    end
                    
                    % Var estimation. 
                    % var_reduction = var(y-yhat)/var(y);                 % no more used 
                    wvar_res = residual'*diag(weight)*residual;
                    wvar_y_no_drift = y_no_drift'*diag(weight)*y_no_drift;
                    wvar_reduction  = wvar_res / wvar_y_no_drift;    % Weighted norm ratio
   
                    % Generate a temporary glitch entry
                    [~, ampEx] = find_extremum(yhat);
                    Glitch_cat(idx_glitch).index                                                            = index;
                    Glitch_cat(idx_glitch).time                                                                = struct_data(i_ax).time(index);
                    Glitch_cat(idx_glitch).lmst                                                                = struct_data(i_ax).lmst(index);
                    Glitch_cat(idx_glitch).(strcat("data_",struct_data(i_ax).label))               = struct_data(i_ax).data(index);
                    Glitch_cat(idx_glitch).(strcat("maxAmpFit_",struct_data(i_ax).label))     = ampEx;
                    Glitch_cat(idx_glitch).(strcat("thetaG_",struct_data(i_ax).label))            = theta(1);
                    Glitch_cat(idx_glitch).(strcat("delayG_",struct_data(i_ax).label))            = theta(2)/theta(1); 
                    Glitch_cat(idx_glitch).(strcat("var_tot_ratio_",struct_data(i_ax).label))     = wvar_reduction;
                    if typefit == "gp"
                        Glitch_cat(idx_glitch).(strcat("thetaP_",struct_data(i_ax).label))       = theta(3);   
                        Glitch_cat(idx_glitch).(strcat("delayP_",struct_data(i_ax).label))       = theta(4)/theta(3);
                    end
                    
                    if wvar_reduction <= var_thres
                        Glitch_cat(idx_glitch).(strcat("isglitch_",struct_data(i_ax).label)) = 1;
                        Glitch_cat(idx_glitch).comment = "PolyGlitch";
                        struct_data(i_ax).drif(index-ND2:index+NDS) = drift;
                        struct_data(i_ax).data_dg(index-ND2:index+NDS) = residual;
                        struct_data(i_ax).synth_glitch(index-ND2:index+NDS) = yhat;      
                    else
                        Glitch_cat(idx_glitch).(strcat("isglitch_",struct_data(i_ax).label)) = 0;
                        Glitch_cat(idx_glitch).comment = "PolyGlitch";
                        struct_data(i_ax).drif(index-ND2:index+NDS) = drift;
                        %struct_data(i_ax).data_dg(index-ND2:index+NDS) = struct_data(i_ax).data(index);
                        %struct_data(i_ax).synth_glitch(index-ND2:index+NDS) = zeros(length(struct_data(i_ax).data(index)), 1);   

                    end % end if sur la variance
                    
                else
                    disp(["Glitch # ",index, " rejected, out of bounds."])
                end % end sur la longueur du glitch
            end
            %disp("Next cluster")
            if plot_fit_clusters == 1
                cpt = cpt+1;
            end
        end
        iter_idx = iter_idx+length(sort_idx);
    end
end