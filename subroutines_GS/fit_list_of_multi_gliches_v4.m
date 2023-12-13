function [struct_data, Glitch_cat] =  fit_list_of_multi_gliches(struct_data, clusters_of_glitches, ...
                                                        max_glitch_per_axis, green_fn_list, dt, typefit, lagrangemult, weight, ND2, NDS, ...
                                                        Glitch_cat, var_thres, verbose,  plot_fit_clusters)
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

    disp("In fit_list_of_multi_glitches function.")

    length_results_fit_init = length(Glitch_cat);
    iter_idx = length_results_fit_init;
    %disp(strcat("Index at start", num2str(iter_idx)))
    for clust_idx=2:length(clusters_of_glitches) % Loop over all the clusters
        
        % Sort the glitches by descending amplitude
        min_idx = clusters_of_glitches{clust_idx}(1);
        max_idx = clusters_of_glitches{clust_idx}(length(clusters_of_glitches{clust_idx}));
        tempvec = min_idx-ND2:max_idx+NDS;
        
         if plot_fit_clusters == 1
            figure (1000+clust_idx)
             cpt = 1;
         end
       
        %% NEW VERSION OF FIT
        % Date: 2021-05-31 
        % Sort amplitude inside a cluster in descending order
        %
        merge_glitches_2_vec = [];
        keep_axes = [];
        keep_idx   = [];
        
        for i_ax=1:3
            keep_axes = [keep_axes, ones(1, length(clusters_of_glitches{clust_idx}))*i_ax];
            keep_idx   = [keep_idx, clusters_of_glitches{clust_idx}];
            merge_glitches_2_vec = [merge_glitches_2_vec,  struct_data(i_ax).data_dg(clusters_of_glitches{clust_idx})'];
        end
         
         %disp(keep_axes);
         %disp(merge_glitches_2_vec)
         %disp(keep_idx)
         
         [absamp_all_Axes, sort_idx_all_axes] = sort(abs(merge_glitches_2_vec), 'descend');
         %disp(absamp_all_Axes);
         %disp(sort_idx_all_axes);
        
         for gl=1:length(sort_idx_all_axes)
             idx_in_ts = keep_idx(sort_idx_all_axes(gl));
             gl_axis    = keep_axes(sort_idx_all_axes(gl));
             
             %disp(strcat("idx = ", num2str(idx_in_ts), " on axis ", num2str(gl_axis)))
             gl_in_cat = find([Glitch_cat.index] == idx_in_ts);
             
             if ~isempty(gl_in_cat)
                 %disp(strcat("Glitch ",num2str(idx_in_ts)," already in catalog at position ", num2str(gl_in_cat)))
                 idx_glitch2 = gl_in_cat;
             else
                 idx_glitch2 = length(Glitch_cat)+1;
                 first_fit_axe = gl_axis;
                 %disp(strcat("New glitch : ", num2str(idx_glitch2), " fitted on axis ", num2str(first_fit_axe)))
             end
             
             if idx_in_ts + NDS <= length(struct_data(gl_axis).data_dg)  && idx_in_ts - ND2 > 0
                if gl_axis == first_fit_axe
                    %disp("Premier fit sur cet index")
                    %disp([idx_glitch2, gl_axis, first_fit_axe])
                    [theta, yhat, drift] = fit_one_glitch_lagrange_v2_ff(struct_data(gl_axis).data_dg, weight, ...
                                                    idx_in_ts, green_fn_list, ND2, NDS, typefit, lagrangemult, gl_axis);    
                    Glitch_cat(idx_glitch2).first_fit_axs = first_fit_axe ; 
                else
                    %if idx_in_ts == 26191
                    %    disp("Paf limace")
                    %end
                    %disp(strcat("Je vais récupérer les paramètres de l'axe ", num2str(first_fit_axe)))
                    sub_glitch_cat = Glitch_cat(idx_glitch2);
                    [theta, yhat, drift] = fit_one_glitch_lagrange_v2_of(struct_data(gl_axis).data_dg, weight, ...
                                                    idx_in_ts, green_fn_list, ND2, NDS, typefit, lagrangemult, gl_axis, sub_glitch_cat);    

                end
                y                = struct_data(gl_axis).data_dg(idx_in_ts-ND2:idx_in_ts+NDS);
                residual      = y - yhat;
                y_no_drift    = y - drift;
       
                wvar_res = residual'*diag(weight)*residual;
                wvar_y_no_drift = y_no_drift'*diag(weight)*y_no_drift;
                wvar_reduction  = wvar_res / wvar_y_no_drift;    % Weighted norm ratio

                % Save teh results of the fits
                [~, ampEx] = find_extremum(yhat);
                Glitch_cat(idx_glitch2).index          = idx_in_ts;
                Glitch_cat(idx_glitch2).time            = struct_data(gl_axis).time(idx_in_ts);
                Glitch_cat(idx_glitch2).lmst            = struct_data(gl_axis).lmst(idx_in_ts);
                %Glitch_cat(idx_glitch2).first_fit_axs = first_fit_axe;

                Glitch_cat(idx_glitch2).(strcat("data_",struct_data(gl_axis).label))               = struct_data(gl_axis).data(idx_in_ts);
                Glitch_cat(idx_glitch2).(strcat("maxAmpFit_",struct_data(gl_axis).label))     = ampEx;
                Glitch_cat(idx_glitch2).(strcat("thetaG_",struct_data(gl_axis).label))            = theta(1);
                Glitch_cat(idx_glitch2).(strcat("thetaGp_",struct_data(gl_axis).label))            = theta(2); 
                Glitch_cat(idx_glitch2).(strcat("delayG_",struct_data(gl_axis).label))            = theta(2)/theta(1); 

                Glitch_cat(idx_glitch2).(strcat("var_tot_ratio_",struct_data(gl_axis).label))     = wvar_reduction;
                if typefit == "gp"
                    Glitch_cat(idx_glitch2).(strcat("thetaP_",struct_data(gl_axis).label))       = theta(3);  
                    Glitch_cat(idx_glitch2).(strcat("thetaPp_",struct_data(gl_axis).label))       = theta(4);  
                    Glitch_cat(idx_glitch2).(strcat("delayP_",struct_data(gl_axis).label))       = theta(4)/theta(3);
                end

                if wvar_reduction <= var_thres
                    Glitch_cat(idx_glitch2).(strcat("isglitch_",struct_data(gl_axis).label)) = 1;
                    Glitch_cat(idx_glitch2).comment = "PolyGlitch";
                    struct_data(gl_axis).drif(idx_in_ts-ND2:idx_in_ts+NDS) = drift;
                    struct_data(gl_axis).data_dg(idx_in_ts-ND2:idx_in_ts+NDS) = residual;
                    struct_data(gl_axis).synth_glitch(idx_in_ts-ND2:idx_in_ts+NDS) = yhat;      
                else
                    Glitch_cat(idx_glitch2).(strcat("isglitch_",struct_data(gl_axis).label)) = 0;
                    Glitch_cat(idx_glitch2).comment = "PolyGlitch";
                    struct_data(gl_axis).drif(idx_in_ts-ND2:idx_in_ts+NDS) = drift;
                end
            else
                disp(["Glitch # ",idx_in_ts, " rejected, out of bounds."])
            end % end sur la longueur du glitch
        end
        iter_idx = iter_idx+length(clusters_of_glitches{clust_idx});
        %disp(strcat("Next start index: ", num2str(iter_idx)))
    end
end