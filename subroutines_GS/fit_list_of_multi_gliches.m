function [struct_data, Glitch_cat] =  fit_list_of_multi_gliches(struct_data, clusters_of_glitches, ...
                                                        max_glitch_per_axis, green_fn_list, dt, typefit, lagrangemult, ND2, NDS, ...
                                                        Glitch_cat, verbose)
% Function to fit a list of cluster of gliches 
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
    cpt =1;
    for i_ax=1:3 
        %disp(strcat("i_ax = ", num2str(i_ax)))
        iter_idx = length_results_fit_init;
        
        for clust_idx=2:length(clusters_of_glitches) % Loop over all the clusters
            %disp( clusters_of_glitches{clust_idx})
            
            [theta, yhat] = fit_multi_glitches_lagrange(struct_data(i_ax).data_dg, ...
                                    clusters_of_glitches{clust_idx}, green_fn_list, dt, ND2, NDS, typefit, lagrangemult);
            
           index_min = min(clusters_of_glitches{clust_idx});
           index_max = max(clusters_of_glitches{clust_idx});
           y = struct_data(i_ax).data_dg(index_min-ND2:index_max+NDS);
           struct_data(i_ax).data_dg(index_min-ND2:index_max+NDS) = y - yhat;
           struct_data(i_ax).synth_glitch(index_min-ND2:index_max+NDS) = yhat;                       

           for gl=0:length(clusters_of_glitches{clust_idx})-1         % Need to decompose each cluster fit
                subtheta = [theta(1+4*gl),  theta(2+4*gl), theta(3+4*gl), ...
                                     theta(4+4*gl), theta(end-1), theta(end)];
                iter_idx = iter_idx+1;
        
%                 x5 = ones(ND2+NDS+1, 1);
%                 x6 = [0:ND2+NDS]'*1.0;
%                 
%                 X = [green_fn_list', x5, x6];
%                 subyhat = X*subtheta';      
%                 [~, ampEx] = find_extremum(subyhat); 
                idxsub = clusters_of_glitches{clust_idx}(gl+1); 
                subyhat = struct_data(i_ax).synth_glitch(idxsub-ND2:idxsub+NDS);
                [~, ampEx] = find_extremum(subyhat); 
                
                ysub = struct_data(i_ax).data_dg(idxsub-ND2:idxsub+NDS);
                
                Glitch_cat(iter_idx).index            = clusters_of_glitches{clust_idx}(gl+1);
                Glitch_cat(iter_idx).time              = struct_data(i_ax).time(clusters_of_glitches{clust_idx}(gl+1));
                Glitch_cat(iter_idx).lmst              = struct_data(i_ax).lmst(clusters_of_glitches{clust_idx}(gl+1));
                Glitch_cat(iter_idx).(strcat("data_", struct_data(i_ax).label))      = struct_data(i_ax).data(clusters_of_glitches{clust_idx}(gl+1));
                Glitch_cat(iter_idx).(strcat("maxAmpFit_",struct_data(i_ax).label))    = ampEx;
                Glitch_cat(iter_idx).(strcat("thetaG_",struct_data(i_ax).label))    = subtheta(1);
                Glitch_cat(iter_idx).(strcat("delayG_",struct_data(i_ax).label))  = subtheta(2)/subtheta(1); 
                % Restimate variance of the glitch only
                Glitch_cat(iter_idx).(strcat("var_",struct_data(i_ax).label)) = var(ysub);  
                Glitch_cat(iter_idx).(strcat("var_res_",struct_data(i_ax).label)) = var(ysub-subyhat);
                Glitch_cat(iter_idx).(strcat("var_ratio_",struct_data(i_ax).label)) = var(ysub-subyhat)/var(ysub);
                % Keep the variance of the whole fit under the 
                Glitch_cat(iter_idx).(strcat("var_tot_",struct_data(i_ax).label)) = var(y);  
                Glitch_cat(iter_idx).(strcat("var_tot_res_",struct_data(i_ax).label)) = var(y-yhat); 
                Glitch_cat(iter_idx).(strcat("var_tot_ratio_",struct_data(i_ax).label)) = var(y-yhat)/var(y);
                if typefit == "gp"
                    Glitch_cat(iter_idx).(strcat("thetaP_",struct_data(i_ax).label))    = subtheta(3);
                    Glitch_cat(iter_idx).(strcat("delayP_",struct_data(i_ax).label))  = subtheta(4)/subtheta(3);  
                    
                 end
                 cpt=cpt+1;
                 Glitch_cat(iter_idx).comment = strcat("MultiGlitch_", num2str(clust_idx));
                 
           end
        end
        
    end
    if verbose == 1
        disp(["Nombre de glitchs multiple: ", (cpt-1)/3])
    end
end