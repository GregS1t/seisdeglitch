function  struct_data = find_local_maxima2(struct_data, prom_factor, cut_fact, plot_extrema, verbose)
    % Function to find the local maxima for each the 3 input vectors.
    % For each time series, the function is looking for the local maxima and the local minima 
    % The function create a unique vector with all the maxima and minima
    % for the 3 axis. 
    %
    % INPUT:
    % ------
    %       @struct_data    : array with the 3 axis signa data
    %       @prom_factor  : prominence factor
    %       @plot_extrema : 0/1 - 1 if user wants to plot the signal with
    %       their maxima
    %       @verbose         : 0/1 - 1 if user wants some explanations
    %
    % OUTPUT:
    % -------
    %        @struct_data : struct type. It contains all the input data 
    %                              + the peaks list for each axis and the 
    %                                  number of max per axis.   
    %
    
    
    % Cut to eliminate "rebound" of the glitches considered as
    % extrema by the functions islocalmin and islocalmax
    epsG = 5;
    gapG = 20;
    
    
    
    if verbose == 1
        disp("In find_local_maxima function : ")
        disp("_______________________________")
    end
       
    for i_ax=1:3                                                                   %loop over axes
        
        prominence = prom_factor*std(struct_data( i_ax).data);            % Relative amplitude of a maxima according to the 2 neighbours
        mean_sig = mean(struct_data( i_ax).data);
        std_sig = std(struct_data( i_ax).data);
       
        
        %First cut
        %TF_localmin = islocalmin(struct_data( i_ax).data, 'SamplePoints', struct_data( i_ax).time, "MinProminence", prominence);
        %TF_localmax = islocalmax(struct_data( i_ax).data, 'SamplePoints', struct_data( i_ax).time, "MinProminence", prominence);
        TF_localmin = islocalmin(struct_data( i_ax).data, 'SamplePoints', struct_data( i_ax).time);
        TF_localmax = islocalmax(struct_data( i_ax).data, 'SamplePoints', struct_data( i_ax).time);
    
        IdxlocalMax = find(TF_localmax == 1);
        IdxlocalMin = find(TF_localmin == 1);
        
        %Second cut on the amplitude to remove 
         idx_sub_max = find((struct_data(i_ax).data(IdxlocalMax) >= mean_sig + cut_fact*std_sig));
         idx_sub_min = find((struct_data(i_ax).data(IdxlocalMin) <= mean_sig - cut_fact*std_sig));
         
         IdxlocalMax = IdxlocalMax(idx_sub_max);
         IdxlocalMin = IdxlocalMin(idx_sub_min);
        
        mergeAxis = union(IdxlocalMax , IdxlocalMin);
     
       for i=1:length(mergeAxis)-1  
            if (abs(mergeAxis(i)-mergeAxis(i+1)) >= gapG - epsG) && ...                                % select extrema inside range gapG +/- epsG 
               (abs(mergeAxis(i)-mergeAxis(i+1))< gapG + epsG)                   
                if struct_data(i_ax).data(mergeAxis(i))*struct_data(i_ax).data(mergeAxis(i+1))<0  % rebound must have opposite sign with 
                                                                                                                                      % primary signal
                    mergeAxis(i+1)=NaN;                                                                                 % set to NaN because convienient to remove
                end
            end
       end
       mergeAxis = mergeAxis(~isnan(mergeAxis));                                                              % Finally remove the NaN created above

       struct_data(i_ax).peaks = mergeAxis;
        
       struct_data(i_ax).nbpeaks = length(struct_data(i_ax).peaks);
      
        % Plot option
        if plot_extrema ==1
            figure(100+10*i_ax)
            %plot(struct_data( i_ax).time, struct_data(i_ax).data, 'DisplayName', struct_data(i_ax).label)
            plot(struct_data(i_ax).data, 'DisplayName', struct_data(i_ax).label)
            hold on
            
            %plot(struct_data( i_ax).time(IdxlocalMin), struct_data( i_ax).data(IdxlocalMin), 'gO',  'DisplayName', "local minima")
            %plot(struct_data( i_ax).time(mergeAxis), struct_data( i_ax).data(mergeAxis), 'gO',  'DisplayName', "Extrema")
            plot(mergeAxis, struct_data( i_ax).data(mergeAxis), 'gO',  'DisplayName', "Extrema")
            title(["Location of the extrema on axis ", i_ax])
            grid
            legend
            
        end
    end
    clear tmp_vec
    if verbose == 1
        disp("End of execution find_local_maxima")
    end
end