function  struct_data = find_local_maxima(struct_data, prom_factor, plot_extrema, verbose)
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
        
    %loop over axes
    for i_ax=1:3
        
        prominence = prom_factor*std(struct_data( i_ax).data);            % Relative amplitude of a maxima according to the 2 neighbours
      
        TF_localmin = islocalmin(struct_data( i_ax).data, 'SamplePoints', struct_data( i_ax).time, "MinProminence", prominence);
        TF_localmax = islocalmax(struct_data( i_ax).data, 'SamplePoints', struct_data( i_ax).time, "MinProminence", prominence);

        NblocalMax = find(TF_localmax == 1);
        NblocalMin = find(TF_localmin == 1);
        
        % merge the min and max in a same array
        struct_data(i_ax).peaks = union(NblocalMax , NblocalMin);
        struct_data(i_ax).nbpeaks = length(struct_data(i_ax).peaks);
        
        if verbose == 1
            disp(["Axis: " struct_data(i_ax).label])
            disp(["   -> Nb local Min: " length(NblocalMin) " with prominence " prominence])
            disp(["   -> Nb local Max: " length(NblocalMax) " with prominence " prominence])
            disp(["   -> A  total of " num2str(struct_data(i_ax).nbpeaks) " local maxima with a prominence = " num2str(prominence)])
        end

        % Plot option
        if plot_extrema ==1
            figure(100+10*i_ax)
            subplot (3,1,1)
            plot(struct_data( i_ax).time, struct_data(i_ax).data, 'DisplayName', struct_data(i_ax).label)
            hold on
            plot(struct_data( i_ax).time(TF_localmin), struct_data( i_ax).data(TF_localmin), 'gO',  'DisplayName', "local minima")
            title("Look for local minima")
            grid
            legend

            subplot (3,1,2)
            plot(struct_data( i_ax).time, struct_data(i_ax).data, 'DisplayName', struct_data(i_ax).label)
            hold on
            plot(struct_data(i_ax).time(TF_localmax), struct_data( i_ax).data(TF_localmax),'b*',  'DisplayName', "local maxima")
            title("Look for local maxima")
            grid
            legend

            subplot (3,1,3)
            tmp_vec = abs(struct_data( i_ax).data);
            plot(struct_data(i_ax).time, tmp_vec, 'k', 'DisplayName', strcat("abs(",struct_data(i_ax).label, ")"))
            hold on
            plot(struct_data( i_ax).time(TF_localmin), tmp_vec(TF_localmin), 'gO',  'DisplayName', "local minima")
            plot(struct_data( i_ax).time(TF_localmax), tmp_vec(TF_localmax), 'b*',  'DisplayName', "local maxima"')
            title("Absolute value of the signal with the local maxima")
            grid on
            legend 
            
        end
    end
    clear tmp_vec
    if verbose == 1
        disp("End of execution find_local_maxima")
    end
end