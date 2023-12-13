function reduced_glitch_list = merge_close_glitches(input_vec, threshold)
    % Function to merge the closest extrema in a single one
    % It take the center of mass of the extrema positions (and the floor to get an integer) 
    % according to the input threshold.
    % ------
    % INPUT:
    %      @input_vec : one dimensional array with the indexes of the
    %      extrema
    %      @threshold : it's a number to group the extrema into a single one
    % 
    %-----
    % OUTPUT:
    %      @reduced_glitch_list : final list of extrema after merging of
    %      the extrema
    %
 
    sdiff = diff([1; input_vec]);

    % creat mask with threshold
    mask = abs(sdiff(:)) < threshold;

    % assign indices to "similar" values
    subs = cumsum(~mask);
    if subs(1) == 0
        subs = subs + 1;
    end
        
    % calculate the means of similar values
    xmean = accumarray(subs, input_vec, [], @mean);
    
    % create a new vector
    reduced_glitch_list = floor(xmean(subs));
    reduced_glitch_list = sortrows(reduced_glitch_list);
    reduced_glitch_list = unique(reduced_glitch_list);
end