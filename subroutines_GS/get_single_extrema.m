function [single_extrema_list] = get_single_extrema(input_vec, radius)
    % Function to get the list of the extrema with no extrema around in a
    % given radius. 
    %  
    % -----
    % INPUT:
    %      @input_vec : one dimensional array with the indexes of the extrema
    %      @radius : it's a number to identify single glitches
    %
    % ------
    % OUTPUT:
    %      @single_extrema_list : list of index of extrema with no
    %      neighbourgh less than the radius
    %
    % 
    single_extrema_list = [];
    input_vec = [0; input_vec; input_vec(end)+2*radius];
    for i=2:length(input_vec)-1
        d_l = abs(input_vec(i-1)-input_vec(i));
        d_r = abs(input_vec(i)-input_vec(i+1));
        if d_l >= radius && d_r >=radius
            single_extrema_list = [single_extrema_list; input_vec(i)];
        end
    end
end