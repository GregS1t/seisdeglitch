function clusters_of_glitches = clusterize_extrema(glitch_list, radius_min, radius_max)
    % Function to clusterize a 1D list of extrema using a minimum 
    % radius and a maximum radius
    % 
    % --------
    % INPUT:
    %       @glitch_list    : array of glitch indexes to cluserize
    %       @radius_min : Integer - Minimum radius between indexes
    %       @radius_max: Integer - Maximum radius between indexes
    %
    % --------
    % OUTPUT:
    %        @clusters_of_glitches 
    
    % calculate the matrix distance of every point to each others
    dist = pdist2(glitch_list, glitch_list);
    
    % select the submatrix corresponding to the distance interval [radius_min, radius_max]
    B = (dist <=radius_max & dist>radius_min);
    
    % Keep the upper triangle and remove the diagonal
    U = triu(B); V = U.*~eye(size(U));

    % Get back to the index coordinates 
    [row,col]=find(V==1);
    idx_r = glitch_list(row);
    idx_c = glitch_list(col);

    stack_rc = union(idx_c, idx_r);

    st = [];

    diff_gl = diff(stack_rc);
    cell=1;
    % Clusterize the glitches
    % Go through each couple of index
    %   1> First one is save in st
    %   2> Second couple -> check is one of the two index is already in the st
    %          - if yes, add it to the current st, reduce st to avoid
    %               multiple times the same index and and save this array
    %               in the current cell clusters_of_glitches
    %          - if not, create a new st and save it in the next cell clusters_of_glitches
    for i=1:length(diff_gl)
        dist = diff_gl(i);
        if dist<= radius_max
            if ismember(stack_rc(i), st) ==1 || ismember(stack_rc(i+1), st) ==1
                st = [st, stack_rc(i), stack_rc(i+1)];
                st =unique(st);
                clusters_of_glitches{cell} = st;
            else
                cell = cell+1;
                st = [stack_rc(i), stack_rc(i+1)];
                clusters_of_glitches{cell} = st;
            end  
        end
    end
end