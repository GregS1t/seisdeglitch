
r2 = -5 + (5+5)*rand(10,1);

[idxEx, ampEx] = find_extremum(r2)

function [idxEx, ampEx] = find_extremum(inputdata)
    % Function to find extremum of an input signal.
    % It can be either a maximum or a minimum: the function return  
    % the maximum of the both absolute value of min and max
    % -------
    % INPUT
    %   @inputdata: input 1D array 
    % 
    % -------
    % OUTPUT
    %   @idxEx   : index of the extremum value
    %   @ampEx : value of the extremum 
    
    [minSig, minIdx]  = min(inputdata)
    [maxSig, maxIdx] = max(inputdata)
    
    if abs(minSig) > abs(maxSig)
        ampEx = minSig ; idxEx = minIdx;
    else
        ampEx = maxSig ;  idxEx = maxIdx;
    end
end
