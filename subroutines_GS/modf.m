
function [frac, integ] = modf(number)
    % modf() returns the fractional and integer parts of number 
    % in a two-item array
    % based on python math.modf(x)
    
    integ=floor(number);
    frac=number-integ;
end
