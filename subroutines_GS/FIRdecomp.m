function [fir2apply] = FIRdecomp(out_sampling, input_sampling, FIR_coeff_list) 
    % Function to find the FIR factors to apply in case of decimation from
    % an input signal with a input_sampling to an output signal with an
    % output_sampling using a list of decimation factors
    % Example, for InSight : FIR_coeff_list = [2, 4, 5]
    %
    % author  : Grégory Sainton (IPGP) on behalf NASA InSight Collaboration
    % date     : 2020/12/04
    % version : 1.0
    %
    % --------
    % INPUT:
    %       @output_sampling : float - final sampling of the signal after decimation
    %       @input_sampling   : float - initial sampling of the signal.
    %       @FIR_coeff_list       : array of integer with the list of decimation filter available  
    %--------- 
    % OUPUT:
    %       @fir2apply : array of coefficients to apply to decimate from input to output sampling 
    %
    % EXAMPLE: 
    %       >> fir2apply = FIRdecomp(20, 2, [2,4,5])
    %       fir2apply =
    %                   2     5
    %
    assert(input_sampling>out_sampling, "Final sampling can't be higher of input signal after decimation.")
    
    if ~isempty(FIR_coeff_list) 
       FIR_coeff_list = sort(FIR_coeff_list);                                  % Sort to ascending order
       fir2apply = [];
       i=1;
       ratio = input_sampling/out_sampling;
       
       while i<length(FIR_coeff_list) || ratio~=1
           tmp = ratio/FIR_coeff_list(i);
           [frac, ~] = modf(tmp);
           if  frac == 0
                fir2apply= [fir2apply, FIR_coeff_list(i)];
                ratio = tmp;
            else
                i=i+1;     
            end
        end
    end
end