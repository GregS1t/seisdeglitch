%% Definition of the function get_green_func.
function green_fn_list = get_green_func2(sensor, mode, resp, PFO_list, dt, tb_peak, ta_peak, plot_green, verbose)
    % Function which estimate the green functions of a step function
    % 
    % We are assuming that a glitch can be modeled as a step in
    % acceleration. 
    % 1> Creation of a transfert function (TF )using zpk matlab function : 
    %                https://www.mathworks.com/help/control/ref/zpk.html
    % 2> Apply this TF to a step function
    %
    % 3> Depending on the content of type_Green the zero list is modified
    %
    % author  : Grégory Sainton (IPGP) on behalf NASA InSight Collaboration
    % version : 1.1 - 2020/12/04 - update with variable output sampling.
    %              1.0 - 2020/11/15 - initial version. 
    % ----------
    % Parameters:
    %       @sensor       : VBB or SP
    %       @mode        : "VEL" or "POS"
    %       @resp          : array with the path of the dataless (contaning the transfert functions)
    %       @PFO_list      : array FIR coeff with gain 2, 4 and 5
    %       @dt              : sampling rate of the signal 
    %       @ND2           : # of point before max of the GF
    %       @NDS           : # of point after the max of the GF     
    %       @plot_green  : 0 or 1 (1 if you want to plot the green functions)
    %       @verbose      : 0 or 1 (1 if you want informations about transfer functions and so on... )
    % -----------
    % Outputs
    %       @green_fn_list: array of array where each array is a green function 
    %                              according to the list type_green
    % ------------
    % Dependance
    %       decimeFPGA : decimation using the FIR coefs filters
    
    % Parameters
    % We calculte the Green's func, its derivative, the "Precursor" and its
    % derivative
    type_green = ["G", "Gprime", "P", "Pprime"];   
    out_sps      = 1/dt;
    in_sps        = 20;   %The input sampling rate is for the moment supposed to be 20Hz  
    NS = 1024;
    dti=dt/10;              % Output sampling rate
    NSi=10*NS;
    
    ND2 = tb_peak/dt;
    NDS = ta_peak/dt;
    
    
    % Read of the transfert functions
    %        TF are saved in the directory "metadata"
   
     if sensor == "VBB" && mode == "VEL"
            [ugain, uzero_init, upole_init, units, uA0]=read_resp_v2(resp.vbb_vel_hg);
     elseif sensor == "VBB" && mode == "POS"
          [ugain, uzero_init, upole_init, units, uA0]=read_resp_v2(resp.vbb_pos_hg);
     elseif sensor == "SP"
          [ugain, uzero_init, upole_init, units, uA0]=read_resp_v2(resp.sp_hg);
     end
      
    if verbose == 1
        disp(strcat("For the sensor ", sensor, " with mode ", mode, ": "))
        disp(["ugain= " ugain])
        disp(["uA0" uA0])
        disp(["uzero= " uzero_init]) 
        disp(["upole= " upole_init])
        disp(["Units = " units])
    end
     
     green_fn_list = [];                                 % Output var with the list of the green functions
     time_20sps=[0:NSi-1]*dti;                    % Time vector 
    %disp(strcat("length(time_20sps)", num2str(length(time_20sps))))
     % Modify pole list if VEl or POS
     if mode == "VEL"
         upole = [0 upole_init];
     else 
         upole = upole_init;
     end
   
   % Estimate the first Green function  
   % 
   sys = zpk(uzero_init, upole, ugain*uA0);                             % Ca ne marche pas en introduisant le Gain dans zpk
   %sys = zpk(uzero, upole, ugain, dti);                                  % Les gains sont incohérents                            
   sys.TimeUnit = "seconds";
   
    if plot_green == 1 
       figure(300)
       bode(sys)
       h = gcr;
       setoptions(h,'FreqUnits','Hz');
       grid on 
   end
  
   % Apply a step function to the TF to create a Glitch @ 20sps
   green_20sps= step(sys, time_20sps);   
   
   % Center the function 
   NF=floor(NSi/2+1);                                                           % middle point of the vector
   green_20sps_shift = [zeros(1, NF), green_20sps'];              % Shift the vector to center the max of the Green func
   green_20sps_shift = green_20sps_shift(1:NSi);
   
   % Append in a list
   green_fn_list_20 = [green_20sps_shift];
   
   % Calculation of the derivates of G
   % To avoid causality issues, calculation is made using the TF derivation
   % f^{(k)}_m = (imw)^kf_m$ avec $f_m = \frac{1}{T}\int_0^T e^{imwt}f(t)dt  
   for deriv=1:3
       omega(1:NF)=[0:NF-1]*1i/NS/dti*2*pi;
       tf=fft(green_20sps_shift);  
       tf(1:NF)=tf(1:NF).*(omega(1:NF));
       green_20sps_shift=ifft(tf,'symmetric');
       green_fn_list_20 = [green_fn_list_20; green_20sps_shift];
   end
   
  [numRows,~] = size(green_fn_list_20);
  
   if plot_green == 1 
         figure(301);
   end 
   
   % Calculus of the decimation factors between the original frequency to
   % the output frequency.
   FIR_coeff_list = [2, 4, 5];                                                    % Decimation filters available
   fir2apply = FIRdecomp(out_sps, in_sps, FIR_coeff_list); 
   
   green_fn_list_decim = green_fn_list_20;
   
   for i=1:numRows    
       % decimation from 20sps to 2sps
       %    TODO : add parametrization to avoid such hard coding
       if length(fir2apply)>=1 
           for fc=1:length(fir2apply)
                if fir2apply(fc) == 2
                    PFO_filt = PFO_list(1,:);
                elseif fir2apply(fc) == 4
                    PFO_filt = PFO_list(2,:);
                elseif fir2apply(fc) == 5
                    PFO_filt = PFO_list(3,:);
                end
                green_decim  =  decimeFPGA(green_fn_list_decim(i,:), fir2apply(fc), PFO_filt, sum(PFO_filt));
           end
       else
            green_fn_list = green_fn_list_20;
       end
        time=downsample(time_20sps, in_sps/out_sps);
        disp(length(time))

        [~, idx_max_2sps] = max(abs(green_decim));
        if i==1
        idx_ref = idx_max_2sps;
        end
        green = green_decim(idx_ref-ND2: idx_ref+NDS);           % normalisation avec le maximum
                                                                                      % On garde 15 points avant et 135 après
        % Append the function to the output list                                                                                      
        green_fn_list = [green_fn_list ; green];

        if plot_green == 1
            cmap = hsv(numRows); 
            sgtitle("Synthetic functions")
            subplot(numRows, 1, i )
            plot(green, "g", "DisplayName" , type_green(i), 'Color', cmap(i,:))
            ylabel("D.U./m.s^{-2}")
            xlabel("Time in seconds")
            %xlim([0 150])
            legend
            grid on
            title(type_green(i))

        end            % end plot option 
    end                 % loop over function 
 
   clear PFO2 PFO5
end