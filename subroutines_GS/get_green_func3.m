%% Definition of the function get_green_func.
function [green_fn_list,  time_output] = get_green_func3(sensor, mode, resp, PFO_list, dt, ND2, NDS)
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
    %              2.0 - 2021/06/21 - 
    % ----------
    % Parameters:
    %       @sensor       : VBB or SP
    %       @mode        : "VEL" or "POS"
    %       @resp          : array with the path of the dataless (contaning the transfert functions)
    %       @PFO_list      : array of cell with FIR coeffs with gain 2, 4 and 5
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
    in_sps        = 20;        %The input sampling rate is for the moment supposed to be 20Hz  
    NS = 4096;
    dti=1/in_sps;              % Output sampling rate
    NSi=10*NS;

    % Read of the transfert functions
    % Transfert functions are saved in the directory "metadata"
     
    [ugain_vbbvel, uzero_init_vbbvel, upole_init_vbbvel, unitsvbbvel, uA0_vbbvel]=read_resp_v2(resp.vbb_vel_hg);
    [ugain_vbbpos, uzero_init_vbbpos, upole_init_vbbpos, units_vbbpos, uA0_vbbpos]=read_resp_v2(resp.vbb_pos_hg);

     
     green_fn_list = [];                                 % Output var with the list of the green functions
     time_20sps=[0:NSi-1]*dti;                    % Time vector 
     
     % Modify pole list if VEL or POS
     upole_vbbvel = [0 upole_init_vbbvel];
     upole_vbbpos = upole_init_vbbpos;
  
   
   %% Estimate the first Green function  for VEL data
   
    %[ugain, uzero_init, upole_init, units, uA0]=read_resp_v2(resp.vbb_vel_hg);
    %sys = zpk(uzero_init, upole, ugain*uA0);  
   
   sysvel = zpk(uzero_init_vbbvel, upole_vbbvel, ugain_vbbvel*uA0_vbbvel);
   sysvel.TimeUnit = "seconds";
  
   figure(3000)
   bode(sysvel);
   grid on
   
   
   
   % Apply a step function to the TF to create a Glitch @ 20sps
   green_20sps_vel= step(sysvel, time_20sps);
   
  
   %% Estimate the first Green function  for POS data
    
   syspos = zpk(uzero_init_vbbpos, upole_vbbpos, ugain_vbbpos*uA0_vbbpos);
   syspos.TimeUnit = "seconds";
  
   figure(3001)
   bode(syspos);
   grid on
   
   % Apply a step function to the TF to create a Glitch @ 20sps
   green_20sps_pos= step(syspos, time_20sps);
   
   
   %% Initiate the time vector and shift the data
   % It is first center than shifted to slide the glitch at the beginning
   % It must be done for both VEL and POS even if only POS is resquested
   % because the code is shifting the green function of the POS according
   % to the max of the green function of the VEL
   
   NF=floor(NSi/2+1);                                                           % middle point of the vector
   green_20sps_shift_vel = [zeros(1, NF), green_20sps_vel'];   % Shift the vector to center the max of the Green func
   green_20sps_shift_vel = green_20sps_shift_vel(1:NSi);
   
   green_20sps_shift_pos = [zeros(1, NF), green_20sps_pos']; % Shift the vector to center the max of the Green func
   green_20sps_shift_pos = green_20sps_shift_pos(1:NSi);
   
   % Same time vector for POS and VEL
   time_20sps_shift   = [zeros(1, NF), time_20sps];
   time_20sps_shift  = time_20sps_shift(1:NSi);
   
   % Append in a list
   green_fn_list_20_vel = [green_20sps_shift_vel];
   green_fn_list_20_pos = [green_20sps_shift_pos];
   
   % Calculation of the derivates of G
   % To avoid causality issues, calculation is made using the TF derivation
   % f^{(k)}_m = (imw)^kf_m$ avec $f_m = \frac{1}{T}\int_0^T e^{imwt}f(t)dt  
   for deriv=1:3
       omega(1:NF)=[0:NF-1]*1i/NS/dti*2*pi;
       % Fourier Transform applied to VEL 
       tf_vel=fft(green_20sps_shift_vel);  
       tf_vel(1:NF)=tf_vel(1:NF).*(omega(1:NF));
       green_20sps_shift_vel=ifft(tf_vel,'symmetric');
       green_fn_list_20_vel = [green_fn_list_20_vel; green_20sps_shift_vel];
       
       % Fourier Transform applied to POS 
       tf_pos=fft(green_20sps_shift_pos);  
       tf_pos(1:NF)=tf_pos(1:NF).*(omega(1:NF));
       green_20sps_shift_pos=ifft(tf_pos,'symmetric');
       green_fn_list_20_pos = [green_fn_list_20_pos; green_20sps_shift_pos];  
   end

    [numRows, ~] = size(green_fn_list_20_vel);
   
    % Calculus of the decimation factors between the original frequency to
    % the output frequency.
    % Here again, decimation is applied on both green's function on VEL and
    % POS 
    FIR_coeff_list = [2, 4, 5];         % Decimation filters available
   
    if out_sps < in_sps
       fir2apply = FIRdecomp(out_sps, in_sps, FIR_coeff_list);
       fir2apply = sort(fir2apply, 'descend');
       if length(fir2apply)>=1 
           for i=1:numRows
               green_decim_f_vel = green_fn_list_20_vel(i, :);
               green_decim_f_pos = green_fn_list_20_pos(i, :);
                    % Decimate data by applying FIR
                    for fc=1:length(fir2apply)
                        if fir2apply(fc) == 2
                            PFO2apply = PFO_list{1};
                        elseif fir2apply(fc) == 4
                            PFO2apply = PFO_list{2};
                        elseif fir2apply(fc) == 5
                            PFO2apply = PFO_list{3};
                        end
                        green_decim_f_vel  =  decimeFPGA(green_decim_f_vel, fir2apply(fc), ... 
                                                                        PFO2apply, sum(PFO2apply));

                        green_decim_f_pos  =  decimeFPGA(green_decim_f_pos, fir2apply(fc), ... 
                                                                        PFO2apply, sum(PFO2apply));
                    end
                disp(strcat("length of POS @ i= ", num2str(i), " ->", num2str(length(green_decim_f_pos))))
                disp(strcat("length of VEL @ i= ", num2str(i), " ->", num2str(length(green_decim_f_vel))))
                % Step to shift all the function the same way
                
                [~, idx_max_2sps] = max(abs(green_decim_f_vel));
                
                disp(strcat("Index of max: ", num2str(idx_max_2sps)))
                
                if i==1
                    idx_ref = idx_max_2sps;
                end

                if mode == "VEL"
                    green = green_decim_f_vel(idx_ref-ND2: idx_ref+NDS);
                elseif mode == "POS"
                    green = green_decim_f_pos(idx_ref-ND2: idx_ref+NDS);
                end

                 disp(strcat("length of green : ", num2str(length(green))))
                
                % Append the function to the output list                                                                                      
                green_fn_list = [green_fn_list ; green];
            end  % loop over function 
            time_output= [0:length(green(1,:))-1]*dt;
       end 
    elseif out_sps == in_sps
        
        for i=1:numRows
               green_decim_f_vel = green_fn_list_20_vel(i, :);
               green_decim_f_pos = green_fn_list_20_pos(i, :);
       
               [~, idx_max_2sps] = max(abs(green_decim_f_vel));
                if i==1
                    idx_ref = idx_max_2sps;
                end

                if mode == "VEL"
                    green = green_decim_f_vel(idx_ref-ND2: idx_ref+NDS);
                elseif mode == "POS"
                    green = green_decim_f_pos(idx_ref-ND2: idx_ref+NDS);
                end
               
                green_fn_list = [green_fn_list ; green];
                time_output= [0:length(green(1,:))-1]*dt;
        end
        
    else
        disp("ERROR: Output frequency is higher or equal to the input frequency.")
    end % End test over sps
end