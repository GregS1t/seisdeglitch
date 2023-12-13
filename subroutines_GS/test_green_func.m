clear all
close all 
clc

%% Parameters

sday = 502;
dt     = 1;

dirDROP = '/Users/Greg/Dropbox (IPGP)/';
dir_FIR         = [dirDROP  'DeglitchPkg_PL_GS/FIR/'];
PFO5            = load([dir_FIR 'PFO_div5.txt']);      % FIR coefs to decimate by 5
PFO4            = load([dir_FIR 'PFO_div4.txt']);      % FIR coefs to decimate by 5
PFO2            = load([dir_FIR 'PFO_div2.txt']);      % FIR coefs to decimate by 2
PFO_list        = {PFO2, PFO4, PFO5};

% Dir to the metadata
dir_resp        = [dirDROP  'DeglitchPkg_PL_GS/metadata/'];
resp.vbb_vel_hg = [dir_resp  'RESP.XB.ELYSE.02.BHU_new'];
resp.vbb_pos_hg = [dir_resp  'RESP.XB.ELYSE.00.LMU'];
resp.sp_hg      = [dir_resp  'RESP.XB.ELYSE.65.EHU'];


ND2 =30*floor(0.5/dt)               % Number of points before the max amplitude of the GF 
NDS = floor(6*ND2)                    % Number of points after the max (Philippe's value is %NDsum=6*ND2;)
                               
sensor = "VBB"; 
mode   = "VEL";

verbose = 1;
%%


%% Green function generation
% Call get_green_func function to create synthetic functions
% See inside the function for more details
plot_green = 1;

[green_fn_list_vel, green_timevec_vel] = get_green_func3(sensor, "VEL", resp, PFO_list, dt, ... 
                                                    ND2, NDS, plot_green, verbose);

[green_fn_list_pos, green_timevec_pos] = get_green_func3(sensor,"POS", resp, PFO_list, dt, ...
                                                    ND2, NDS, plot_green, verbose);

     
figure(1)
for cpt=1:4
    subplot(4,1, cpt)
    plot(green_timevec_vel, green_fn_list_vel(cpt, :), "DisplayName", "VEL data")
    ylabel("VEL in DU")
    yyaxis right
    plot(green_timevec_pos, green_fn_list_pos(cpt, :), "DisplayName", "POS data")
    ylabel("POS in DU")
    legend
    grid on
end
sgtitle(strcat("VEL vs POS Green fn and derivatives @ ", num2str(1/dt), " sps"))