% NDeriv est toujours égal à 3
% icase = 2
% uz = 
% NF=floor(NS/2+1); 

if icase == 1
    dti=dt;
    %%%%% built Green response from Pole and Zero ( remove the two first null zero)
    green(1:NS)=0.;  % Vecteur de 1024 points
    green(NS/2)=1.;  % Ca c'est un heaviside sur NS/2 
    time=[0:NS-1]*dti;
    
    tf0=fft(green); 
    
    tf=tf0;
    tf1=tf;
    tf2=tf;
    tf3=tf;
    NF=floor(NS/2+1);
    omega(1:NF)=[0:NF-1]*1i/NS/dti*2*pi;

    %%%% remove the two first zero and keep the other ones for Nderiv=0;
    for ii=Nderiv:length(uz)
        tf(1:NF)=tf(1:NF).*(omega(1:NF)-uz(ii));
    end
    %%%% take all poles
    for ii=1:length(up)
        tf(1:NF)=tf(1:NF)./(omega(1:NF)-up(ii));
    end

    %%%%
    ii1=length(up)-1;
    ii2=length(up);
    
    
    for i=1:NF
        tf3(i)=tf(i)*(complex(0.,1.)*imag(up(ii1))/(omega(i)-up(ii1))+complex(0.,1.)*imag(up(ii2))/(omega(i)-up(ii2)) );
        tf2(i)=tf(i)*(real(up(ii1))/(omega(i)-up(ii1))+real(up(ii2))/(omega(i)-up(ii2)) );
    end

    %%%%% get the time function of the Green response as well as time derivativej
    green=ifft(tf,'symmetric');
    tf1(1:NF)=tf(1:NF).*omega(1:NF);
    %tf2(1:NF)=tf1(1:NF).*omega(1:NF);
    tf2(1:NF)=tf0(1:NF);
    green_1=ifft(tf1,'symmetric');
    green_2=ifft(tf2,'symmetric');
    green_3=ifft(tf3,'symmetric');

    figure(30+Nderiv)
    disp(["NDeriv: " Nderiv])
    
    [adelay delay]=max(abs(green_1));
    subplot(1,3,1)
    plot(time-delay*dt, green/adelay, 'k')
    xlabel('Time sec')
    ylabel('green')
    hold on
    subplot(1,3,2)
    plot(time-delay*dt,green_1/adelay,'k')
    xlabel('Time sec')
    ylabel('green_1')
    hold on
    subplot(1,3,3)
    plot(time-delay*dt, green_2/adelay, 'k')
    xlabel('Time sec')
     ylabel('green_2')
    hold on
    title("Green functions")
    
%% icase == 2    
elseif icase == 2
    %%%%% decimate by 10, e.g. 2 sps from 20 sps
    dti=dt/10;
    NSi=10*NS;                              % NS = 1024
    %%%%% built Green response from Pole and Zero ( remove the two first null zero)
    green20sps(1:NSi)=0.;
    green20sps(NSi/2)=1.;
    time20sps=[0:NSi-1]*dti;
    tf020sps=fft(green20sps);
    tf20sps=tf020sps;
    tf120sps=tf20sps;
    tf220sps=tf20sps;
    tf320sps=tf20sps;
    tf420sps=tf20sps;
    tf520sps=tf20sps;
    tf620sps=tf20sps;
    tf720sps=tf20sps;
    tf820sps=tf20sps;
    NF=floor(NSi/2+1);
    omega20sps(1:NF)=[0:NF-1]*1i/NSi/dti*2*pi;
    tf720sps(1:NF)=tf020sps(1:NF).*omega20sps(1:NF);
    tf820sps(1:NF)=tf720sps(1:NF).*omega20sps(1:NF);

    %%%% pass through the response for the impulse

    for ii=1:length(uz)  % 4 "zeros" dans la TF des VBB
        tf020sps(1:NF)=tf020sps(1:NF).*(omega20sps(1:NF)-uz(ii));
    end

    %%%% all pole for TF5 
    for ii=2:length(uz)         % 2 à 4 
        tf520sps(1:NF)=tf520sps(1:NF).*(omega20sps(1:NF)-uz(ii));
    end
    %%%% remove the two first zero and keep the other ones
    for ii=Nderiv:length(uz)  % 3 à 5
        tf20sps(1:NF)=tf20sps(1:NF).*(omega20sps(1:NF)-uz(ii));
    end
    
    %%%% take all poles
    for ii=1:length(up)
        tf520sps(1:NF)=tf520sps(1:NF)./(omega20sps(1:NF)-up(ii));
        tf20sps(1:NF)=tf20sps(1:NF)./(omega20sps(1:NF)-up(ii));
        tf020sps(1:NF)=tf020sps(1:NF)./(omega20sps(1:NF)-up(ii));
    end
    %%%%
    ii1=length(up)-1;
    ii2=length(up);
    for i=1:NF
        tf320sps(i)=tf20sps(i)*(complex(0.,1.)*imag(up(ii1))/(omega20sps(i)-up(ii1))+complex(0.,1.)*imag(up(ii2))/(omega20sps(i)-up(ii2)) );
        tf420sps(i)=tf20sps(i)*(real(up(ii1))/(omega20sps(i)-up(ii1))+real(up(ii2))/(omega20sps(i)-up(ii2)) );
    end

    tf620sps(1:NF)=tf520sps(1:NF).*omega20sps(1:NF);

    %%%%% get the time function of the Green response as well as time derivativej
    green20sps=ifft(tf20sps,'symmetric');
    tf120sps(1:NF)=tf20sps(1:NF).*omega20sps(1:NF);
    green20sps_1=ifft(tf120sps,'symmetric');
    tf220sps(1:NF)=tf120sps(1:NF).*omega20sps(1:NF);
    green20sps_2=ifft(tf220sps,'symmetric');
    green20sps_3=ifft(tf320sps,'symmetric');
    green20sps_4=ifft(tf420sps,'symmetric');
    green20sps_5=ifft(tf520sps,'symmetric');
    green20sps_6=ifft(tf620sps,'symmetric');
    green20sps_7=ifft(tf020sps,'symmetric');
    green20sps_8=ifft(tf720sps,'symmetric');
    green20sps_9=ifft(tf820sps,'symmetric');

    %%%%% integrate the green response
    green20sps_10(1)=0.;
    for i=2:length(green20sps)
    green20sps_10(i)=green20sps_10(i-1)+green20sps(i);
    end


    %%%% decimate by 10
    %%%%% decimate with PFO filter %%%%%%%%%%%%%%%
    PFO5=load('PFO_div5.txt');
    PFO2=load('PFO_div2.txt');
    gain5=sum(PFO5);
    gain2=sum(PFO2);
    N2=length(PFO2);
    N5=length(PFO5);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps=decimeFPGA(green20sps,5,PFO5,gain5);
    green2sps=decimeFPGA(green4sps,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_1=decimeFPGA(green20sps_1,5,PFO5,gain5);
    green2sps_1=decimeFPGA(green4sps_1,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_2=decimeFPGA(green20sps_2,5,PFO5,gain5);
    green2sps_2=decimeFPGA(green4sps_2,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_3=decimeFPGA(green20sps_3,5,PFO5,gain5);
    green2sps_3=decimeFPGA(green4sps_3,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_4=decimeFPGA(green20sps_4,5,PFO5,gain5);
    green2sps_4=decimeFPGA(green4sps_4,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_5=decimeFPGA(green20sps_5,5,PFO5,gain5);
    green2sps_5=decimeFPGA(green4sps_5,2,PFO2, gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_6=decimeFPGA(green20sps_6,5,PFO5,gain5);
    green2sps_6=decimeFPGA(green4sps_6,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_7=decimeFPGA(green20sps_7,5,PFO5,gain5);
    green2sps_7=decimeFPGA(green4sps_7,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_8=decimeFPGA(green20sps_8,5,PFO5,gain5);
    green2sps_8=decimeFPGA(green4sps_8,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_9=decimeFPGA(green20sps_9,5,PFO5,gain5);
    green2sps_9=decimeFPGA(green4sps_9,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green4sps_10=decimeFPGA(green20sps_10,5,PFO5,gain5);
    green2sps_10=decimeFPGA(green4sps_10,2,PFO2,gain2);

    %%%%% 1 et 4 sont lisses
    time=downsample(time20sps, 10);
    [adelay delay]=max(abs(green2sps));
    [agdelay gdelay]=max(abs(green2sps_5));
    [ag2delay g2delay]=max(abs(green2sps_6));
    [ag3delay g3delay]=max(abs(green2sps_6));
    [ag4delay g4delay]=max(abs(green2sps_7));
    figure(30+Nderiv)
    subplot(1,4,1)
    plot(time-delay*dt,green2sps/adelay,'b')
    hold on
    plot(time-delay*dt,green2sps_5/agdelay,'k') 
    subplot(1,4,2)
    plot(time-delay*dt,green2sps_1/adelay,'b')
    hold on
    subplot(1,4,3)
    plot(time-delay*dt,green2sps_2/adelay,'b')
    hold on
    subplot(1,4,4)
    plot(time-delay*dt,green2sps_3/adelay,'k') 
    hold on
    %%%% weight
    dtweight=4.;
    ntweight=floor(dtweight/dt);
    weight(1:length(green2sps_5))=1.;
    weight(gdelay-ntweight:gdelay+ntweight)=100.;


    green=green2sps;
    green_1=green2sps_1;
    Ng2=length(green2sps_2);
    %green_2(1:Ng2-1)=green2sps_2(2:Ng2); %green_2(Ng2)=green_2(Ng2-1);
    green_2=green2sps_2;
    green_3=green2sps_3;
    green_4=green2sps_4;
    green_5=green2sps_5*adelay/agdelay;
    green_6=green2sps_6*adelay/ag2delay;
    green_7=green2sps_7*adelay/ag3delay;
    green_8=green2sps_8*adelay/ag4delay;
    green_9=green2sps_9;
    green_10=green2sps_10;
%% icase == 3
elseif icase == 3

    %%%% decimate by 2, e.g. 20 sps for 10 sps

    dti=dt/2;
    NSi=2*NS;
    %%%%% built Green response from Pole and Zero ( remove the two first null zero)
    
    clear green20sps green20sps_1 green20sps_2 tf20sps tf20sps_2 tf120sps tf220sps omega20sps
    green20sps(1:NSi)=0.;
    green20sps(NSi/2)=1.;
    time20sps=[0:NSi-1]*dti;
    tf20sps=fft(green20sps);
    tf20sps_2=tf20sps;
    tf120sps=tf20sps;
    tf220sps=tf20sps;
    NF=floor(NSi/2+1);
    omega20sps(1:NF)=[0:NF-1]*1i/NSi/dti*2*pi;

    %%%% remove the two first zero and keep the other ones
    for ii=Nderiv:length(uz)
        tf20sps(1:NF)=tf20sps(1:NF).*(omega20sps(1:NF)-uz(ii));
    end
    %%%% take all poles
    for ii=1:length(up)
        tf20sps(1:NF)=tf20sps(1:NF)./(omega20sps(1:NF)-up(ii));
    end
    %%%%
    %%%%
    ii1=length(up)-1;
    ii2=length(up);
    for i=1:NF
    tf20sps_2(i)=tf20sps(i)*(-complex(0.,1.)/(omega20sps(i)-up(ii1))+complex(0.,1.)/(omega20sps(i)-up(ii2)) );
    end

    %%%%% get the time function of the Green response as well as time derivativej
    green20sps=ifft(tf20sps,'symmetric');
    tf120sps(1:NF)=tf20sps(1:NF).*omega20sps(1:NF);
    green20sps_1=ifft(tf120sps,'symmetric');
    green20sps_2=ifft(tf20sps_2,'symmetric');

    %%%% decimate by 2
    %%%%% decimate with PFO filter %%%%%%%%%%%%%%%
    PFO2=load('PFO_div2.txt');
    gain2=sum(PFO2);
    N2=length(PFO2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green10sps=decimeFPGA(green20sps,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green10sps_1=decimeFPGA(green20sps_1,2,PFO2,gain2);
    %%%%  simulate the FPGA             %%%%%%%%%
    green10sps_2=decimeFPGA(green20sps_2,2,PFO2,gain2);


    time=downsample(time20sps,2);
    NN=min([length(time),length(green10sps_1),length(green10sps_2)]);
    [adelay delay]=max(abs(green10sps_1));
    figure(30+Nderiv)
    subplot(2,1,1)
    plot(time(1:NN)-delay*dt,green10sps_1(1:NN)/adelay,'g')
    hold on
    subplot(2,1,2)
    plot(time(1:NN)-delay*dt,green10sps_2(1:NN)/adelay,'g')
    hold on

    green=green10sps;
    green_1=green10sps_1;
    green_2=green10sps_2;

end
