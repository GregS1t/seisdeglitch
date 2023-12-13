Nu=length(u);
Nv=length(v);
Nw=length(w);

%%%%% put all sensors with the TF of the u component
NF=floor(Nu/2)+1;
dtu=(ut(2)-ut(1))*24*3600;
dfu=1./(Nu*dtu);
Frequ=[0:NF-1]*dfu;
omega(1:NF)=[0:NF-1]*1i*2*pi*dfu;
%%%% remove the two first zero and keep the other ones for Nderiv=0;
tfu(1:NF)=complex(1.,0.);
tfv(1:NF)=complex(1.,0.);
tfw(1:NF)=complex(1.,0.);
%%%% 
for ii=1:length(Uz)
	if Uz(ii) ~= 0.
    		tfu(1:NF)=tfu(1:NF).*(omega(1:NF)-Uz(ii));
	end
end
for ii=1:length(Up)
    tfu(1:NF)=tfu(1:NF)./(omega(1:NF)-Up(ii));
end
%%%% 
for ii=1:length(Vz)
	if Vz(ii) ~= 0.
    		tfv(1:NF)=tfv(1:NF).*(omega(1:NF)-Vz(ii));
	end
end
for ii=1:length(Vp)
    tfv(1:NF)=tfv(1:NF)./(omega(1:NF)-Vp(ii));
end
%%%% 
for ii=1:length(Wz)
	if Wz(ii) ~= 0.
    		tfw(1:NF)=tfw(1:NF).*(omega(1:NF)-Wz(ii));
	end
end
for ii=1:length(Wp)
    tfw(1:NF)=tfw(1:NF)./(omega(1:NF)-Wp(ii));
end
%%%%%%
if plot_data==1
    figure(60)
    subplot(4,1,1)
    plot(ut,u,'k', "DisplayName", "u")
    subplot(4,1,2)
    plot(vt,v,'k', "DisplayName", "v")
    hold on
    subplot(4,1,3)
    plot(wt,w,'k', "DisplayName", "w")
    hold on
end
if Nv == Nu
	sv=fft(v);
	for i=1:NF
		sv(i)=sv(i)*uA0*ugain/vA0/vgain*tfu(i)/tfv(i);
	end
	v=ifft(sv,'symmetric');
	Vz=Uz;
	Vp=Up;
end
if Nw == Nw
	sw=fft(w);
	for i=1:NF
		sw(i)=sw(i)*uA0*ugain/wA0/wgain*tfu(i)/tfw(i);
	end
	w=ifft(sw,'symmetric');
	Wz=Uz;
	Wp=Up;
end


if plot_data ==1
    figure(60)
    subplot(4,1,1)
    plot(ut,u,'k', "DisplayName", "u_{eq}")
    ylabel("u")
    legend
    grid on
    subplot(4,1,2)
    plot(vt,v,'r',  "DisplayName", "v_{eq}"')
    ylabel("v")
    legend
    grid on
    subplot(4,1,3)
    plot(wt,w,'r',  "DisplayName", "w_{eq}")
    ylabel("w")
    legend

    grid on
end
