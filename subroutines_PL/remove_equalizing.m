function [u_noequal]=remove_equalizing(u,Uz,Up,uA0,ugain,Rz,Rp,rA0,rgain,dtu)

Nu=length(u);
%%%%% put all sensors with the a slepian TF
NF=floor(Nu/2)+1;
dfu=1./(Nu*dtu);
Frequ=[0:NF-1]*dfu;
omega=2.*pi*1i*Frequ;
%%%% remove the two first zero and keep the other ones
tfu(1:NF)=complex(1.,0.);
%%%%
for ii=1:length(Uz)
        if Uz(ii) ~= 0.
                tfu(1:NF)=tfu(1:NF).*(omega(1:NF)-Uz(ii));
        end
end
for ii=1:length(Up)
    tfu(1:NF)=tfu(1:NF)./(omega(1:NF)-Up(ii));
end
tfr(1:NF)=complex(1.,0.);
%%%%
for ii=1:length(Rz)
        if Rz(ii) ~= 0.
                tfr(1:NF)=tfr(1:NF).*(omega(1:NF)-Rz(ii));
        end
end
for ii=1:length(Rp)
    tfr(1:NF)=tfr(1:NF)./(omega(1:NF)-Rp(ii));
end
su=fft(u);
for i=1:NF
        su(i)=su(i)*uA0*ugain/(rA0*rgain)*tfu(i)/tfr(i);
end
u_noequal=ifft(su,'symmetric');
end
