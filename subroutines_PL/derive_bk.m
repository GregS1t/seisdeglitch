function [ ts_new ] = derive_bk(time,ts_or,dt )
%This script derives by tqking into account gaps in the data
%%%Find gaps
ts_new=[];
N1=[];
N2=[];
N1=1;
tol=.5;
for ii=2:length(time)
    if 86400*(time(ii)-time(ii-1))>((1+tol)*dt)
        N2=[N2;ii-1];
        N1=[N1;ii];
    end
end
    N2=[N2; length(time)];

NGAPS=length(N2); %Number of gaps+1


%%%Loop on data without gaps
ts_c2=[];
clear w
for NT=1:NGAPS
    ts=ts_or(N1(NT):N2(NT));
    NS=length(ts);
    w=[0:NS-1]'/NS/dt*1i*2*pi;
    ft=fft(ts);
    ft=ft.*w;
    ts=ifft(ft,'symmetric');
    ts_new=[ts_new; ts];
   
    

end
end

