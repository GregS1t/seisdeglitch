function [TT,bbTT]=Fourierleastsquare(time,u,Nhar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% correct for the major dayly waves with least square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sol=24*3600+39*60+35.244;
Nhar1=Nhar+1;
n=length(time);
for in=1:Nhar
	wtide=2.*pi*in;
for i=1:n
	ccT(i,in)=cos(wtide*double(time(i)-time(1)));
	ssT(i,in)=sin(wtide*double(time(i)-time(1)));
end
end
% mean and slope
for i=1:n
	ccT(i,Nhar1)=1.;
	ssT(i,Nhar1)=time(i);
end
% 
% invert thermal amplitude every new sol
clear a
clear b
a(2*Nhar1,2*Nhar1)=0.;
b(2*Nhar1)=0.;
% check if the sol has data
for i=1:n
for in1=1:Nhar1
	for in2=1:Nhar1
		a(2*(in1-1)+1,2*(in2-1)+1)= a(2*(in1-1)+1,2*(in2-1)+1)+ccT(i,in1)*ccT(i,in2);
		a(2*(in1-1)+2,2*(in2-1)+1)= a(2*(in1-1)+2,2*(in2-1)+1)+ssT(i,in1)*ccT(i,in2);
		a(2*(in1-1)+1,2*(in2-1)+2)= a(2*(in1-1)+1,2*(in2-1)+2)+ccT(i,in1)*ssT(i,in2);
		a(2*(in1-1)+2,2*(in2-1)+2)= a(2*(in1-1)+2,2*(in2-1)+2)+ssT(i,in1)*ssT(i,in2);
	end
		b(2*(in1-1)+1)= b(2*(in1-1)+1)+ccT(i,in1)*u(i);
		b(2*(in1-1)+2)= b(2*(in1-1)+2)+ssT(i,in1)*u(i);
end
end
% invert
bbTT=b/a;
TT(1:n)=0.;
for in=1:Nhar1
	TT=TT+bbTT(2*(in-1)+1)*ccT(:,in)'+bbTT(2*(in-1)+2)*ssT(:,in)';
end
end
