function [wdgD,wdgDD,gwD0,gwD1,var0,var,reduc]=betershift3(green2sps_7,green2sps_8,wdg,idx0, dt, ND, DeltaN, g4delay)
%%%%%%%%% better shift %%%%%%%%%%%%%%%%%%%%%%
gwD=green2sps_7(g4delay-ND:g4delay+ND);
gw1D=green2sps_8(g4delay-ND:g4delay+ND);
aA(1:4,1:4)=0.;
bb(1:4)=0.;
ig=ND+1;
idx1=idx0+DeltaN;
for i=-ND:ND
    aA(1,1)=aA(1,1)+gwD(i+ig)*gwD(i+ig);
    aA(1,2)=aA(1,2)+gwD(i+ig)*gw1D(i+ig);
    aA(1,3)=aA(1,3)+gwD(i+ig);
    aA(1,4)=aA(1,4)+gwD(i+ig)*(i+ig);
    aA(2,2)=aA(2,2)+gw1D(i+ig)*gw1D(i+ig);
    aA(2,3)=aA(2,3)+gw1D(i+ig);
    aA(2,4)=aA(2,4)+gw1D(i+ig)*(i+ig);
    aA(3,3)=aA(3,3)+1;
    aA(3,4)=aA(3,4)+(i+ig);
    aA(4,4)=aA(4,4)+(i+ig)^2;
    bb(1)=bb(1)+gwD(i+ig)*wdg(idx1+i);
    bb(2)=bb(2)+gw1D(i+ig)*wdg(idx1+i);
    bb(3)=bb(3)+wdg(idx1+i);
    bb(4)=bb(4)+(i+ig)*wdg(idx1+i);
end
aA(2:4,1)=aA(1,2:4);
aA(3:4,2)=aA(2,3:4);
aA(4,3)=aA(3,4);
cc=inv(aA)*bb';

for i=-ND:ND
gwD0(i+ig)=cc(3)+cc(4)*(i+ig);
gwD1(i+ig)=cc(1)*gwD(i+ig)+cc(2)*gw1D(i+ig)+cc(3)+cc(4)*(i+ig);
wdgDD(idx1+i)=wdg(idx1+i)-cc(1)*gwD(i+ig)-cc(2)*gw1D(i+ig);
wdgD(idx1+i)=cc(1)*gwD(i+ig)+cc(2)*gw1D(i+ig);
end
%%%%%%%%%%%%%
%for i=-2:2
%wdgDD(idx1+i)=wdgDD(idx1-3)+(wdgDD(idx1+3)-wdgDD(idx1-3))*(i+3)/6;
%end
%%%%%%%%%%%%%%
var=0.;
var0=0.;
for i=-ND:ND
var0=var0+(wdg(idx1+i)-gwD0(i+ig))^2;
var =var+ (wdg(idx1+i)-gwD1(i+ig))^2;
end
reduc=var/var0;
%%%%%%%%%%%%%%
