function E=kepler_E(e,M)
%E-偏近点角
%e-偏心率
%M-平近点角
%pi-3.1415926
%误差设计
error=1.e-8;
if M<pi
    E=M+e/2;
else
    E=M-e/2;
end
ratio=1;
while abs(ratio)>error
    ratio=(E-e*sin(E)-M)/(1-e*cos(E));
    E=E-ratio;
end