function F=kepler_F(e,M)
%E-偏近点角
%e-偏心率
%M-平近点角
%pi-3.1415926
%误差设计
error=1.e-8;
% if M<pi
%     F=M+e/2;
% else
%     F=M-e/2;
% end
F = M;
ratio=1;
while abs(ratio)>error
    ratio=(e*sinh(F)-F-M)/(e*cosh(F)-1);
    F=F-ratio;
end