function F=kepler_F(e,M)
%E-ƫ�����
%e-ƫ����
%M-ƽ�����
%pi-3.1415926
%������
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