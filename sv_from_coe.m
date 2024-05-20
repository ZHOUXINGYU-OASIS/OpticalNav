function sv=sv_from_coe(M,a,e,i,w,om,mu)

i=i*pi/180;
om=om*pi/180;
w=w*pi/180;
if e<1
    E=kepler_E(e,M);%牛顿法解开普勒方程（已知e,平近点角M,求偏近点角E）
    sita=2*atan(((1+e)/(1-e))^0.5*tan(E/2));
    if E>=pi
        sita=sita+2*pi;
    end
else
    F=kepler_F(e,M);%牛顿法解双曲线轨道的开普勒方程
    sita=2*atan(((e+1)/(e-1))^0.5*tanh(F/2));
    if F>=pi
        sita=sita+2*pi;
    end
end

p=a*(1-e*e);
r=p/(1+e*cos(sita));
Rp=[r*cos(sita);r*sin(sita);0];
vr=(mu/p)^(1/2)*e*sin(sita);
vu=(mu/p)^(1/2)*(1+e*cos(sita));
Vp=[vr*cos(sita)-vu*sin(sita);vr*sin(sita)+vu*cos(sita);0];
Lpi=[cos(w)*cos(om)-sin(w)*cos(i)*sin(om) cos(w)*sin(om)+sin(w)*cos(i)*cos(om) sin(w)*sin(i);
-sin(w)*cos(om)-cos(w)*cos(i)*sin(om) -sin(w)*sin(om)+cos(w)*cos(i)*cos(om) cos(w)*sin(i);
sin(i)*sin(om) -sin(i)*cos(om) cos(i);];
R=Lpi'*Rp;
V=Lpi'*Vp;
sv=[R; V];