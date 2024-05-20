function value = mypoly(x, y,small)
small_x = floor((x+1)/2);
small_y = floor((y+1)/2);
if small_x==(x+1)/2 && small_y==(y+1)/2
    value = small(small_x,small_y);
    return;
end
f1=int32(0);
f2=int32(0);
f3=int32(0);
f4=int32(0);
f5=int32(0);
f6=int32(0);
for localx=-2:1:2
    for localy=-2:1:2
        f1=f1+int32(localx^2*int32(small(uint32(small_x+localx),uint32(small_y+localy))));
        f2=f2+int32(localy^2*int32(small(uint32(small_x+localx),uint32(small_y+localy))));
        f3=f3+int32(localx*localy*int32(small(uint32(small_x+localx),uint32(small_y+localy))));
        f4=f4+int32(localx*int32(small(uint32(small_x+localx),uint32(small_y+localy))));
        f5=f5+int32(localy*int32(small(uint32(small_x+localx),uint32(small_y+localy))));
        f6=f6+int32(small(uint32(small_x+localx),uint32(small_y+localy)));
    end
end
A = [170 100 0 0 0 50;
     100 170 0 0 0 50;
     0 0 100 0 0 0;
     0 0 0 50 0 0;
     0 0 0 0 50 0;
     50 50 0 0 0 25];
D = [double(f1);double(f2);double(f3);double(f4);double(f5);double(f6)];
C = A\D;
a = C(1,1);
b = C(2,1);
c = C(3,1);
d = C(4,1);
e = C(5,1);
h = C(6,1);
if small_x<(x+1)/2 && small_y==(y+1)/2
    value = a/4-d/2+h;
    return;
end
if small_x==(x+1)/2 && small_y<(y+1)/2
    value = b/4-e/2+h;
    return;
end
if small_x<(x+1)/2 && small_y<(y+1)/2
    value = a/4+b/4+c/4-d/2-e/2+h;
    return;
end
end