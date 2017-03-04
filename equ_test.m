a= linspace(0,10e-9,1000);  %Imem
b=1.5e-9;   %Ith
c=20e-9;    %Iin
d=7e-9;     %Itau
equ=b*c+1*a.^2-a.*d;
figure(3)
plot(a,equ)