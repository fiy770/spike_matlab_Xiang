a= linspace(0,200e-9,1000);  %Imem Imem_max=1.9774e-7
b=6e-9;   %Ith
c=20e-9;    %Iin
d=2e-9;     %Itau
e=0.01;     %slope of equation
I0= 2.178e-12; 
K = 0.6777;            %0.7054;            % Subthreshold slope factor
Ut= 25.9e-3;
equ=b*c+(e*(a.^2))-a.*d;
lowest_value= b*c-e*(d/(2*e))^2;
delta=(0.6777/(5e-13*0.0259))*5e-7;
mim_delta_Imem=delta*lowest_value;
center = d/(2*e);
Vthr=(Ut/K)*log((b)/I0);
Vtau=(Ut/K)*log((d)/I0);
figure(3);
plot(a,equ);