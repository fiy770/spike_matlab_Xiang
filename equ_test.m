clear all
close all
a= linspace(0,20e-9,1000);  %Imem Imem_max=1.9774e-7
b=5e-9;   %Ith
c=5e-9;    %Iin
d=2e-9;     %Itau
e=0.1;     %slope of equation
I0= 2.178e-12; 
K = 0.6777;            %0.7054;            % Subthreshold slope factor
Ut= 25.9e-3;
Ka		=(K*K)/(K+1);
for i=1:length(a)   
    equ(i)=5e-7*((b*c-d*b-d*a(i)+e*a(i)*a(i))*(K/(5e-13*0.0259*(1+(b/a(i))))));%b*c+(e*(a(i)^2))-a(i)*d;
    normalize_equ(i)=5e-7*((b*c-d*b-d*a(i)+e*a(i)*a(i))*(K/(5e-13*0.0259*(1+(b/a(i))))))*(1/a(i));
end
lowest_value= b*c-e*(d/(2*e))^2;
delta=(0.6777/(5e-13*0.0259))*5e-7;
mim_delta_Imem=delta*lowest_value;
center = (d-e*b)/(2*e);
Vthr=(Ut/K)*log((b)/I0);
Vtau=(Ut/K)*log((d)/I0);
Vmem_c=log(center/I0)*(1/Ka)*Ut;
figure(1);
plot(a,equ);
 title('dImem VS Imem')
 xlabel('Imem')	
 ylabel('dImem')	
figure(2);
plot(a,normalize_equ);
title('dImem VS Imem')
ylabel('dImem/Imem')	
xlabel('Imem')	