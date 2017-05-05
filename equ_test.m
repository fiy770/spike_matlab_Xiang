 clear all
close all
a= linspace(0,20e-9,1000);  %Imem Imem_max=1.9774e-7
b=0.5e-9;   %Ith
c=6.5e-9;    %Iin
d=2e-9;     %Itau
e=0.5;     %slope of equation
f= zeros(1,1000);   % i_fb
gain = 1; % abs of diff pair gain
vdd = 1.8;
vthp = 0.5; %Pmos Vth
I0= 2.178e-12; 
K = 0.6777;            %0.7054;            % Subthreshold slope factor
Ut= 25.9e-3;
Ka		=(K*K)/(K+1);
const=(0.1)*ones(1,1000)
Vmem = zeros(1,1000)
vout = zeros(1,1001)

for i=1:length(a)   
    if vout(i) < vthp 
        f(i) = I0*exp((K/Ut)*vout(i))
    end
    equ(i)=5e-7*((b*c-d*b-d*a(i)+f(i)*(a(i)+b))*(K/(5e-13*0.0259*(1+(b/a(i)))))); %b*c+(e*(a(i)^2))-a(i)*d;
    normalize_equ(i)=5e-7*((b*c-d*b-d*a(i)+e*a(i)*a(i))*(K/(5e-13*0.0259*(1+(b/a(i))))))*(1/a(i));
    diffequ(i)=e-((b*c)/((a(i)+b)^2));
    Vmem(i)=log(a(i)/I0)*(1/Ka)*Ut;
    if Vmem(i)> 0.4 
       vout(i+1)=(Vmem(i)-0.4)*gain;
       if vout(i+1) > vdd
           vout(i+1) = vdd;
       end
    end
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
 title('dImem*dt VS Imem')
 xlabel('Imem')	
 ylabel('dImem')	
figure(2);
plot(a,normalize_equ);
title('dImem*dt/Imem VS Imem')
ylabel('dImem*dt/Imem')	
xlabel('Imem')	
figure(3)
plot(a,diffequ,a,const)
title('d(dImem*dt/Imem*delta)/dImem VS Imem')
ylabel('d(dImem*dt/Imem*delta)/dImem')	
xlabel('Imem')	
figure(4)
plot(a,Vmem)
figure(5)
plot(a,vout(1,1000))