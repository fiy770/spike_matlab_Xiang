
%for paper: neuromorphic electronic circuit for building autonomous cognitive systems
%P4-5 equ (6) (8) fig (2)
% no adaption 
clear all;
close all;

I0      = 2.178e-12;         %
I0mem   = 8e-12;             %I0 reset current 8pA
K       = 0.6777;            %0.7054;            % Subthreshold slope factor
Ut      = 25.9e-3;           % Thermal voltage 
Ka		=(K*K)/(K+1);
% model parameters
Vp      = 1.1;               %peak Vmem                                   
Vr      = 0;              %reset Vmem(
Vs      = 0.6;
% Time settings
dt   = 5e-8;            % time step 5us
tmax = 1;               % simulation time
t    = 0:dt:tmax;       % time step vector
% Circuit voltages
Vlk		=0.2607%0.2607%0.2073;
Vthr	=0.2762%0.3027%0.293;
Cp      = 5e-13;  %-----         % Capacitor of adaptation circuit'
Cmem    = 1e-12;  %-----         % Membrane capacitance of the neuron cell

sptimes    = zeros(1,length(t)-1);
%Iahp 	   = zeros(1,length(t));
Vahp 	   = zeros(1,length(t));
DImem      = zeros(1,length(t));
Imem       = zeros(1,length(t)); 
Vmem       = zeros(1,length(t)); 
Ia         = zeros(1,length(t));
NPR        = zeros(1,length(t));
%Frequ 	   = zeros(1,length(t));
%initial 
Vmem(1) 	=Vr;
Imem(1)		=I0mem;
Tspike_ref	=1;
%current define
%percent =linspace(0.2,0.3515,100);%0.4
%percent2=linspace(0.2,0.8,100);%linspace(1,1.7,100);%1.3
%Iahp    =linspace(1e-10,1.1515e-10,31);
%[pp,pp2]=meshgrid(percent, percent2);
%Frequ =zeros(length(percent),length(percent2) );
%Vthr_sim =zeros(length(percent),length(percent2) );
%Vlk_sim =zeros(length(percent),length(percent2) );
% for a=1:length(percent)
%      for b=1:length(percent2)
%Iin2 	=linspace(15e-9,21e-9,100);
% Vthr_sim(a,b)=0;
% Vlk_sim(a,b)=0;
Ith 	=I0*exp((K*Vthr)/Ut) %I0*exp((K*Vthr)/Ut)%Iin*percent*percent2;Iin*percent2(b);%Iin*percent(a)*percent2(b);	%
% Vthr_sim(a,b)=(Ut/K)*log((Ith)/I0);
Itau	=I0*exp((K*Vlk)/Ut)%I0*exp((K*Vlk)/Ut)%Iin*0.1*percent;%Iin*0.1*percent(a);	%7e-9%
% Vlk_sim(a,b)=(Ut/K)*log((Itau)/I0);
Tau		=(Cmem*Ut)/(K*Itau);
% %adaptation para:
% Ica		=1e-8;
% Vlkahp  =0.1;;
% Vthrahp	=0.3;
% Itau_ahp=I0*exp((K/Ut)*Vlkahp);
% Ith_ahp	=I0*exp((K/Ut)*Vthrahp);
% Tau2	=(Cp*Ut)/(K*Itau_ahp);
% Frequ(a,b)=0;
% for in=1:length(Iin2)
    Iin=4.45e-9%Iin2(in)
%for k=1:length(Iahp)
    Iahp=0;
for j=2:length(t)
	if Imem(j-1)<I0mem
            Imem(j) = I0mem;
            Vmem(j) = Vr;
%             Iahp(j)	=0;
	end
	if Vmem(j-1) < Vp
		Ia(j-1)=0.1*Imem(j-1);%I0*exp((Ka*Vmem(j-1))/Ut);	%positive feedback current
% 		Iahp(j)=Iahp(j-1)-(dt/Tau2)*(Iahp(j-1));
        DImem(j)=(((Ith/Itau)*(Iin-Iahp-Itau)+((Ia(j-1)/Itau)-1-(Iahp/Itau))*Imem(j-1)+(Ia(j-1)/Itau)*Ith)*(dt/(Tau*(1+(Ith/Imem(j-1))))));
        Imem(j)=Imem(j-1)+DImem(j);%(((Ith/Itau)*(Iin-Iahp-Itau)+((Ia(j-1)/Itau)-1-(Iahp/Itau))*Imem(j-1)+(Ia(j-1)/Itau)*Ith)*(dt/(Tau*(1+(Ith/Imem(j-1))))));
% 		Imem(j)=Imem(j-1)+(((Ith/Itau)*(Iin-Iahp(j-1)-Itau)+((Ia(j-1)/Itau)-1-(Iahp(j-1)/Itau))*Imem(j-1)+(Ia(j-1)/Itau)*Ith)*(dt/(Tau*(1+(Ith/Imem(j-1))))));
%		Imem(j)=Imem(j-1)+(((Ith/Itau)*(Iin-0-Itau)+((Ia(j-1)/Itau)-1-(0/Itau))*Imem(j-1)+(Ia(j-1)/Itau)*Ith)*(dt/(Tau*(1+(Ith/Imem(j-1))))));
%       Vmem(j)=(Ut/K)*log((Imem(j))/I0);
        Vmem(j)=log(Imem(j)/I0)*(1/Ka)*Ut;
%  		Vahp(j)=(Ut/K)*log((Iahp(j))/I0);
	elseif  Vmem(j-1) >Vp
        Vmem(j-1)   = Vp;
        Vmem(j)		= Vr;    
		Imem(j)		= I0mem;
		Ia(j-1) 	= 0.1*Imem(j);%I0*exp((Ka*Vmem(j-1))/Ut);
%  		Iahp(j)=Iahp(j-1)+(dt/Tau2)*((Ith_ahp/Itau_ahp)*Ica-Iahp(j-1));
%  		Vahp(j)=(Ut/K)*log((Iahp(j))/I0);
		Tspike 		=j-1- Tspike_ref;
		Tspike_ref 	=j-1;
		sptimes(j-1) = Tspike*dt;
%  		Frequ(a,b)		 = 1/sptimes(j-1);
        Frequ		 = 1/sptimes(j-1);
	end
end
%end
% end
spikes   = sptimes(sptimes>0); 
%      end
%  end
 figure(1)
 subplot(8,1,1:3)
 plot(t, Vmem*10^3)
 title('Exponential IF neuron')
 ylabel('Vmem (mV)')	
 xlabel('Time (s)')	
 figure(2)
 plot(t, Imem*10^9);
 ylabel('nA')
 axis([0,30e-4,0,inf])
% figure(3)
%plot(t, Iahp)
%figure(2)
%peaks
%shading interp;
%meshz(pp2,pp,Frequ);
%axis([-inf,inf,-inf,inf,0,inf])
%rotate3d on



