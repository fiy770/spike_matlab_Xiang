%% -------- Simultation of exponential integrate an fire neuron --------- %
%                                                                         %
% Author: Quentin AUGEREAU                                                %
% Date: 15/04/2015                                                        %
%                                                                         %
% Description: This matlab file aims to simulate and compute the          %
% differential neuron's behaviour based on exponential integrate and Fire %        
% model                                                                   % 
%                                                                         %
% ----------------------------------------------------------------------- %

%% Initialisation
clear all;
close all;
clc;
%behaviour flag

% Circuit parameters
%I0      = 1.285e-11;        % Dark current 
I0      = 2.178e-12          %3.604e-12;
I0mem   = 8e-12;             %I0 reset current 8pA
K       = 0.6777;            %0.7054;            % Subthreshold slope factor
Ut      = 25.9e-3;           % Thermal voltage 
% model parameters
Vp      = 1.1;               %peak Vmem                                   
Vr      = 0.05;              %reset Vmem
Vs      = 0.6;
% Time settings
dt   = 5e-5;            % time step
tmax = 1;               % simulation time
t    = 0:dt:tmax;       % time step vector
% Circuit voltages
Vth     =0.3;           %1.56e-2;                % Threshold voltage
%Vlk     = 5.48e-3; %-----         % Leakage voltage
Vahp    = 3e-1;                % Adaptation voltage
%Vlka    = 1.56e-2; %-----        % Leakage votage for adaptation circuit
Vlka    = 0;
Vtha    = 1.23e-1;                % Threshold voltage for adaptation circuit
% Circuit capacitors
Cp      = 1e-12;  %-----         % Capacitor of adaptation circuit'
%Cp      = 5e-12;  %-----         % Capacitor of adaptation circuit'
Cmem    = 1e-12;  %-----         % Membrane capacitance of the neuron cell
% Input current
Icons   = 11e-9;                               % constant input current
jstart           = round(0.1*tmax/dt);        % start current subscript
jend             = length(t)-jstart;          % end current subscript
Iin              = zeros(1,length(t));
Iin(jstart:jend) = Icons;                     % initialize the current

%Ith     = I0*exp(K*Vth/Ut);       % Current induced by Vth
Ith =5.74e-9;                       %11.6e-11;

%Itau    = I0*exp(K*Vlk/Ut);       % Current induced by Vlk
%Itau = 13e-10;
Vtau    = 0.32;
Itau    = I0*exp(K*Vtau/Ut);%55.4e-10;
Itaua   = I0*exp(K*Vlka/Ut);      % Current induced by Vlka
Itha    = I0*exp(K*Vtha/Ut);      % Current induced by Vtha
Ica     = I0*exp(K*Vahp/Ut);      % Current induced by Vahp
tau     = (Cmem*Ut/(K*Itau));     % constant tau 
taua    = (Cp*Ut/(K*Itaua));      % constant tau for ahp circuit

%y = flicker(length(t))/(1*10^8); % ne faudrait il pas générer le bruit ?partir de Imem afin d'avoir un NPR stable?
%y = flicker(length(t))/(1*10^10);
y = flicker(length(t))/(4*10^8);
%% Main program

% preallocating
sptimes    = zeros(1,length(t)-1);
sptimes_n  = zeros(1,length(t)-1);
Imem       = zeros(1,length(t)); 
Imem_n     = zeros(1,length(t)); 
Iahp       = zeros(1,length(t));
Iahp_n     = zeros(1,length(t));
Vmem       = zeros(1,length(t)); 
Vmem_n     = zeros(1,length(t)); 
Fimem      = zeros(1,length(t));
Fimem_n    = zeros(1,length(t));
Ia         = zeros(1,length(t));
Ia_n       = zeros(1,length(t));
Iahp_inf   = Itha/Itaua*Ica;
Imem_inf   = zeros(1,length(t));
Imem_inf_n = zeros(1,length(t));
NPR        = zeros(1,length(t));
% initial values
Vmem(1)    = Vr;
Vmem_n(1)  = Vr;
Imem(1)    = I0mem*exp((K*Vr-Vs)/Ut);
Imem_n(1)  = I0mem*exp((K*Vr-Vs)/Ut);

N = 50;
% Main loop
%for k=1:N
%    f_spike = zeros(1,N);
    for j=2:length(t)
        %y = flicker(length(t))/(5*10^10);
        % Time response
        if Imem(j-1)<I0mem
            Imem(j-1) = I0mem;
            %Vmem(j) = (Ut/K)*log((Imem(j))/I0mem)+Vs/K;
            Vmem(j-1) = Vr;
        end
        if Imem_n(j-1)<I0mem
            Imem_n(j-1) = I0mem;
            Vmem_n(j-1) = (Ut/K)*log((Imem_n(j))/I0mem)+Vs/K;
        end
        
        Ia(j-1)       = 1*Imem(j-1);
        Ia_n(j-1)     = 1*(Imem_n(j-1)+y(j-1));
        %Ia_n(j-1)       = 1*Imem_n(j-1);
        Fimem(j-1)    = (Ia(j-1)/Itau)*(Imem(j-1)+Ith);
        Fimem_n(j-1)  = (Ia_n(j-1)/Itau)*(Imem_n(j-1)+Ith);
        Imem_inf(j)   = Ith/Itau*(Iin(j-1)-Iahp(j-1)-Itau);
        Imem_inf_n(j) = Ith/Itau*(Iin(j-1)-Iahp_n(j-1)-Itau);
        %NPR(j-1)      = y(j-1)^2/Imem(j-1)^2;
        NPR(j-1)      = y(j-1)^2/Ia_n(j-1)^2;
        if Vmem(j-1)<Vp 
            % memomory current
            Imem(j) = Imem(j-1)+dt/(tau*(1+Ith/Imem(j-1)))*(Ith/Itau*...
            (Iin(j-1)-Iahp(j-1)-Itau)-Imem(j-1)*(1+Iahp(j-1)/Itau)+ ...
            Fimem(j-1));
            % Adaptation current
            Iahp(j)    = Iahp(j-1) + dt/taua*(-Iahp(j-1));
            Vmem(j)    = (Ut/K)*log((Imem(j))/I0)+Vs/K;
        elseif abs(Vmem(j-1))>Vp % spike behaviour
            Vmem(j-1)    = Vp;    
            Imem(j-1)    = I0mem*exp((K*Vp-Vs)/Ut);
            Ia(j-1)      = Imem(j-1);
            Vmem(j)      = Vr;                     % reset voltage
            Imem(j)      = I0mem*exp((K*Vr-Vs)/Ut); 
            Iahp(j)      = Iahp(j-1) + dt/taua*(-Iahp(j-1)+Iahp_inf);
            indp         = j-1 + (Vp-Vmem(j-1))/(Vmem(j)-Vmem(j-1));
            sptimes(j-1) = indp*dt;       
        end    
        if Vmem_n(j-1)<Vp
            % memory current with noise
            Imem_n(j) = Imem_n(j-1)+dt/(tau*(1+Ith/Imem_n(j-1)))*(Ith/Itau*...
                (Iin(j-1)-Iahp(j-1)-Itau)-Imem_n(j-1)*(1+Iahp(j-1)/Itau)+ ...
                Fimem_n(j-1));
            %Imem_n(j)  = -y(j)+Imem_n(j-1)+y(j-1)+dt/(tau*(1+Ith/       ...
            %    (Imem_n(j-1)+y(j-1))))*(Ith/Itau*(Iin(j-1)-Iahp(j-1)-   ...
            %    Itau)-(Imem_n(j-1)+y(j-1))*(1+Iahp(j-1)/Itau)+Fimem_n(j-1));
            %Vmem_n(j)  = (Ut/K)*log((Imem_n(j))/I0);
            Vmem_n(j)  = (Ut/K)*log((Imem_n(j))/I0)+Vs/K;
            % adaptation current
            Iahp_n(j)    = Iahp_n(j-1) + dt/taua*(-Iahp_n(j-1)); 
        elseif abs(Vmem_n(j-1))>Vp
            Vmem_n(j-1) = Vp;
            Imem_n(j-1) = I0mem*exp((K*Vp-Vs)/Ut);
            Vmem_n(j)   = Vr;
            Imem_n(j)   = I0mem*exp((K*Vr-Vs)/Ut); 
            
%             Imem_n(j-1)    = I0*exp(K*Vp/Ut)-y(j-1);
%             Vmem_n(j-1)    = Ut/K*log(Imem_n(j-1)/I0);
%             Imem_n(j)      = I0*exp(K*Vr/Ut)-y(j);
%             Vmem_n(j)      = Ut/K*log(Imem_n(j)/I0);
            
            Iahp_n(j)      = Iahp_n(j-1) + dt/taua*(-Iahp_n(j-1)+Iahp_inf);
            indp_n         = j-1 + (Vp-Vmem_n(j-1))/(Vmem_n(j)-Vmem_n(j-1));
            sptimes_n(j-1) = indp_n*dt;
        end
    end
    
   spikes   = sptimes(sptimes>0);
   spikes_n = sptimes_n(sptimes_n>0);
%    f_spike(k) = spikes_n(5);
%end

    if (length(spikes)==length(spikes_n))
        result=spikes-spikes_n;
    else
        result = sptimes(5000:end)-sptimes_n(5000:end);
        result = result(result>0|result<0);
        %result=spikes(10:end)-spikes_n(10:end);
    end
cpt   = 0;
cpt_n = 0;
for k = jstart:(2*jstart);
    if (Vmem_n(k))== Vp
        cpt_n = cpt_n + 1;
    elseif (Vmem(k) == Vp)
        cpt = cpt + 1;
    end
end


scale = 1:1:length(result);
% 
% 
% spikes = sptimes(sptimes>0);
pif = 10*log(NPR(jstart:jend));
NPR_mean  = sum(pif)/length(NPR(jstart:jend));


freq = zeros(1, length(spikes)-1);

%frequency
for k=1:length(spikes)-1
    freq(k) = 1/(spikes(k+1)-spikes(k));
end


%% plots
figure(1)
subplot(8,1,1:3)
plot(t, Vmem*10^3)
title('Exponential IF neuron')
ylabel('Vmem (mV)')
%axis([0.35 0.4 -10 1020])
%axis([0.2 0.4 -10 1020])
%axis([0 tmax 50 1120])
grid on;
%grid on
subplot(8,1,4:5)
plot(t, Imem*1e9, 'red')
ylabel('Imem (nA)'); xlabel('t (s)')
%axis([0.2 0.4 min(Imem)*1e9-0.02 max(Imem)*1e9+0.02])
%axis([0.2 0.4 -1 1.5*10^9])
axis([0 tmax -1 250])
subplot(8,1,6:7)
plot(t,Iahp*1e9, 'black')
ylabel('Iahp (nA)')
%axis([0.2 0.4 min(Iahp)*1e9-0.02 max(Iahp)*1e9+0.02])
axis([0 tmax min(Iahp)*1e9-0.02 max(Iahp)*1e9+0.02])
%grid on
subplot(8,1,8)
plot(t,Iin*1e9, 'magenta')
xlabel('t (s)'); ylabel('Iin (nA)')
%axis([0.2 0.4 min(Iin)*1e9-10 max(Iin)*1e9+10])
axis([0 tmax min(Iin)*1e9-10 max(Iin)*1e9+10])
%grid on;


figure(2)
subplot(8,1,1:3)
plot(t, Vmem_n*10^3)
title('Exponential IF neuron')
ylabel('Vmem_n (mV)')
%axis([0.35 0.4 -10 1020])
axis([0 tmax 800 1120])
grid on;
%grid on
subplot(8,1,4:5)
plot(t, Imem_n*1e9, 'red')
ylabel('Imem_n (nA)'); xlabel('t (s)')
%axis([0 tmax -1e+9-0.02 1e9+0.02])
subplot(8,1,6:7)
plot(t,Iahp_n*1e9, 'black')
ylabel('Iahp (nA)')
%axis([0 tmax min(Iahp)*1e9-0.02 max(Iahp)*1e9+0.02])
%grid on
subplot(8,1,8)
plot(t,Iin*1e9, 'magenta')
xlabel('t (s)'); ylabel('Iin (nA)')
axis([0 tmax min(Iin)*1e9-10 max(Iin)*1e9+10])
%grid on;

figure(3)
plot(t,y)
title('Noise current (A)')
xlabel('time (s)'); ylabel('Noise current (A)');


x=-3*10^-4:2*10^-5:3*10^-4;
figure (4)
subplot(1,2,1)
plot(scale, result)
subplot(1,2,2)
hist(result,x)
title('histogram of space between noisy spikes and normal spikes');

NPR = NPR(1:end-1);
figure(5)
plot(t(1:end-1), 20*log10(NPR))
axis([0 tmax -100 0])
hold on;
plot(t(1:end-1), (NPR_mean*ones(1,length(t)-1)), 'black');
title('NPR dB');

figure
subplot(2,1,1)
plot(t, Vmem*10^3)
title('Exponential IF neuron')
ylabel('Vmem (mV)')
xlabel('time (s)')
axis([0.15 0.4 800 1120])
%axis([0 tmax -10 1020])
grid on;
subplot(2,1,2)
plot(t, Vmem_n*10^3)
title('Exponential IF neuron')
ylabel('Vmem_n (mV)')
xlabel('time (s)')
axis([0.15 0.4 800 1120])
%axis([0 tmax -10 1020])
grid on;

% [f,x1]=hist(10*log(NPR(jstart:jend)),200);%# create histogram from a normal distribution.
% x1 = x1+100;
% g=1/sqrt(2*pi)*exp(-0.5*x1.^2);%# pdf of the normal distribution
% 
% %#METHOD 1: DIVIDE BY SUM
% figure
% bar(x1(78:end),f(78:end)/sum(f(78:end)));hold on
% plot(x1(78:end),g(78:end),'r');hold off
% 
% %#METHOD 2: DIVIDE BY AREA
% figure
% bar(x1(78:end),f(78:end)/trapz(x1(78:end),f(78:end)));hold on
% plot(x1(78:end),g(78:end),'r');hold off
% 
% 
% xx = -1000:10:0;
% NPR_h = hist(10*log(NPR(jstart:jend)),xx);
% NPR_h = NPR_h./sum(NPR_h);
% 
% figure
% subplot(2,1,1)
% bar(NPR_h, 'yellow');
% axis([-200 0 0 1])
% %hold on;
% subplot(2,1,2)
% histfit(10*log(NPR(jstart:jend)));
% axis([-200 0 0 4000])
% % figure
% % hist(f_spike)

figure
plot(freq)
