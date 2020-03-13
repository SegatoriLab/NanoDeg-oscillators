%% Compare fast degradation with slow degradation of activator

clear all; 
close all; 
parameter_AR;

%time step
t0 = 0;
tf = 600;
dt = .01;
t = t0:dt:tf;

z = 0; kn = 0; 
inicondd = [1e-3 1e-3 z 0];
deltaA = log(2)/4;  
deltaB = log(2)/4;  
[td, xd] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);

deltaA = log(2)/.876; %deltaA = log(2)/.3 to simulate half-life of 18 mins
deltaB = log(2)/4;
[td2, xd2] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);

z = 5; kn = deltaN*z; 
inicondd = [1e-3 1e-3 z 0];
deltaA = log(2)/4;  
deltaB = log(2)/4;  
deltaN = log(2)/.9;
[td3, xd3] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);

figure(1); subplot(311);  plot(td,xd(:,1)); hold all; plot(td2,xd2(:,1)); plot(td3,xd3(:,1));
title('Activator'); 
legend('\tau_{1/2} = 4hr', '\tau_{1/2} = 18 min','\tau_{1/2} = 4hr, \tau_{NB} = 54 min');
set (gca, 'YScale', 'log');
xlabel('time'), ylabel('nM');

figure(1); subplot(312);  plot(td,xd(:,2)); hold all; plot(td2,xd2(:,2)); plot(td3,xd3(:,2))
title('Repressor'); 
%legend('\tau_{1/2} = 4hr', '\tau_{1/2} = 18 min','\tau_{1/2} = 4hr, \tau_{NB} = 54 min');

set (gca, 'YScale', 'log');
xlabel('time'), ylabel('nM');

figure(1); subplot(313);  plot(xd(:,1),xd(:,2)); hold all; plot(xd2(:,1),xd2(:,2)); plot(xd3(:,1),xd3(:,2)) 
title('Activator-Repressor Phase Diagram'); 
%legend('\tau_{1/2} = 4hr', '\tau_{1/2} = 18 min','\tau_{1/2} = 4hr, \tau_{NB} = 54 min');
xlabel('Activator (nM)'); ylabel('Repressor (nM)');
set (gca, 'YScale', 'log'); set(gca,'XScale', 'log');

%% Generate Amplitude Chart x Nanobody concentration for Activator attack (switchable)

clear all; 
close all; 
parameter_AR;

%time step
t0 = 0;
tf = 150;
dt = .01;
t = t0:dt:tf;

zvec = [0:.1:10];

for i =1:length(zvec)
    z = zvec(i); 
    kn = deltaN*z;
    inicondd = [1 10 z 0];
    [td, xd] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);
    
    activAmp(i) = range(xd(td>50,1));
    reprAmp(i) = range(xd(td>50,2));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-mean(xd(td>50,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
        
end

figure(2); subplot(221);  plot(zvec, activAmp); 
ylabel('nM');title ('Activator Amplitude');
subplot(222); plot (zvec, period); title ('Oscillation Period'); 
ylabel('hr'); xlabel('Total Nanobody (nM)');

%% Measuring half-life effect of nanobody

clear all; 
close all; 
parameter_AR;

hl = 4*10.^[-2:.1:0];
dnvec = log(2)./hl;


%time step
t0 = 0;
tf = 150;
dt = .01;
t = t0:dt:tf;

for i =1:length(dnvec)
    z = 5; 
    deltaN = dnvec(i);
    kn = deltaN*z;
    inicondd = [1 10 z 0];
    [td, xd] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);
    
    activAmp2(i) = range(xd(td>50,1));
    reprAmp2(i) = range(xd(td>50,2));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-mean(xd(td>50,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period2(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
        
end

figure(2); subplot(223);  plot(hl, activAmp2); 
ylabel('nM');title ('Activator Amplitude');
subplot(224); plot (hl, period2); title ('Oscillation Period'); 
ylabel('hr'); xlabel('Nanobody half-life (hrs)');






%% Comparing weak binding with strong binding keeping kd
close all; clear all;

parameter_AR;

%time step
t0 = 0;
tf = 200;
dt = .01;
t = t0:dt:tf;

zvec = 0:.1:10;

for i =1:length(zvec)
    z = zvec(i); 
    kn = deltaN*z;
    kon = 2.7648; koff = .6264;
    inicondd = [1 10 z 0];
    [td, xd] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);
    
    activAmp(i) = range(xd(td>50,1));
    reprAmp(i) = range(xd(td>50,2));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-mean(xd(td>50,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
        
end

for i =1:length(zvec)
    z = zvec(i); 
    G = 3;
    kon = G*2.7648; koff = G*.6264;
    kn = deltaN*z;
    inicondd = [1 10 z 0];
    [td, xd] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);
    
    activAmp2(i) = range(xd(td>50,1));
    reprAmp2(i) = range(xd(td>50,2));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-mean(xd(td>50,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period2(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
        
end

for i =1:length(zvec)
    z = zvec(i); 
    G = 10;
    kon = G*2.7648; koff = G*.6264;
    kn = deltaN*z;
    inicondd = [1 10 z 0];
    [td, xd] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);
    
    activAmp3(i) = range(xd(td>50,1));
    reprAmp3(i) = range(xd(td>50,2));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-mean(xd(td>50,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period3(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
        
end

figure(2); subplot(211);  plot(zvec, activAmp); hold on; plot(zvec, activAmp2); plot(zvec, activAmp3);
ylabel('nM');title ('Effect of "retroactivity"'); legend('G=1','G=3','G=10');


for i =1:length(zvec)
    z = zvec(i); 
    G = 3;
    kon = G*2.7648; koff = .6264;
    kn = deltaN*z;
    inicondd = [1 10 z 0];
    [td, xd] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);
    
    activAmp4(i) = range(xd(td>50,1));
    reprAmp4(i) = range(xd(td>50,2));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-mean(xd(td>50,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period4(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
        
end

for i =1:length(zvec)
    z = zvec(i); 
    G = 10;
    kon = G*2.7648; koff = .6264;
    kn = deltaN*z;
    inicondd = [1 10 z 0];
    [td, xd] = ode23s(@(t, x) ActivatorRepressorDegNB_act(t,x), t, inicondd);
    
    activAmp5(i) = range(xd(td>50,1));
    reprAmp5(i) = range(xd(td>50,2));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-mean(xd(td>50,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period5(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
        
end

figure(2); subplot(212);  plot(zvec, activAmp); hold on; plot(zvec, activAmp4); plot(zvec, activAmp5); 
ylabel('nM');title ('Effect of affinity'); legend('kon = 2.8', 'kon = 8.3', 'kon = 27.6');
    
   



