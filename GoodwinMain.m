%% Standalone Goodwin with oscillation
close all;
clear;
parameter_Goodwin;

%time step
t0 = 0;
tf = 500;
dt = .01;
tt = t0:dt:tf;

delvect = [.5 10 15];

for i=1:length(delvect)

    delay = delvect(i);

    sol = dde23(@GoodwinOsc, delay, @history1, tt);
    
    td = tt;
    yd = deval(sol,tt);
    ydd(:,tt>delay) = deval(sol, tt(tt>delay)-delay);
    ydd(:,tt<= delay) = history1(tt);
    
    figure(1); %hold off;
    subplot(2,1,1);
    plot(td, yd(1,:)); hold all;
    title('Repressor');
    xlabel('Time');
    ylabel('nM');
    set (gca, 'YScale','log');
    
    
    figure(1);
    subplot(2,1,2);
    plot(ydd(1,:), yd(1,:)); hold all;
    xlabel('Repressor (t-\tau)');
    ylabel('Repressor (t)');
    set (gca, 'XScale','log');
    set (gca, 'YScale','log');
    


end

legend('delay = .5h','delay = 10h', 'delay = 15h');

%% Standalone Goodwin with oscillation, change in degradation rate
close all;
clear;
parameter_Goodwin;

%time step
t0 = 0;
tf = 50;
dt = .01;
tt = t0:dt:tf;

degvec = log(2)./[11 2/3];

for i=1:length(degvec)

    deltaA = degvec(i);
    delay = 0.5;

    sol = dde23(@GoodwinOsc, delay, @history1, tt);
    
    td = tt;
    yd = deval(sol,tt);
    ydd(:,tt>delay) = deval(sol, tt(tt>delay)-delay);
    
    figure(1); %hold off;
    subplot(2,1,1);
    plot(td, yd(1,:)); hold all;
    title('Repressor');
    xlabel('Time');
    ylabel('nM');
    set (gca, 'YScale','log');
    
    
    figure(1);
    subplot(2,1,2);
    plot(ydd(1,:), yd(1,:)); hold all;
    xlabel('Repressor (t-\tau)');
    ylabel('Repressor (t)');
    set (gca, 'XScale','log');
    set (gca, 'YScale','log');
    


end

legend('half-life = 11h','delay = 38min');

%% Test for Goodwin. Trigger oscillation by changing rates

close all;
clear;
parameter_Goodwin;


%time step
t0 = 0;
tf = 50;
dt = .01;
tt = t0:dt:tf;

%initial conditions
G = [1/100 1/30 1];
kn = 30;


for i =1:length(G)
    koff = .6264*G(i);
    kon = 2.7648*G(i);
    
    delay = .5;

    sol = dde23(@GoodwinOscNB, delay, @history3, tt);
    
    td = tt;
    yd = deval(sol,tt);
    ydd(:,tt>delay) = deval(sol, tt(tt>delay)-delay);
    
    %calculating equilibrium
    AA = 10.^(-2:.001:4);
    AA2 = (1+koff/deltaN)*kn./(deltaN+kon*AA+koff);
    AA3 = kn/deltaN - AA2;
    F_Aa = pA_T*((beta1*AA.^m + beta2*ka^m)./(AA.^m + ka^m));

    rr = (F_Aa + koff*AA3)./(deltaA+kon*AA2)-AA;
    zerocrossing = rr(1:end-1).*rr(2:end)<0;
    equil = AA(zerocrossing);
    
    disp(equil);
    
    figure(1); %hold off;
    subplot(2,2,1);
    plot(td, yd(1,:)); hold all;
    title('Repressor');
    xlabel('Time');
    ylabel('nM');
    set (gca, 'YScale','log');
    
    figure(1);
    subplot(2,2,2);
    plot(td, yd(2,:)); hold all;
    title('nanobody');
    xlabel('Time');
    ylabel('nM');
    set (gca, 'YScale','log');
    
    figure(1);
    subplot(2,2,3);
    plot(td, yd(2,:)+ yd(3,:)); hold all;
    title('Total Nanobody');
    xlabel('Time');
    ylabel('nM');
    set (gca, 'YScale','log');
    
    figure(1);
    subplot(2,2,4);
    plot(ydd(1,:), yd(1,:)); hold all;
    title('Repressor-Nanobody Complex');
    xlabel('Repressor (t-\tau)');
    ylabel('Repressor (t)');
    set (gca, 'XScale','log');
    set (gca, 'YScale','log');
    
end

%% Test for Goodwin. Trigger oscillation by crossing the kd threshold, phase diagram illustration 
close all;
clear;
parameter_Goodwin;

%plot Hill function and identity for z=0
figure(1); 
subplot(2,1,1);
AA = 10.^(-2:.001:4);
F_Aa = pA_T*((beta1*AA.^m + beta2*ka^m)./(deltaA*(AA.^m + ka^m)));
plot (F_Aa, AA,'k--'); hold on;
plot(AA, AA, 'k');
title('Hill Function');
set (gca, 'XScale','log');
set (gca, 'YScale','log');
set (gca, 'XLim', [.5 1e2]);
xlabel('Repressor Concentration');
ylabel('Repressor Expression');

%simulate z=0 to plot phase diagram
%time step
t0 = 0;
tf = 50;
dt = .01;
tt = t0:dt:tf;

z = 0;
kn = deltaN*z;
delay = .5;

sol = dde23(@GoodwinOscNB, delay, @history3, tt);
td = tt;
yd = deval(sol,tt);
ydd(:,tt>delay) = deval(sol, tt(tt>delay)-delay);

plot(ydd(1,:), yd(1,:),'r');
xlabel('Repressor (t-\tau)');
ylabel('Repressor (t)');
set (gca, 'XLim', [1e-2 1e4]);
set (gca, 'YLim', [1e-2 1e4]);
set (gca, 'XScale','log');
set (gca, 'YScale','log');

%plot nullcline to illustrate the point for z = 30
z = 30;
kn = deltaN*z;
subplot(2,1,2);
AA = 10.^(-2:.001:4);
AA2 = (1+koff/deltaN)*kn./(deltaN+kon*AA+koff);
AA3 = kn/deltaN - AA2;
F_Aa = pA_T*((beta1*AA.^m + beta2*ka^m)./(AA.^m + ka^m));

hill = (F_Aa + koff*AA3)./(deltaA+kon*AA2);

plot (hill, AA,'k--'); hold on;
plot(AA, AA, 'k');
title('Hill Function');
set (gca, 'XScale','log');
set (gca, 'YScale','log');
set (gca, 'XLim', [.5 1e2]);
xlabel('Repressor Concentration');
ylabel('Repressor Expression');

%plot simulation for z=30
sol = dde23(@GoodwinOscNB, delay, @history3, tt);
td = tt;
yd = deval(sol,tt);
ydd(:,tt>delay) = deval(sol, tt(tt>delay)-delay);

plot(ydd(1,:), yd(1,:),'r');
xlabel('Repressor (t-\tau)');
ylabel('Repressor (t)');
set (gca, 'XLim', [1e-2 1e4]);
set (gca, 'YLim', [1e-2 1e4]);
set (gca, 'XScale','log');
set (gca, 'YScale','log');


%% Test for Goodwin. Trigger oscillation by crossing the kd threshold, trajectory illustration
 
close all;
clear;
parameter_Goodwin;

%time step
t0 = 0;
tf = 50;
dt = .01;
tt = t0:dt:tf;

%initial conditions
zvec = 0:10:30;

for i =1:length(zvec)
    z = zvec(i); 
    kn = deltaN*z;
    delay = .5;

    sol = dde23(@GoodwinOscNB, delay, @history3, tt);
    
    %calculating equilibrium
    AA = 10.^(-2:.001:4);
    AA2 = (1+koff/deltaN)*kn./(deltaN+kon*AA+koff);
    AA3 = kn/deltaN - AA2;
    F_Aa = pA_T*((beta1*AA.^m + beta2*ka^m)./(AA.^m + ka^m));

    rr = (F_Aa + koff*AA3)./(deltaA+kon*AA2)-AA;
    zerocrossing = rr(1:end-1).*rr(2:end)<0;
    equil = AA(zerocrossing);
    
    
    td = tt;
    yd = deval(sol,tt);
    ydd(:,tt>delay) = deval(sol, tt(tt>delay)-delay);
        
    
     figure(1); %hold off;
     subplot(4,1,i); 
     plot(td, yd(1,:),'b'); hold on;
     plot(td, ydd(1,:),'r--'); hold on;
     title(sprintf('%g nM  Nanobody ',z));
     ylabel('Repressor (nM)');
     plot(td, ones(1,length(td))*equil,'k--'); hold on;
     xlim([0 30]);    
    
end
     figure (1); 
     subplot(4,1,4); xlabel('Time');
     legend('Immature Repressor', 'Mature Repressor','Equilibrium');

%% Generate Amplitude Chart x Nanobody concentration 

clear all; 
close all; 
parameter_Goodwin;

%time step
t0 = 0;
tf = 150;
dt = .01;
tt = t0:dt:tf;
tl = 100;

zvec = 0:1:100;

for i =1:length(zvec)

    z = zvec(i); 
    kn = deltaN*z;
    delay = .5;

    sol = dde23(@GoodwinOscNB, delay, @history3, tt);
    
    disp(i);
    
    td = tt;
    yd = deval(sol,tt);
    
    figure(1);  plot(td, yd(1,:)); ylim([0 100]);
    
    repAmp(i) = range(yd(1,td>tl));
    nanobodyAmp(i) = range(yd(2,td>tl));
    
    ttd = td(td>tl);
    osci_per = yd(1,td>tl)-mean(yd(1,td>tl));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
        
end

figure(2); subplot(221);  plot(zvec, repAmp); 
ylabel('nM');title ('Repressor Amplitude');
subplot(222); plot (zvec, period); title ('Oscillation Period'); 
ylabel('hr'); xlabel('Total Nanobody (nM)');


