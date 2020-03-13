%% Generate Amplitude Chart x Nanobody concentration for Repressilator with independent nanobodies

close all;
clear;

parameter_REP;

%time step
t0 = 0;
tf = 500;
dt = .01;
t = t0:dt:tf;

%initial conditions
zvec = 10.^(0:.05:3);


for i =1:length(zvec)
    
    kna = zvec(i)*dna;
    knb = zvec(i)*dnb;
    knc = zvec(i)*dnc;
    
    inicondd = 1e-3*ones (1,9);
    inicondd(1) = 10;
    inicondd(2) = 5;
    
    [td, xd] = ode23s(@(t, x) Repressilator_NB(t,x), t, inicondd);
    
    R1(i) = range(xd(td>50,1));
    DR (i) = max(xd(td>50,1))/min(xd(td>50,1));
    
    R2(i) = range(xd(td>50,2));
    R3(i) = range(xd(td>50,3));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-mean(xd(td>50,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
    
    disp(zvec(i)); 

    
    figure(2); 
    subplot(2,2,1);
    plot (td, xd(:,1)); 
    set (gca, 'YScale','log');
    title('Repressor 1');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,2,2);
    plot (td, xd(:,2)); hold all;
    set (gca, 'YScale','log');
    title('Repressor 2');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,2,3);
    plot (td, xd(:,3)); hold all;
    set (gca, 'YScale','log');
    title('Repressor 3');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,2,4); 
    plot(xd(:,1),xd(:,2)); hold all;
    xlabel('Repressor 1'); ylabel('Repressor 2');
    set (gca, 'YScale','log');
    set (gca, 'XScale','log');
    
    
end

figure(4); subplot(221);  plot(zvec, R1); 
ylabel('nM');title ('Repressor 1 Amplitude'); set(gca,'XScale','log');
subplot(222); plot (zvec, R2); title ('Repressor 2 Amplitude'); 
ylabel('nM'); xlabel('Nanobody DNA (nM)');set(gca,'XScale','log');
subplot(223); plot (zvec, R3); title ('Repressor 3 Amplitude'); 
ylabel('nM'); xlabel('Nanobody DNA (nM)');set(gca,'XScale','log');
subplot(224); plot (zvec, period); title ('Oscillation Period'); 
ylabel('hr'); xlabel('Nanobody DNA (nM)');set(gca,'XScale','log');

%% Generate sample trajectories for the common Nanobody

close all;
clear;

parameter_REP

%time step
t0 = 0;
tf = 400;
dt = .01;
t = t0:dt:tf;

%initial conditions
zvec = [0 30];


for i =1:length(zvec)
    
    z = zvec(i);
    kna = z*dna;
    
    inicondd = 1e-3*ones (1,7);
    inicondd(1) = 1;
    inicondd(2) = 10;
    inicondd(3) = 100;
    inicondd(4) = z;
    
    [td, xd] = ode23s(@(t, x) Repressilator_common_NB(t,x), t, inicondd);
    
    figure(1); 
    subplot(2,1,i);
    plot (td, xd(:,1:4));
    set (gca, 'YScale','log');
    legend('Repressor 1','Repressor 2', 'Repressor 3', 'Free Nanobody');
    xlabel('Time');
    ylabel('nM');
    xlim([50 100]);
   
    figure(2);
    plot(xd(:,1),xd(:,4)); hold all;
    set (gca, 'XScale','log'); set (gca,'YScale','log');
    xlabel('Repressor 1'); ylabel('Free Nanobody');
 
end

figure(1); subplot(211); title("No nanobody"); subplot(212); title("Common nanobody");

%% Generate sample trajectories for the common Nanobody, different rates

close all;
clear;

parameter_REP

%time step
t0 = 0;
tf = 400;
dt = .01;
t = t0:dt:tf;

%initial conditions
zvec = [0 30];
g1 = 1.2;
g2 = .8;

k1 = k1*g1;
k2 = k2*g1;
k5 = k5*g2;
k6 = k6*g2;

for i =1:length(zvec)
    
    z = zvec(i);
    kna = z*dna;
    
    inicondd = 1e-3*ones (1,7);
    inicondd(1) = 1;
    inicondd(2) = 10;
    inicondd(3) = 100;
    inicondd(4) = z;
    
    [td, xd] = ode23s(@(t, x) Repressilator_common_NB(t,x), t, inicondd);
    
    figure(1); 
    subplot(2,1,i);
    plot (td, xd(:,1:4));
    set (gca, 'YScale','log');
    legend('Repressor 1','Repressor 2', 'Repressor 3', 'Free Nanobody');
    xlabel('Time');
    ylabel('nM');
    xlim([50 100]);
   
    figure(2);
    plot(xd(:,1),xd(:,4)); hold all;
    set (gca, 'XScale','log'); set (gca,'YScale','log');
    xlabel('Repressor 1'); ylabel('Free Nanobody');
 
end

figure(1); subplot(211); title("No nanobody"); subplot(212); title("Common nanobody");

%% Generate sample trajectories for the common Nanobody - transient

close all;
clear;

parameter_REP

%time step
t0 = 0;
tf = 400;
dt = .01;
t = t0:dt:tf;

%initial conditions
zvec = [5 10 20];


for i =1:length(zvec)
    
    z = zvec(i);
    kna = z/dna;
    
    disp(kna);
    
    inicondd = 1e-3*ones (1,7);
    inicondd(1) = 1;
    inicondd(2) = 10;
    inicondd(3) = 100;
    inicondd(4) = z;
    
    [td, xd] = ode23s(@(t, x) Repressilator_common_NB(t,x), t, inicondd);
    
    figure(1); 
    subplot(2,1,1);
    plot (td, xd(:,1)); hold all
    set (gca, 'YScale','log');
    title('Repressor 1');
    xlabel('Time');
    ylabel('nM');
   
    figure(1); 
    subplot(2,1,2);
    plot (td, xd(:,4)); hold all;
    set (gca, 'YScale','log');
    title('Free Nanobody');
    xlabel('Time');
    ylabel('nM');
    
 
end

legend('5nM','10nM','20nM');


%% Generate Amplitude Chart x Nanobody concentration for Repressilator with common Nanobody

close all;
clear;

parameter_REP

%time step
t0 = 0;
tf = 400;
dt = .01;
t = t0:dt:tf;
tl = 350;

%initial conditions
zvec = 10.^(1:.05:3);


for i =1:length(zvec)
    
    z = zvec(i);
    kna = z*dna;
   
    if z< 40 || z >100
        %time step
        t0 = 0;
        tf = 500;
        dt = .01;
        t = t0:dt:tf;
        tl = 400;
    else
        %time step
        t0 = 0;
        tf = 200;
        dt = .01;
        t = t0:dt:tf;
        tl = 150;
    end
    
    inicondd = 1e-3*ones (1,7);
    inicondd(1) = 1;
    inicondd(2) = 10;
    inicondd(3) = 100;
    inicondd(4) = z;
    
    [td, xd] = ode23s(@(t, x) Repressilator_common_NB(t,x), t, inicondd);
    
    R1(i) = range(xd(td>tl,1));
    R2(i) = range(xd(td>tl,2));
    R3(i) = range(xd(td>tl,3));
    R4(i) = range(xd(td>tl,4));
    
    ttd = td(td>50);
    osci_per = xd(td>tl,1)-mean (xd(td>tl,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
    
    disp(zvec(i)); 

    
    figure(2); 
    subplot(2,3,1);
    plot (td, xd(:,1)); hold all;
    set (gca, 'YScale','log');
    title('Repressor 1');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,3,2);
    plot (td, xd(:,2)); hold all;
    set (gca, 'YScale','log');
    title('Repressor 2');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,3,3);
    plot (td, xd(:,3)); hold all;
    set (gca, 'YScale','log');
    title('Repressor 3');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,3,4); 
    plot(xd(:,1),xd(:,2)); hold all;
    xlabel('Repressor 1'); ylabel('Repressor 2');
    set (gca, 'YScale','log');
    set (gca, 'XScale','log');
    
    figure(2);
    subplot(2,3,5);
    plot (td, xd(:,3)); hold all;
    set (gca, 'YScale','log');
    title('Free Nanobody');
    xlabel('Time');
    ylabel('nM');
    
    
end

figure(4); subplot(221);  plot(zvec, R1); 
ylabel('nM');title ('Repressor 1 Amplitude');set(gca,'XScale','log');
subplot(222); plot (zvec, R2); title ('Repressor 2 Amplitude'); 
ylabel('nM'); xlabel('Nanobody DNA (nM)');set(gca,'XScale','log');
subplot(223); plot (zvec, R3); title ('Repressor 3 Amplitude'); 
ylabel('nM'); xlabel('Nanobody DNA (nM)'); set(gca,'XScale','log');
subplot(224); plot (zvec, period); title ('Oscillation Period'); 
ylabel('hr'); xlabel('Nanobody DNA (nM)'); set(gca,'XScale','log');

%% Generate Amplitude Chart x Nanobody degradation rate for Repressilator with common Nanobody

close all;
clear;

parameter_REP

%time step
t0 = 0;
tf = 400;
dt = .01;
t = t0:dt:tf;
tl = 250;

%initial conditions
hlvec = 10.^(-1:.05:2);



for i =1:length(hlvec)

    
    
    dna = log(2)/hlvec(i);
    
    z = 30;
    kna = z*dna;
    
%     if z< 20 || z >60
        %time step
        t0 = 0;
        tf = 500;
        dt = .01;
        t = t0:dt:tf;
        tl = 400;
%     else
%         %time step
%         t0 = 0;
%         tf = 200;
%         dt = .01;
%         t = t0:dt:tf;
%         tl = 150;
%     end
    
    inicondd = 1e-3*ones (1,7);
    inicondd(1) = 1;
    inicondd(2) = 10;
    inicondd(3) = 100;
    inicondd(4) = 1;
    
    [td, xd] = ode23s(@(t, x) Repressilator_common_NB(t,x), t, inicondd);
    
    R1(i) = range(xd(td>tl,1));
    R2(i) = range(xd(td>tl,2));
    R3(i) = range(xd(td>tl,3));
    R4(i) = range(xd(td>tl,4));
    
    ttd = td(td>50);
    osci_per = xd(td>tl,1)-mean (xd(td>tl,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
    
    disp(hlvec(i)); 

    
    figure(2); 
    subplot(2,3,1);
    plot (td, xd(:,1)); hold all;
    set (gca, 'YScale','log');
    title('Repressor 1');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,3,2);
    plot (td, xd(:,2)); hold all;
    set (gca, 'YScale','log');
    title('Repressor 2');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,3,3);
    plot (td, xd(:,3)); hold all;
    set (gca, 'YScale','log');
    title('Repressor 3');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,3,4); 
    plot(xd(:,1),xd(:,2)); hold all;
    xlabel('Repressor 1'); ylabel('Repressor 2');
    set (gca, 'YScale','log');
    set (gca, 'XScale','log');
    
    figure(2);
    subplot(2,3,5);
    plot (td, xd(:,3)); hold all;
    set (gca, 'YScale','log');
    title('Free Nanobody');
    xlabel('Time');
    ylabel('nM');
    
    
end

figure(4); subplot(221);  plot(hlvec, R1); 
ylabel('nM');title ('Repressor 1 Amplitude');set(gca,'XScale','log');
subplot(222); plot (hlvec, R2); title ('Repressor 2 Amplitude'); 
ylabel('nM'); xlabel('Nanobody DNA (nM)');set(gca,'XScale','log');
subplot(223); plot (hlvec, R3); title ('Repressor 3 Amplitude'); 
ylabel('nM'); xlabel('Nanobody DNA (nM)'); set(gca,'XScale','log');
subplot(224); plot (hlvec, period); title ('Oscillation Period'); 
ylabel('hr'); xlabel('Nanobody DNA (nM)'); set(gca,'XScale','log');


