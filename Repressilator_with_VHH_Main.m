%% Generate activation example to plot the 3 trajectories

close all;
clear;

parameter_REP;

%time step
t0 = 0;
tf = 200;
dt = .01;
t = t0:dt:tf;

%initial conditions
pvec = [0 10];


for i =1:length(pvec)
    
    PNB = pvec(i);
 
    inicondd = zeros (1,4);

    
    [td, xd] = ode23s(@(t, x) Repressilator_NB_node(t,x), t, inicondd);
       
    figure(2); hold off;
    subplot(2,1,i);
    plot (td, xd(:,1)); hold all;
    plot(td, xd(:,2));
    if i == 1
        legend('Repressor 1','Repressor 2');
    else
        plot(td, xd(:,3)+xd(:,4));
        legend('Repressor 1','Repressor 2','Total Nanobody');
    end
    set (gca, 'YScale','log');
    ylim([1e-1 1e3]);

    xlabel('Time');
    ylabel('nM');
    

   
    
    
end

%% Generate Amplitude Chart x Nanobody degradation rate for Repressilator

close all;
clear;

parameter_REP;

%time step
t0 = 0;
tf = 300;
dt = .01;
t = t0:dt:tf;
tl = 250;


hl = 10.^[-1:.1:3];
dnvec = log(2)./hl;

for i =1:length(hl)
    
    dna = dnvec(i);
  
    PNB = 10;
    
    inicondd = 1e-3*zeros (1,4);
    
    [td, xd] = ode23s(@(t, x) Repressilator_NB_node(t,x), t, inicondd);
    
    R1(i) = range(xd(td>tl,1));
    R2(i) = range(xd(td>tl,2));
    R3(i) = range(xd(td>tl,3));
    R4(i) = range(xd(td>tl,4));
    R5(i) = range(xd(td>tl,3)+xd(td>tl,4));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-mean(xd(td>tl,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
    
    disp(hl(i)); 

    
    figure(2); hold off;
    subplot(2,2,1);
    plot (td, xd(:,1)); %hold all;
    set (gca, 'YScale','log');
    title('Repressor 1');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,2,2);
    plot (td, xd(:,2)); %hold all;
    set (gca, 'YScale','log');
    title('Repressor 2');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,2,3);
    plot (td, xd(:,3)); %hold all;
    title('Free Nanobody');
    set (gca, 'YScale','log');

    
    figure(2);
    subplot(2,2,4);
    plot(xd(:,1),xd(:,2)); %hold all;
    xlabel('Repressor 1'); ylabel('Repressor 2');
    title ('Phase Plot');
    set (gca, 'YScale','log');
    set (gca, 'XScale','log');
    
    
end

figure(4); 
subplot(221);  plot(hl, R1); title ('Repressor 1 Amplitude');
ylabel('nM'); xlim([1e-1 11]); xlabel('Half-life (hrs)'); set (gca, 'XScale','log');
subplot(222); plot (hl, R2); title ('Repressor 2 Amplitude'); 
ylabel('nM'); xlim([1e-1 11]); xlabel('Half-life (hrs)');set (gca, 'XScale','log');
subplot(223); plot (hl, R3); title ('Free Nanobody Amplitude'); 
ylabel('nM'); xlim([1e-1 11]); xlabel('Half-life (hrs)');set (gca, 'XScale','log');
subplot(224); plot (hl, period); title ('Oscillation Period'); 
ylabel('nM'); xlim([1e-1 11]); xlabel('Half-life (hrs)');set (gca, 'XScale','log');



%% Generate Amplitude Chart x Nanobody concentration for Repressilator

% Repressor 2 maximum expression rate should be smaller than repressor 1
% maximum expression rate

close all;
clear;

parameter_REP;

%time step
t0 = 0;
tf = 300;
dt = .01;
t = t0:dt:tf;

%initial conditions
pvec = [10.^(-1:.05:4)];


for i =1:length(pvec)
    
    PNB = pvec(i);

    if PNB < 3 || PNB > 300
        t0 = 0;
        tf = 1000;
        dt = .01;
        t = t0:dt:tf;
        tl = 800;
    else
        t0 = 0;
        tf = 100;
        dt = .01;
        t = t0:dt:tf;
        tl = 50;
    end
    
    inicondd = 1e-3*ones (1,4);
    inicondd(1) = 10;
    inicondd(2) = 5;
    
    [td, xd] = ode23s(@(t, x) Repressilator_NB_node(t,x), t, inicondd);
    
    R1(i) = range(xd(td>tl,1));
    R2(i) = range(xd(td>tl,2));
    R3(i) = range(xd(td>tl,3));
    R4(i) = range(xd(td>tl,4));
    R5(i) = range(xd(td>tl,3)+xd(td>tl,4));
    
    ttd = td(td>50);
    osci_per = xd(td>tl,1)-mean(xd(td>tl,1));
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period(i) = mean(timesofcross(3:end)- timesofcross(1:end-2));
    
    disp(i); 

    
    figure(2); hold off;
    subplot(2,2,1);
    plot (td, xd(:,1)); %hold all;
    set (gca, 'YScale','log');
    title('Repressor 1');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,2,2);
    plot (td, xd(:,2)); %hold all;
    set (gca, 'YScale','log');
    title('Repressor 2');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,2,3);
    plot (td, xd(:,3)); %hold all;
    title('Free Nanobody');
    set (gca, 'YScale','log');

    
    figure(2);
    subplot(2,2,4);
    plot(xd(:,1),xd(:,2)); %hold all;
    xlabel('Repressor 1'); ylabel('Repressor 2');
    title ('Phase Plot');
    set (gca, 'YScale','log');
    set (gca, 'XScale','log');
    
    
end

figure(4); subplot(221);  plot(pvec, R1); 
ylabel('nM');title ('Repressor 1 Amplitude');set (gca, 'XScale','log');
subplot(222); plot (pvec, R2); title ('Repressor 2 Amplitude'); 
ylabel('nM'); xlabel('Nanobody DNA (nM)');set (gca, 'XScale','log');
subplot(223); plot (pvec, R3); title ('Free Nanobody Amplitude'); 
ylabel('nM'); xlabel('Nanobody DNA (nM)');set (gca, 'XScale','log');
subplot(224); plot (pvec, period); title ('Oscillation Period'); 
ylabel('hr'); xlabel('Nanobody DNA (nM)');set (gca, 'XScale','log');



%% Generate Amplitude Chart x sensitivity of hill function

% Repressor 2 maximum expression rate should be smaller than repressor 1
% maximum expression rate

close all;
clear;

parameter_REP;

%time step
t0 = 0;
tf = 300;
dt = .01;
t = t0:dt:tf;
tl = 250;

PNB = 10;

%initial conditions
hillvecm = 2.9:.1:10;
hillvecn = 3:.1:10;


for i =1:length(hillvecm)
 
  for j =1:length(hillvecn)
    
    m = hillvecm(i);
    n = hillvecn(j);
    
    inicondd = 1e-3*ones (1,4);
    inicondd(1) = 10;
    inicondd(2) = 5;
    
    [td, xd] = ode23s(@(t, x) Repressilator_NB_node(t,x), t, inicondd);
    
    R1(i,j) = range(xd(td>tl,1));
    R2(i,j) = range(xd(td>tl,2));
    R3(i,j) = range(xd(td>tl,3));
    R4(i,j) = range(xd(td>tl,4));
    R5(i,j) = range(xd(td>tl,3)+xd(td>tl,4));
    
    ttd = td(td>50);
    osci_per = xd(td>50,1)-10;
    zerocross = osci_per(1:end-1).*osci_per(2:end) < 0;
    timesofcross = ttd(zerocross);
    period(i,j) = mean(timesofcross(3:end)- timesofcross(1:end-2));
    
    fprintf("m=%d n=%d\n", m, n); 

    
    figure(2); hold off;
    subplot(2,2,1);
    plot (td, xd(:,1)); %hold all;
    set (gca, 'YScale','log');
    title('Repressor 1');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,2,2);
    plot (td, xd(:,2)); %hold all;
    set (gca, 'YScale','log');
    title('Repressor 2');
    xlabel('Time');
    ylabel('nM');
    
    figure(2);
    subplot(2,2,3);
    plot (td, xd(:,3)); %hold all;
    title('Free Nanobody');
    set (gca, 'YScale','log');

    
    figure(2);
    subplot(2,2,4);
    plot(xd(:,1),xd(:,2)); %hold all;
    xlabel('Repressor 1'); ylabel('Repressor 2');
    title ('Phase Plot');
    set (gca, 'YScale','log');
    set (gca, 'XScale','log');
    
  end   
end

figure(4); subplot(211);  plot(R1(:,5:9)); 
ylabel('nM');title ('Repressor 1 Amplitude');
subplot(212); plot (period(:,5:9)); title ('Oscillation Period'); 
ylabel('hr'); xlabel('Hill Coefficient of Repressor 1 (m)');
for i = 6:10
    leg{i-5} =  sprintf("n=%d",i);
end
legend(leg);

RR1 = R1;
RR1 (R1<.1) = 0;
pperiod = period;
pperiod (R1<.1) = NaN;

figure; h1 = heatmap(hillvecn(11:end),hillvecm(end:-1:12), pperiod(12:end,end:-1:11)');
h1.MissingDataLabel = 'No oscillation';
h1.MissingDataColor = 'w';
h1.GridVisible = 'off';
figure; h2 = heatmap(hillvecn(11:end),hillvecm(end:-1:12),RR1(12:end,end:-1:11)');
h2.GridVisible = 'off';


   