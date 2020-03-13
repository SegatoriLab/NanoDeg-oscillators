global pA_T pB_T pO_T m n ka kb k5 k6 k1 k2 k3 alpha1 alpha2 beta1 beta2;
global k4 k7 k8 deltaA deltaB kon koff kd;  

global deltaN kn z;

alpha1 = 112.5;
alpha2 = 1;
beta1 = .04;
beta2 = 1.8;

k1 = alpha2*beta2; % 1.8; %36;  %1.8
k2 = alpha1*beta2; %200; 
k3 = alpha2*beta1; %0.04; %1.8; %0.05
k4 = alpha1*beta1; %4.5; %1.8; %0.05

pA_T = 1;
m = 2;
n = 2;
ka = 3; 
kb = 3;
k5 = 36; %217;  %36
k6 = .05; %36;   %.05
pB_T = 1;
k7 = 217;
k8 = 1.8; %36;  %1.8
pO_T =1;

deltaA = log(2)/4;  %.9; %2.2 needs to increase
deltaB = log(2)/4;  % 11

deltaN = log(2)/.9;
kn = 5;
z = kn/deltaN;

kon = 2.7648;
koff = .6264;
kd = koff/kon;
