global pA_T  m ka beta1 beta2;
global deltaA  kon koff kd;  

global deltaN kn z;

beta1 = 1.8;
beta2 = 181;

pA_T = 1;
m = 2;

ka = 3; 

deltaA = log(2)/11;  %.9; %2.2 needs to increase

deltaN = log(2)/.9;
kn = 5;
z = kn/deltaN;

kon = 2.7648;
koff = .6264;
kd = koff/kon;
