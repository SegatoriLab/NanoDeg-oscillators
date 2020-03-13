global k1 k2 k3 k4 k5 k6  ka kb kc
global kna knb knc dna dnb dnc S_A S_B S_C kon koff

global PA_T PB_T PC_T PO_T PNB_A PNB_B PNB_C

global m r n

PNB_A = 1;
PNB_B = 1;
PNB_C = 1;

k1 = 1.8;
k2 = 181; 
k3 = k1;
k4 = k2;
k5 = k1;
k6 = k2;

m = 4;
r = 4; 
n = 4;

ka = 3;
kb = 3;
kc = 3;

kna = 1; 
knb = 1;
knc = 1;

%degradation rate
S_A = log(2)/11;
S_B = log(2)/11;
S_C = log(2)/11;

dna = log(2)/.9;
dnb = log(2)/.9;
dnc = log(2)/.9;

%total promoters
PA_T = 1;
PB_T = 1;
PC_T = 1;
PO_T = 1;

kon = 2.7648;
koff = .6264;
