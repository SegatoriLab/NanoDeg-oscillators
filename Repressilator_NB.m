function P_dot = Repressilator_NB(~,sol)

global k1 k2 k3 k4 k5 k6  ka kb kc
global kna knb knc dna dnb dnc S_A S_B S_C kon koff
global PA_T PB_T PC_T  PNB_A PNB_B PNB_C

global m n r

%{
A = sol(1)
B = sol(2)
C = sol(3)
O = sol(4)
%}

F_Ac = PA_T*((k1*sol(3)^r + k2*kc^r)/(sol(3)^r + kc^r));

F_Ba = PB_T*((k3*sol(1)^m + k4*ka^m)/(sol(1)^m + ka^m));

F_Cb = PC_T*((k5*sol(2)^n + k6*kb^n)/(sol(2)^n + kb^n));

P_dot = zeros(9,1); 

P_dot(1,1) = F_Ac - S_A*sol(1) - kon*sol(1)*sol(4) + koff*sol(5);

P_dot(2,1) = F_Ba - S_B*sol(2)  - kon*sol(2)*sol(6) + koff*sol(7);

P_dot(3,1) = F_Cb - S_C*sol(3) - kon*sol(3)*sol(8) + koff*sol(9);

P_dot(4,1) = PNB_A*kna - dna*sol(4) - kon*sol(1)*sol(4) + koff*sol(5);

P_dot(5,1) = -koff*sol(5) - dna*sol(5) + kon*sol(1)*sol(4);

P_dot(6,1) = PNB_B*knb - dnb*sol(6) - kon*sol(2)*sol(6) + koff*sol(7);

P_dot(7,1) = -koff*sol(7) - dnb*sol(7) + kon*sol(2)*sol(6);

P_dot(8,1) = PNB_C*knc - dnc*sol(8) - kon*sol(3)*sol(8) + koff*sol(9);

P_dot(9,1) = -koff*sol(9) - dnc*sol(9) + kon*sol(3)*sol(8);

end