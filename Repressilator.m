function P_dot = Repressilator(t,sol,k1,k2,k3,k4,k5,k6,k7,k8,ka,kb,kc)

%degradation rate
S_A = log(2)/11;
S_B = log(2)/11;
S_C = log(2)/11;
S_O = log(2)/0.633;

%total promoters
PA_T = 10;
PB_T = 1;
PC_T = 1;
PO_T = 1;

%exponent
r = 2;
m = 2;
n = 2;

%{
A = sol(1)
B = sol(2)
C = sol(3)
O = sol(4)
%}

F_Ac = PA_T*((k1*sol(3)^r + k2*kc^r)/(sol(3)^r + kc^r));%Corrected an error here and in the PDF with equations. The arithemetic between terms in num and denom should be addition. Was subtraction.

F_Ba = PB_T*((k3*sol(1)^m + k4*ka^m)/(sol(1)^m + ka^m));

F_Cb = PC_T*((k5*sol(2)^n + k6*kb^n)/(sol(2)^n + kb^n));

F_Oa = PO_T*((k7*sol(1)^m + k8*ka^m)/(sol(1)^m + ka^m));

P_dot = zeros(4,1); 

P_dot(1,1) = F_Ac - S_A*sol(1);

P_dot(2,1) = F_Ba - S_B*sol(2);

P_dot(3,1) = F_Cb - S_C*sol(3);

P_dot(4,1) = F_Oa - S_O*sol(4);

end