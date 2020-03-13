function P_dot = Repressilator_NB_node(~,sol)

global k1 k2 k3 k4  ka kb 
global  dna   S_A S_B  kon koff
global PA_T PB_T  PNB

global m n 

%{
Repressor A = sol(1)
Repressor B = sol(2)
Free Nanobody = sol(3)
Nanobody complex = sol(4)
%}

F_Ba = PB_T*((k1*sol(1)^m + k2*ka^m)/(sol(1)^m + ka^m));

F_NBb = PNB*((k3*sol(2)^n + k4*kb^n)/(sol(2)^n + kb^n));

P_dot = zeros(4,1); 

P_dot(1,1) = PA_T*k2 - S_A*sol(1) - kon*sol(1)*sol(3) + koff*sol(4);

P_dot(2,1) = F_Ba - S_B*sol(2);  

P_dot(3,1) = F_NBb - dna*sol(3) - kon*sol(1)*sol(3) + koff*sol(4);

P_dot(4,1) = -koff*sol(4) - dna*sol(4) + kon*sol(1)*sol(3);


end