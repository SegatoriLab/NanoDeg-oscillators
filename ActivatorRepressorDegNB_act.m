
function P_dot = ActivatorRepressorDegNB_act(~,sol)
% times in hours
% Activator = P_dot(1)
% Repressor = P_dot(2)
% Free nanobody = P_dot(3)
% Nanobody-Activator Complex = P_dot(4)


global pA_T pB_T m n ka kb k5 k6 k1 k2 k3;
global k4 deltaA deltaB kon koff;   

global deltaN kn;


F_Ba = pB_T*((k5*sol(1)^m + k6*ka^m)/(sol(1)^m + ka^m));
F_Aab = (pA_T)*( (k1*ka^m*kb^n + k2*kb^n*sol(1)^m + k3*ka^m*(sol(2))^n + ...
    k4*(sol(2))^n*sol(1)^m)/((kb^n+(sol(2))^n)*(ka^m+sol(1)^m)) );

P_dot = zeros(4,1); 

P_dot(1,1) = F_Aab - deltaA*sol(1) - kon*sol(1)*sol(3) + koff *sol(4) ;
P_dot(2,1) = F_Ba - deltaB*sol(2) ;
P_dot(3,1) = kn - deltaN*sol(3) - kon*sol(1)*sol(3) + koff*sol(4);
P_dot(4,1) = -koff*sol(4) - deltaN*sol(4) + kon*sol(1)*sol(3);

end






