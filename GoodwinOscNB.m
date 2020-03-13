function P_dot = GoodwinOscNB(~,sol, lag)
% degradation rate in hours
% A = sol(1)
% A(t-tau) = lag(1)
% Free nanobody = sol(2)
% nanobody - A(t-tau) complex = sol(3)


global pA_T beta1 beta2  m  ka;
global deltaA kon koff;   

global deltaN kn;

F_Aa = pA_T*((beta1*lag(1)^m + beta2*ka^m)/(lag(1)^m + ka^m));

P_dot = zeros(3,1); 

P_dot(1,1) = F_Aa - deltaA*sol(1) - kon*lag(1)*sol(2) + koff*sol(3);

P_dot(2,1) = kn - deltaN*sol(2) - kon*lag(1)*sol(2) + koff*sol(3);

P_dot(3,1) = -koff*sol(3) - deltaN*sol(3) + kon*lag(1)*sol(2);

end






