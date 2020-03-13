
function P_dot = GoodwinOsc(~,sol, lag)
%degradation rate

global pA_T beta1 beta2  m  ka;
global deltaA;   


%{
A = sol(1)
A(t-tau) = lag(1)
%}

F_Aa = pA_T*((beta1*lag(1)^m + beta2*ka^m)/(lag(1)^m + ka^m));

P_dot = zeros(1,1); 

P_dot(1,1) = F_Aa - deltaA*lag(1);


end






