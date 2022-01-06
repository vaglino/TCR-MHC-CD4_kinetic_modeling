% syms u(t) v(t) a b
% % Define the equations using == and represent differentiation using the diff function.
% 
% ode1 = diff(u) == 3*u + 4*v;
% ode2 = diff(v) == -4*u + 3*v;
% odes = [ode1; ode2]
% 
% S = dsolve(odes)
% 
% S.u
% S.v
% 
% cond1 = u(0) == a;
% cond2 = v(0) == b;
% conds = [cond1; cond2];
% [uSol(t), vSol(t)] = dsolve(odes,conds)

syms R(t) A(t) kf ks ka k_a uR uA

ode1 = diff(R) == -kf*R(t) - ka*R(t) + k_a*A(t);
ode2 = diff(A) == -ks*A(t) + ka*R(t) - k_a*A(t);
odes = [ode1; ode2]
cond1 = R(0) == uR;
cond2 = A(0) == uA;
% cond2 = A(0) == 1 - uR;
conds = [cond1; cond2];
[RSol(t), ASol(t)] = dsolve(odes,conds)


syms RX(t) k2 k_2 uRX

% ode1 = diff(R) == -kf*R(t) - ka*R(t) + k_a*A(t)  + k_2*RX(t);
% ode2 = diff(A) == -ks*A(t) + ka*R(t) - k_a*A(t) - k2*A(t);
% ode3 = diff(RX) == + k2*A(t) - k_2*RX(t);

% ode1 = diff(R) == -kf*R(t) - ka*R(t) + k_a*A(t) - k2*R(t) + k_2*RX(t);
% ode2 = diff(A) == -ks*A(t) + ka*R(t) - k_a*A(t);
% ode3 = diff(RX) == + k2*R(t) - k_2*RX(t);

ode1 = diff(R) == -kf*R(t) - ka*R(t) + k_a*A(t) + k_2*RX(t);
ode2 = diff(A) == -ks*A(t) + ka*R(t) - k_a*A(t) ;
ode3 = diff(RX) == - k_2*RX(t);
 
odes = [ode1; ode2; ode3]
cond1 = R(0) == uR;
cond2 = A(0) == 0;
% cond3 = RX(0) == uRX;
cond3 = RX(0) == 1 - uR - uA;

conds = [cond1; cond2; cond3];
[RSol(t), ASol(t), RXSol(t)] = dsolve(odes,conds)
% 

R = simplify(RSol(t))
A = simplify(ASol(t))
RX = simplify(RXSol(t))

pdfR = - diff(RSol(t),t)
pdfA = - diff(ASol(t),t)
pdfRX = - diff(RXSol(t),t)

% simplify(simplify(pdfR + pdfA))
simplify(simplify(pdfR + pdfA + pdfRX))


