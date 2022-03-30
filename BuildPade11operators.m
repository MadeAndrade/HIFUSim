%% Authored by Joshua Soneson 2018
function[P1,P2] = BuildPade11operators(A,kk,dz,k,JJ)

I = speye(JJ);
kkk = k*kk;
A = A/kkk/kkk;
s = 1i*kkk*dz;
P1 = I + 0.25*(1-s)*A;
P2 = I + 0.25*(1+s)*A;
