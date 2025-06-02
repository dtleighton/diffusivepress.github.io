function out = miss2(x)
% This function takes in a guess for the derivative of the temperature at
% y = 0 for an impulsively heated semi-infinite domain, but with heat flux
% at the wall this time.

fdot = @(eta,f) [f(2); 0.5*(f(1)-eta*f(2))]; %The differential equation

f0 = [x,-1]; %The initial value

[etaout fout] = ode45(fdot,[0 10],f0);

out = fout(end,1);
