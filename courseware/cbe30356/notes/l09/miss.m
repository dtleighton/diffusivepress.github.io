function out = miss(x)
% This function takes in a guess for the derivative of the temperature at
% y = 0 for an impulsively heated semi-infinite domain.

fdot = @(eta,f) [f(2); -0.5*eta*f(2)]; %The differential equation

f0 = [1,x]; %The initial value

[etaout fout] = ode45(fdot,[0 10],f0);

out = fout(end,1);
