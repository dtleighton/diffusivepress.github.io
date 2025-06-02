function out = miss4(x)
% This function takes in a guess for the unknown wall temperature for the
% constant heat flux entrance length problem and returns the temperature
% far away.

fdot = @(eta,f) [f(2); 1/3*eta*(f(1) - eta*f(2))];

f0 = [x,-1];

[etaout fout] = ode23(fdot,[0,10],f0);

out = fout(end,1);