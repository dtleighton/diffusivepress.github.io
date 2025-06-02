%% Heated Semi infinite slab: Constant Wall Heat Flux
%
% The miss2 program is:
%
% function out = miss2(x)
% % This function takes in a guess for the derivative of the temperature at
% % y = 0 for an impulsively heated semi-infinite domain, but with heat flux
% % at the wall this time.
% 
% fdot = @(eta,f) [f(2); 0.5*(f(1)-eta*f(2))]; %The differential equation
% 
% f0 = [x,-1]; %The initial value
% 
% [etaout fout] = ode45(fdot,[0 10],f0);
% 
% out = fout(end,1);

x = 1; % Our initial guess for the temperature at the wall

x = fzero('miss2',x) % Our solution!

%% Plotting things up
% We just cut and paste from the miss.m routine to get the profile:

fdot = @(eta,f) [f(2); 0.5*(f(1)-eta*f(2))]; %The differential equation

f0 = [x,-1]; %The initial value

[etaout fout] = ode45(fdot,[0 10],f0);

out = fout(end,1);

figure(1)
plot(etaout,fout(:,1))
xlabel('eta')
ylabel('dimensionless temperature')
title('Temperature Distribution for Constant Heat Flux')
axis([0 5 0 1.5])
grid on

