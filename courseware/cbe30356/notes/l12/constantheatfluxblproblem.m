%% Nusselt-Graetz Entrance Region Problem
% In this problem we look at the boundary layer solution for the
% Nusselt-Graetz problem for short tubes in the entrance region.  We use
% the constant wall heat flux boundary condition.  The transformed ODE is:
%
% f'' = 1/3 * eta * (f - eta * f')
%
% f'(0) = -1; f(inf) = 0
%
% We want the temperature at the wall f(0), which is our unknown shooting
% parameter.  We solve the differential equation in the function miss4.m
% given below:
% 
% function out = miss4(x)
% % This function takes in a guess for the unknown wall temperature for the
% % constant heat flux entrance length problem and returns the temperature
% % far away.
% 
% fdot = @(eta,f) [f(2); 1/3*eta*(f(1) - eta*f(2))];
% 
% f0 = [x,-1];
% 
% [etaout fout] = ode23(fdot,[0,10],f0);
% 
% out = fout(end,1);

x = fzero('miss4',1)

% We can plot it up too:
fdot = @(eta,f) [f(2); 1/3*eta*(f(1) - eta*f(2))];

f0 = [x,-1];

[etaout fout] = ode23(fdot,[0,10],f0);

figure(1)
plot(etaout,fout(:,1))
xlabel('eta')
ylabel('dimensionless scaled temperature')
title('Plot of Dimensionless Scaled Temperature in Entrance Region')
axis([0 5 0 2])
grid on

%% Conclusion:
% The boundary layer solution converges quite fast, so a limit of 10 is
% more than adequate.  The dimensionless wall temperature is f(0) = 1.5363,
% a nicely O(1) result.