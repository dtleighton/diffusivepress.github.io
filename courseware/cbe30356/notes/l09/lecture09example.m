%% Heated Semi infinite slab: comparison of numerical solution to exact
% In this case we use the shooting method to compare the numerical solution
% for this problem to the exact error function solution.  We use the
% function miss.m to get the degree to which we miss our boundary condition
% at infinity.  Because the solution decreases exponentially as exp(-0.5
% eta^2) we just go out to a value of 10.
%
% The miss program is:
%
% function out = miss(x)
% % This function takes in a guess for the derivative of the temperature at
% % y = 0 for an impulsively heated semi-infinite domain.
% 
% fdot = @(eta,f) [f(2), -0.5*eta*f(2)]; %The differential equation
% 
% f0 = [1,x]; %The initial value
% 
% [etaout fout] = ode45(fdot,[0 10],f0);
% 
% out = fout(end,1);

x = -1; % Our initial guess

x = fzero('miss',x) % Our solution!

%% Plotting things up
% We just cut and paste from the miss.m routine to get the profile:

fdot = @(eta,f) [f(2); -0.5*eta*f(2)]; %The differential equation

f0 = [1,x]; %The initial value

[etaout fout] = ode45(fdot,[0 10],f0);

fexact = 1 - erf(etaout/2); %The exact solution

figure(1)
plot(etaout,fout(:,1),'o',etaout,fexact)
legend('numerical solution','exact solution')
xlabel('eta')
ylabel('dimensionless temperature')
title('Comparison of Exact and Numerical Solutions')
axis([0 5 0 1])

% We can also plot up the deviation:
figure(2)
plot(etaout,fout(:,1)-fexact,'o-')
xlabel('eta')
ylabel('deviation from exact solution')

%% Conclusion
% The deviation between the numerical and exact solutions is less than 5e-5
% for all eta.  This could be further reduced if you set a tighter
% tolerance for fzero (the root finder) or ode45 (the integrator), but it
% illustrates the ease of solving this sort of problem numerically.


