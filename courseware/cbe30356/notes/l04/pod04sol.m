%% Solution to Fin Radiation Problem
% In this script we solve the temperature distribution in a fin where the
% heat transfer from the fin is in the form of radiation.  We assume that
% there is a specified heat flux at the base of the fin.  Because of the
% strong dependence of the radiative heat flux on temperature the solution
% is very dependent on the initial condition - the temperature at z = 0.
% We thus use an iterative solution, progressing to larger values of z to
% achieve convergence.

global zlimpass

zlimpass = 1;
x = 1;
for i = 1:30
    x = fzero('miss',x);
    zlimpass = zlimpass*1.1;
end
initialtemp = x
%% Conclusion
% The final expression converges to the initial condition (temperature) of
% 1.2011, and the temperature of the fin drops to 0.5 at a distance of 2.2,
% both nicely O(1) quantities.  The decrease in temperature slows
% drastically at larger z because of the T^4 dependence.  At this position
% the dimensionless heat flux is only -0.11, which means that 89% of the
% total heat loss occurs for z less than 2.2.  Thus, the rest of the fin is
% of little use: the optimal fin length should be about 1 - comprising 74%
% of the total heat loss.  If we include back radiation (which becomes
% increasingly significant as the fin temperature drops) this further
% shortens the optimal length of the fin.

%% Comparison to Exact Solution
% It turns out that you -can- solve this problem analytically: multiplying
% both sides by dT/dz you can get a perfect differential on both sides of
% the equation.  Integrating and applying the boundary conditions you get a
% simple analytic solution.  The deviation from the numerical solution is
% very small for small z, with the temperature deviating slightly at large
% z. We can add this to our figure:

T = @(z) (3*z/10^.5 + (2/5)^(3/10)).^(-2/3)

dTdz = @(z) -(2/5)^.5*T(z).^2.5

z = [0:.01:zlimpass];

figure(1)
hold on
plot(z,T(z),'k',z,dTdz(z),'g')
hold off
axis([0 10 -1 1.5])
legend('numerical solution','heat flux','exact temp','exact flux')
%% Miss.m function
% The function called by the program is given below.  Uncomment it and save
% it in a file named miss.m

% function out = miss(x)
% %This function takes in the initial value of the temperature at z = 0,
% %performs the integration to zlimpass (a global variable) and returns the
% %value of the derivative at zlimpass.  We then use this with a rootfinder
% %to get the correct value for x.
% 
% global zlimpass
% 
% fdot = @(z,f) [sign(f(1))*f(2);f(1)^4]; %The derivatives
% % Note that we have added a little "fix" to fdot, as things go haywire if
% % the temperature becomes negative as can happen if you have the wrong
% % initial condition.  Flipping the sign of the derivative forces stability
% % of the differential equation for these conditions: it is just a numerical
% % fix dealing with errors caused by the incorrect IC guess.
% 
% [zout,fout] = ode23(fdot,[0 zlimpass],[x,-1]);
% out = fout(end,2); % we require the temperature derivative to be zero at zlimpass.
% 
% % We add in a little graphics to see how we are doing.
% figure(1)
% plot(zout,fout)
% legend('temperature','heat flux')
% xlabel('z')
% ylabel('temperature')
% grid on
% zoom on
% drawnow

