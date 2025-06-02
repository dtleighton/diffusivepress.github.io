function out = miss(x)
%This function takes in the initial value of the temperature at z = 0,
%performs the integration to zlimpass (a global variable) and returns the
%value of the derivative at zlimpass.  We then use this with a rootfinder
%to get the correct value for x.

global zlimpass

fdot = @(z,f) [sign(f(1))*f(2);f(1)^4]; %The derivatives
% Note that we have added a little "fix" to fdot, as things go haywire if
% the temperature becomes negative as can happen if you have the wrong
% initial condition.  Flipping the sign of the derivative forces stability
% of the differential equation for these conditions: it is just a numerical
% fix dealing with errors caused by the incorrect IC guess.

[zout,fout] = ode23(fdot,[0 zlimpass],[x,-1]);
out = fout(end,2); % we require the temperature derivative to be zero at zlimpass.

% We add in a little graphics to see how we are doing.
figure(1)
plot(zout,fout)
legend('temperature','heat flux')
xlabel('z')
ylabel('temperature')
grid on
zoom on
drawnow

