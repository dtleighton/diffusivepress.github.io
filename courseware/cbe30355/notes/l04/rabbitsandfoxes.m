%% The predator-prey problem
% In this script we show how to analyze the predator-prey problem and show
% how you can simulate the periodicity of biological systems by solving the
% dynamics as a system of first order equations.  The classic example is
% rabbits and foxes.  If there is sufficient food (assumed not to be a
% limiting source) then rabbits will increase in population in proportion
% to their density (exponential growth).  Foxes, on the other hand, will
% eat the rabbits (producing more foxes) and in the absence of rabbits will
% wander off to find food elsewhere.  The population dynamics may be
% approximated by the first order equations:
%
% dr/dt = r*(2-f)
%
% df/dt = (r-1)*f
%
% These equations (known as the the Lotka?Volterra equations) can be set up
% as a set of first order equations (a vector).  We let y(1) be rabbits and
% y(2) be foxes.  We define the time derivative as the vector ydot.  We
% need to specify an initial condition (we take these to be (1;.5) - it is
% fun to play with this, and with the constants in the model as well) and a
% time range for the integration (we do it for a time of 10).  We will
% define ydot(t,y) as an anonymous function (a Matlab utility).  Note that
% it isn't an explicit function of time (for this example) but integrators
% are set up to take t in as one of the arguments.  So:

y0 = [1;.5]; % The initial condition
tspan = [0 10];

ydot = @(t,y) [y(1).*(2-y(2)) ; (y(1)-1).*y(2)];

% And we get the output:
[tout yout] = ode23(ydot,tspan,y0);

figure(1)
plot(tout,yout)
xlabel('time')
ylabel('population')
legend('rabbits','foxes')
title('Predator-Prey Simulation')

%% Conclusion
% As you can see from the graph, the rabbit population initially grows
% because there aren't enough foxes around to eat them up.  After awhile,
% however, the fox population grows in response to the food supply and
% eventually the rabbit population crashes.  The fox population then falls
% in response, and the cycle repeats itself.  You find these sorts of
% cyclical population dynamics in many ecological systems.
%
% The key to this example, of course, is the ease with which you can code
% up a set of ODE's and solve them nummerically.  This is easily done in
% any language (and a canned integrator isn't really necessary either).  It
% really just takes a few lines of code, and you should become adept at
% doing this.