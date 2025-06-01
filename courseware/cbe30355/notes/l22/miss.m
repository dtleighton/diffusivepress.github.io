function out=miss(x)
% We define the derivatives (system of first order differential equations).
% Note that they aren't a function of eta in this case, but we add that in
% as integrators usually require it.  We tack on the displacement thickness
% integral as an extra equation.
fdot = @(eta,f) [f(2);f(3);-f(1).*f(3);1-f(2)];

% We have the initial conditions, taking the input variable of the function
% as the unknown initial condition (fpp, or the wall shear stress).
f0 = [0,0,x,0];

% We use the canned integrator ode23 to do the integration.  We specify a
% range of integration from eta = 0 to 10 (rather than infinity) because
% the integrals converge exponentially outside the boundary layer
% (truncation works fine for this equation).
[etaout,fout] = ode23(fdot,[0 10],f0);

% Now we plot it up to see how it does.  Note that the graphics chew up all
% the CPU time and would be commented out to run faster.
figure(1)
plot(etaout,fout(:,[2,4]))
xlabel('eta')
ylabel('velocity, delta*')
grid on
legend('velocity','delta*')
drawnow

% And we spit out the two things we need: the amount we miss the boundary
% condition at "infinity" by (e.g., our miss), and the displacement
% thickness.
out = fout(end,2)-1
delta = fout(end,4)

% The function would be called using a root finder and an initial guess for
% fpp(0), e.g., "fpp = fzero('miss',1)".