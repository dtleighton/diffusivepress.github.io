function out=fsmiss(x)
% This function integrates the Falkner-Skan equation and returns zero when
% the boundary condition at infinity is satisfied.  The unknown argument is
% the stress at the wall f''(0).  It can be called from the command line
% using the function minimization routine fminsearch('fsmiss',1).  We use a
% combined output including both the "miss" at infinity as well as the
% second derivative (which should vanish at infinity as well).  This
% provides a more robust convergence for larger values of beta.

beta = .5; % We pick a value of beta corresponding to a 90° wedge.

fdot = @(t,f) [f(2); f(3); -f(1)*f(3) - beta*(1-f(2)^2)];

f0 = [0;0;x]; %The initial condition

[etaout fout] = ode23(fdot,[0 8],f0);

figure(1)
plot(etaout,fout(:,2))
xlabel('eta')
ylabel('scaled velocity')
drawnow

out = norm([fout(end,2)-1,fout(end,3)]);