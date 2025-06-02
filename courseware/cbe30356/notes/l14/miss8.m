function out=miss8(x)
% This function takes in the unknown wall shear stress and heat flux and
% returns zero if the boundary conditions at infinity are satisfied.
% Because the integral is unstable for large eta unless the correct initial
% conditions are specified, we "sneak up" on the solution.

% call function with:
%
% global etalim;etalim=2;x=[1;-1];
%
% for i=1:5;x=fminsearch('miss8',x);etalim=etalim*1.5;end;x

global etalim

y0 = [0 x(1) 0 1 x(2)]'; %The initial condition

pr = .7; %The prandtl number

ydot = @(eta,y) [y(2);-y(4)+1/pr*((y(1)^2-eta*y(1)*y(2))/2+y(3)*y(2));...
    -(y(1)-eta*y(2))/2;y(5);y(5)*(y(3)-eta*y(1)/4)];

[etaout yout] = ode23(ydot,[0 etalim],y0);

figure(1)
plot(etaout,yout)
xlabel('eta')
ylabel('function values')
legend('f','fp','g','h','hp')
grid on
drawnow

out = norm([yout(end,1),yout(end,4)]);