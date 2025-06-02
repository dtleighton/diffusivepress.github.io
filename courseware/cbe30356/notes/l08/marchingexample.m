%% Finite Difference Marching Solution: Temperature Dependent Conductivity
% In this script we demonstrate the simplest of the techniques for
% numerical solutions to non-linear parabolic PDEs: the finite difference
% approach.  We solve the example problem of a heated slab of width 1 with
% a temperature dependent conductivity.  We use a simple explicit forward
% difference in time and a second order finite difference in space.  We
% take the conductivity to be k = k0*(1+c*(T-T0)).

c = 2;

% The analytic solution to the asymptotic temperature is just:
Tinf = @(x) ((1+2*c*(1-x)).^.5-1)/c

% We discretize the domain in space:
n = 50;
dx = 1/n;
x = [0:dx:1]'; %Note: there are n+1 elements

T = zeros(size(x)); %The initial condition

dt = .49/(c+1)/n^2; %The time discretization we need to avoid instability for c>0

t = 0;

tkeep = [0:dt:2]; %The times we want to evaluate the bottom temperature at
Tbottomkeep = zeros(size(tkeep)); %We initialize that array

pt = [2:n]'; %These are the indices of the interior points.

for i = 2:length(tkeep)
    right = (1+c*(T(pt+1)+T(pt))/2).*(T(pt+1)-T(pt))/dx;
    left = (1+c*(T(pt)+T(pt-1))/2).*(T(pt)-T(pt-1))/dx;
    dTdt = (right-left)/dx;
    T(pt) = T(pt)+dTdt*dt;
    
    % and the boundary conditions:
    T(n+1) = 0;
    T(1) = 4/3*T(2) -T(3)/3 + 2/3*dx/(1+c*(2*T(2)-T(3)));
    
    Tbottomkeep(i) = T(1);
    
    % We add a little graphics (comment out for speed!)
    if i/200==floor(i/200) % We plot up every 200th iteration
    figure(1)
    plot(x,T,x,Tinf(x))
    xlabel('x')
    ylabel('T')
    legend(['t = ',num2str(tkeep(i))],'asymptote')
    title('Temperature Distribution for non-linear conductivity')
    drawnow
    end
end

%% Plot of Bottom Temperature
% We can plot up the bottom temperature vs. time.  We could have kept some
% of the intermediate temperature profiles as well, but not all of them!
% Because of the Neumann condition we have to use a very small
% discretization in time: way too much data!

figure(2)
plot(tkeep,Tbottomkeep)
xlabel('time')
ylabel('Temperature at x = 0')
grid on